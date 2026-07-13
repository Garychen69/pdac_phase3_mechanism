"""
preprocess_singlecell_cohorts.py
Load or simulate scRNA-seq data for each cohort.
Applies QC, normalization, and saves preprocessed AnnData objects.

Real data handling (updated for Phase 3 with downloaded data):
- GSE154778: CSV matrix (genes x cells), primary tumors only
- GSE202051: h5ad already QC'd and annotated (snRNA-seq)
- Peng_et_al: 10x CellRanger output extracted from RAW.tar
"""

import os
import sys
import re
import io
import tarfile
import tempfile
import numpy as np
import pandas as pd
import random
import yaml
import warnings
warnings.filterwarnings("ignore")

np.random.seed(1234)
random.seed(1234)

import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(BASE_DIR, "config")
DATA_DIR = os.path.join(BASE_DIR, "data")
sc.settings.verbosity = 1


def load_config():
    with open(os.path.join(CONFIG_DIR, "singlecell_cohorts.yml")) as f:
        return yaml.safe_load(f)

def load_gene_sets():
    with open(os.path.join(CONFIG_DIR, "gene_sets.yml")) as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------------
# Real-data detection (check for actual files BEFORE checking flag files)
# ---------------------------------------------------------------------------

def has_real_data_gse154778(raw_dir):
    """Check if GSE154778 CSV matrix is present."""
    csv_gz = os.path.join(raw_dir, "GSE154778_dgeMtx.csv.gz")
    return os.path.exists(csv_gz)


def has_real_data_gse202051(raw_dir):
    """Check if GSE202051 h5ad is present. Prefers the full 43-patient
    combined object; the old adata_010nuc_10x.h5ad is a single-sample
    subset (1 usable patient) kept only as a legacy fallback."""
    full = os.path.join(raw_dir, "GSE202051_full.h5ad")
    subset = os.path.join(raw_dir, "GSE202051_adata_010nuc_10x.h5ad")
    return os.path.exists(full) or os.path.exists(subset)


def has_real_data_peng(raw_dir):
    """Check if Peng RAW.tar exists and is > 900 MB (SI units)."""
    tar_path = os.path.join(raw_dir, "GSE155698_RAW.tar")
    if not os.path.exists(tar_path):
        return False
    return os.path.getsize(tar_path) > 900_000_000  # 900 MB in SI


def remove_download_failed_flag(cohort_name, raw_dir):
    """Remove DOWNLOAD_FAILED flag file after successful real-data load."""
    flag = os.path.join(raw_dir, f"{cohort_name}_DOWNLOAD_FAILED.txt")
    if os.path.exists(flag):
        try:
            os.remove(flag)
            print(f"  Removed DOWNLOAD_FAILED flag: {flag}")
        except Exception as e:
            print(f"  Could not remove flag: {e}")


# ---------------------------------------------------------------------------
# Real data loaders
# ---------------------------------------------------------------------------

def load_gse154778(raw_dir):
    """Load GSE154778 CSV matrix and return preprocessed AnnData."""
    csv_path = os.path.join(raw_dir, "GSE154778_dgeMtx.csv.gz")
    print(f"  Loading {csv_path} ...")

    # Rows = genes, columns = cell barcodes (with patient prefix like P03:1)
    df = pd.read_csv(csv_path, index_col=0)
    print(f"  Raw matrix: {df.shape[0]} genes × {df.shape[1]} cells")

    # Filter to primary tumor cells only (not MET)
    primary_mask = [bool(re.match(r'^P\d+:', c)) for c in df.columns]
    primary_cols = df.columns[primary_mask].tolist()
    df = df[primary_cols]
    print(f"  After primary-tumor filter: {len(primary_cols)} cells")

    # Build AnnData: cells × genes
    X_sparse = csr_matrix(df.values.T, dtype=np.float32)
    obs = pd.DataFrame(index=df.columns)
    var = pd.DataFrame(index=df.index)
    adata = ad.AnnData(X=X_sparse, obs=obs, var=var)

    # Patient ID from barcode prefix (e.g., 'P03' from 'P03:1')
    adata.obs["patient_id"] = [c.split(":")[0] for c in adata.obs.index]
    adata.obs["cohort"] = "GSE154778"
    adata.obs["is_simulated"] = False
    adata.var_names_make_unique()

    print(f"  Created AnnData: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def load_gse202051(raw_dir):
    """Load GSE202051 h5ad (already QC'd and annotated). Prefers the full
    43-patient combined object (totaldata-final-toshare) over the legacy
    single-sample subset (adata_010nuc_10x), which only yielded 1 usable
    patient and made this cohort untestable at the patient level."""
    full_path = os.path.join(raw_dir, "GSE202051_full.h5ad")
    subset_path = os.path.join(raw_dir, "GSE202051_adata_010nuc_10x.h5ad")
    h5ad_path = full_path if os.path.exists(full_path) else subset_path
    print(f"  Loading {h5ad_path} ...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"  Raw AnnData: {adata.n_obs} cells × {adata.n_vars} genes")

    # Set required metadata
    if "pid" in adata.obs.columns:
        adata.obs["patient_id"] = adata.obs["pid"].astype(str)
        print(f"  Patient count: {adata.obs['patient_id'].nunique()}")
    else:
        adata.obs["patient_id"] = "GSE202051_P01"

    adata.obs["cohort"] = "GSE202051"
    adata.obs["is_simulated"] = False

    # Check X and normalize if needed
    x_max = float(adata.X.max())
    print(f"  X max = {x_max:.4f}")
    if x_max > 100:
        print("  X appears to be raw counts → applying normalize_total + log1p")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif x_max > 10:
        print("  X appears partially normalized -> applying log1p only")
        sc.pp.log1p(adata)
    else:
        print("  X max < 10 -> likely already log-normalized, using as-is")

    # Store raw
    adata.raw = adata.copy()
    return adata


def load_peng_et_al(raw_dir):
    """Extract PDAC_TISSUE samples from RAW.tar and load 10x data."""
    tar_path = os.path.join(raw_dir, "GSE155698_RAW.tar")
    extract_dir = os.path.join(raw_dir, "extracted")
    os.makedirs(extract_dir, exist_ok=True)

    print(f"  Extracting PDAC_TISSUE samples from {tar_path} ...")
    with tarfile.open(tar_path, "r") as tar:
        for member in tar.getmembers():
            if "PDAC_TISSUE" in member.name and member.name.endswith(".tar.gz"):
                dest = os.path.join(extract_dir, member.name)
                if not os.path.exists(dest):
                    tar.extract(member, extract_dir)

    # Now extract each inner .tar.gz to its own directory
    sample_adatas = []
    inner_tars = sorted([f for f in os.listdir(extract_dir)
                         if f.endswith(".tar.gz") and "PDAC_TISSUE" in f])
    print(f"  Found {len(inner_tars)} PDAC_TISSUE samples")

    for fname in inner_tars:
        sample_name = fname.replace(".tar.gz", "")
        sample_dir = os.path.join(extract_dir, sample_name)
        os.makedirs(sample_dir, exist_ok=True)

        inner_tar_path = os.path.join(extract_dir, fname)
        # Extract if not done yet
        if not any(os.listdir(sample_dir)):
            try:
                with tarfile.open(inner_tar_path, "r:gz") as inner_tar:
                    inner_tar.extractall(sample_dir)
            except Exception as e:
                print(f"  Could not extract {fname}: {e}")
                continue

        # Find the actual directory containing matrix.mtx.gz (may be nested)
        mtx_dir = None
        for root, dirs, files in os.walk(sample_dir):
            if any(f.startswith("matrix") for f in files):
                mtx_dir = root
                break
        if mtx_dir is None:
            mtx_dir = sample_dir  # fallback

        try:
            # Try MTX format first
            if mtx_dir and any(f.startswith("matrix") for f in os.listdir(mtx_dir)):
                sample_adata = sc.read_10x_mtx(
                    mtx_dir, var_names="gene_symbols", cache=True
                )
            else:
                # Try h5 format
                h5_files = []
                for root, dirs, files in os.walk(sample_dir):
                    for f in files:
                        if f.endswith(".h5"):
                            h5_files.append(os.path.join(root, f))
                if not h5_files:
                    raise FileNotFoundError(f"No matrix.mtx.gz or .h5 files found in {sample_name}")
                sample_adata = sc.read_10x_h5(h5_files[0])
            sample_adata.obs["patient_id"] = sample_name
            sample_adata.obs["cohort"] = "Peng_et_al"
            sample_adata.obs["is_simulated"] = False
            sample_adata.var_names_make_unique()
            sample_adatas.append(sample_adata)
            print(f"    Loaded {sample_name}: {sample_adata.n_obs} cells")
        except Exception as e:
            print(f"    Failed to load {sample_name}: {e}")

    if not sample_adatas:
        raise RuntimeError("No Peng_et_al samples loaded from 10x matrices")

    # Concatenate
    adata = ad.concat(sample_adatas, join="outer", fill_value=0, label="sample")
    adata.var_names_make_unique()
    print(f"  Concatenated: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


# ---------------------------------------------------------------------------
# Simulated data fallback
# ---------------------------------------------------------------------------

def collect_all_genes(gene_sets):
    """Collect all unique genes from gene sets and markers."""
    genes = set()
    for key, val in gene_sets.items():
        if isinstance(val, list):
            genes.update(val)
        elif isinstance(val, dict):
            for subval in val.values():
                if isinstance(subval, list):
                    genes.update(subval)
    return sorted(genes)


def simulate_anndata(cohort_name, sim_params, gene_sets, patient_seed_offset=0):
    """Generate synthetic scRNA-seq AnnData with realistic distributions."""
    np.random.seed(1234 + patient_seed_offset)

    n_cells = sim_params.get("n_cells", 2000)
    n_genes_total = sim_params.get("n_genes", 500)
    n_patients = sim_params.get("n_patients", 10)
    marker_boost = sim_params.get("marker_boost_factor", 3.0)
    cell_type_fractions = sim_params["cell_type_fractions"]

    signature_genes = collect_all_genes(gene_sets)
    background_genes = [f"GENE_{i:04d}" for i in range(n_genes_total)]
    all_genes = list(set(signature_genes + background_genes))
    all_genes.sort()
    n_genes = len(all_genes)
    gene_to_idx = {g: i for i, g in enumerate(all_genes)}

    cells_per_patient = n_cells // n_patients
    patient_ids = []
    for p in range(n_patients):
        n = cells_per_patient if p < n_patients - 1 else n_cells - cells_per_patient * (n_patients - 1)
        patient_ids.extend([f"{cohort_name}_P{p+1:02d}"] * n)
    patient_ids = np.array(patient_ids)

    cell_types_list = list(cell_type_fractions.keys())
    cell_type_probs = np.array([cell_type_fractions[ct] for ct in cell_types_list])
    cell_type_probs = cell_type_probs / cell_type_probs.sum()

    np.random.seed(1234 + patient_seed_offset)
    cell_type_assignments = np.random.choice(cell_types_list, size=n_cells, p=cell_type_probs)

    n_aggressive = n_patients // 2
    aggressive_patients = [f"{cohort_name}_P{p+1:02d}" for p in range(n_aggressive)]
    is_aggressive = np.array([pid in aggressive_patients for pid in patient_ids])

    np.random.seed(1234 + patient_seed_offset + 1)
    X_base = np.random.lognormal(mean=0.5, sigma=1.2, size=(n_cells, n_genes)).astype(np.float32)

    markers = gene_sets.get("cell_type_markers", {})
    for ct, marker_list in markers.items():
        ct_mask = (cell_type_assignments == ct)
        if not ct_mask.any():
            continue
        for gene in marker_list:
            if gene in gene_to_idx:
                X_base[ct_mask, gene_to_idx[gene]] *= marker_boost

    mal_mask = (cell_type_assignments == "malignant_epithelial")
    agg_mal_mask = mal_mask & is_aggressive

    hypoxia_genes = gene_sets.get("hypoxia", [])
    lipid_syn_genes = gene_sets.get("lipid_synthesis_srebp", [])
    desat_genes = gene_sets.get("desaturation_elongation", [])
    faox_genes = gene_sets.get("fatty_acid_uptake_oxidation", [])
    emt_genes_list = gene_sets.get("emt", [])
    acinar_genes = gene_sets.get("acinar_identity", [])

    for gene in hypoxia_genes:
        if gene in gene_to_idx:
            X_base[agg_mal_mask, gene_to_idx[gene]] *= 2.5
    for gene in lipid_syn_genes + desat_genes:
        if gene in gene_to_idx:
            X_base[agg_mal_mask, gene_to_idx[gene]] *= 2.0
    for gene in faox_genes:
        if gene in gene_to_idx:
            X_base[agg_mal_mask, gene_to_idx[gene]] *= 0.5
    for gene in emt_genes_list:
        if gene in gene_to_idx:
            X_base[agg_mal_mask, gene_to_idx[gene]] *= 2.0
    for gene in acinar_genes:
        if gene in gene_to_idx:
            X_base[agg_mal_mask, gene_to_idx[gene]] *= 0.4

    ref_mal_mask = mal_mask & ~is_aggressive
    for gene in acinar_genes:
        if gene in gene_to_idx:
            X_base[ref_mal_mask, gene_to_idx[gene]] *= 2.0
    for gene in hypoxia_genes:
        if gene in gene_to_idx:
            X_base[ref_mal_mask, gene_to_idx[gene]] *= 0.6

    caf_mask = (cell_type_assignments == "caf_fibroblast")
    agg_caf = caf_mask & is_aggressive
    caf_genes = gene_sets.get("caf", [])
    for gene in caf_genes + emt_genes_list:
        if gene in gene_to_idx:
            X_base[agg_caf, gene_to_idx[gene]] *= 1.8

    np.random.seed(1234 + patient_seed_offset + 2)
    noise = np.random.normal(0, 0.3, X_base.shape).astype(np.float32)
    X_base = np.clip(X_base + noise, 0, None)
    X_counts = np.round(X_base).astype(np.float32)

    obs = pd.DataFrame({
        "cell_id": [f"{cohort_name}_cell_{i:05d}" for i in range(n_cells)],
        "patient_id": patient_ids,
        "cohort": cohort_name,
        "true_cell_type": cell_type_assignments,
        "is_aggressive_patient": is_aggressive,
        "is_simulated": True,
    })
    obs.index = obs["cell_id"].values

    var = pd.DataFrame({"gene_name": all_genes}, index=all_genes)
    var["is_marker_gene"] = var.index.isin(signature_genes)

    adata = ad.AnnData(X=csr_matrix(X_counts), obs=obs, var=var)
    adata.var_names_make_unique()

    print(f"  Simulated AnnData: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


# ---------------------------------------------------------------------------
# QC + normalization
# ---------------------------------------------------------------------------

def preprocess_real_adata(adata, cohort_name):
    """Apply QC and normalization to real data AnnData."""
    print(f"  Preprocessing {cohort_name} (real data)...")

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Mitochondrial filter
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    if adata.var["mt"].sum() > 0:
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                                   log1p=False, inplace=True)
        before = adata.n_obs
        adata = adata[adata.obs["pct_counts_mt"] < 25].copy()
        print(f"  Mito filter: {before} -> {adata.n_obs} cells")
    else:
        print("  No MT- genes found; skipping mito filter")

    # Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Store raw
    adata.raw = adata.copy()

    print(f"  After QC/norm: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def preprocess_simulated_adata(adata, cohort_name):
    """Preprocess simulated data (already clean)."""
    print(f"  Preprocessing {cohort_name} (simulated)...")
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    print(f"  After QC/norm: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


# ---------------------------------------------------------------------------
# Gene coverage check
# ---------------------------------------------------------------------------

def check_gene_coverage(adata, gene_sets, cohort_name):
    """Check and report how many signature genes are present in the dataset."""
    print(f"\n  Gene coverage report for {cohort_name}:")
    sig_names = ["hypoxia", "acinar_identity", "lipid_synthesis_srebp",
                 "desaturation_elongation", "fatty_acid_uptake_oxidation", "caf", "emt"]
    coverage = {}
    for sig in sig_names:
        genes = gene_sets.get(sig, [])
        found = [g for g in genes if g in adata.var_names]
        pct = 100 * len(found) / len(genes) if genes else 0
        coverage[sig] = {"found": len(found), "total": len(genes), "pct": pct}
        status = "OK" if pct >= 50 else "LOW"
        print(f"    {sig}: {len(found)}/{len(genes)} ({pct:.0f}%) [{status}]")
    return coverage


def check_dissociation_stress(adata, cohort_name):
    """Check for dissociation-stress genes and warn if highly expressed."""
    stress_genes = ["FOS", "JUN", "HSPA1A", "HSPA1B", "DNAJB1", "EGR1", "ATF3"]
    found = [g for g in stress_genes if g in adata.var_names]
    if not found:
        print(f"  Dissociation stress genes: none found in dataset")
        return False
    # Compute mean expression of stress genes. scipy sparse matrices support
    # .mean() natively without densifying first -- densifying the full matrix
    # (previously done via .toarray() below) needs n_cells x n_genes x 4 bytes,
    # which is fine for small cohorts but allocates double-digit GB and crashes
    # on a 225k-cell x 22k-gene cohort like the full GSE202051 object.
    gene_idx = [list(adata.var_names).index(g) for g in found]
    stress_expr = adata.X[:, gene_idx].mean()
    overall_mean = adata.X.mean()
    ratio = stress_expr / (overall_mean + 1e-9)
    if ratio > 1.5:
        print(f"  WARNING: Dissociation stress genes elevated (ratio={ratio:.2f}). "
              f"Hypoxia scores may be confounded by dissociation artifact.")
        return True
    else:
        print(f"  Dissociation stress check OK (ratio={ratio:.2f})")
        return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=== Preprocessing Single-Cell Cohorts (Real Data Mode) ===\n")

    config = load_config()
    gene_sets = load_gene_sets()
    sim_params = config["simulation_params"]
    cohorts = config["cohorts"]

    for i, cohort in enumerate(cohorts):
        name = cohort["name"]
        raw_dir = os.path.join(BASE_DIR, cohort["raw_dir"])
        processed_file = os.path.join(BASE_DIR, cohort["processed_file"])
        os.makedirs(os.path.dirname(processed_file), exist_ok=True)

        print(f"--- Cohort: {name} ---")

        adata = None
        is_real = False

        # ---- GSE154778 ----
        if name == "GSE154778":
            if has_real_data_gse154778(raw_dir):
                print(f"  Real CSV data found.")
                try:
                    adata = load_gse154778(raw_dir)
                    adata = preprocess_real_adata(adata, name)
                    adata.uns["simulated"] = False
                    adata.uns["data_source"] = "REAL"
                    is_real = True
                    remove_download_failed_flag(name, raw_dir)
                except Exception as e:
                    import traceback
                    print(f"  ERROR loading GSE154778: {e}")
                    traceback.print_exc()
                    adata = None
            else:
                print("  No real CSV data found.")

        # ---- GSE202051 ----
        elif name == "GSE202051":
            if has_real_data_gse202051(raw_dir):
                print(f"  Real h5ad data found.")
                try:
                    adata = load_gse202051(raw_dir)
                    # QC already done; just filter genes
                    sc.pp.filter_genes(adata, min_cells=1)
                    adata.uns["simulated"] = False
                    adata.uns["data_source"] = "REAL"
                    is_real = True
                    remove_download_failed_flag(name, raw_dir)
                except Exception as e:
                    import traceback
                    print(f"  ERROR loading GSE202051: {e}")
                    traceback.print_exc()
                    adata = None
            else:
                print("  No real h5ad found.")

        # ---- Peng_et_al ----
        elif name == "Peng_et_al":
            if has_real_data_peng(raw_dir):
                print(f"  Real RAW.tar data found (>900 MB).")
                try:
                    adata = load_peng_et_al(raw_dir)
                    adata = preprocess_real_adata(adata, name)
                    adata.uns["simulated"] = False
                    adata.uns["data_source"] = "REAL"
                    is_real = True
                    remove_download_failed_flag(name, raw_dir)
                except Exception as e:
                    import traceback
                    print(f"  ERROR loading Peng_et_al: {e}")
                    traceback.print_exc()
                    adata = None
            else:
                print("  RAW.tar not found or < 900 MB.")

        # ---- Fallback: simulate ----
        if adata is None:
            print(f"  Falling back to SIMULATED DATA.")
            adata = simulate_anndata(name, sim_params, gene_sets, patient_seed_offset=i * 100)
            adata = preprocess_simulated_adata(adata, name)
            adata.uns["simulated"] = True
            adata.uns["data_source"] = "SIMULATED"
            adata.uns["data_note"] = "SIMULATED DATA — awaiting real dataset download"

        adata.uns["cohort"] = name

        # Gene coverage + dissociation stress check
        if is_real:
            coverage = check_gene_coverage(adata, gene_sets, name)
            adata.uns["gene_coverage"] = coverage
            stress_flag = check_dissociation_stress(adata, name)
            adata.uns["dissociation_stress_elevated"] = bool(stress_flag)

        adata.write_h5ad(processed_file)
        print(f"  Saved: {processed_file}\n")

    print("=== Preprocessing Complete ===")


if __name__ == "__main__":
    main()
