"""
annotate_cell_types.py
Run PCA, neighbors, Leiden clustering, cell type scoring, and annotation.
Saves UMAP, cell type labels, and figures.

Special handling:
- GSE202051: already has UMAP, Leiden, and cell type score columns → use directly
- GSE154778, Peng_et_al: full pipeline (PCA → neighbors → UMAP → Leiden → scoring)
"""

import os
import sys
import numpy as np
import pandas as pd
import random
import yaml
import warnings
warnings.filterwarnings("ignore")

np.random.seed(1234)
random.seed(1234)

import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(BASE_DIR, "config")
PROCESSED_SC_DIR = os.path.join(BASE_DIR, "data", "processed", "singlecell")
FIGURES_DIR = os.path.join(BASE_DIR, "results", "figures")
sc.settings.verbosity = 1

# Canonical markers for de-novo annotation
CANONICAL_MARKERS = {
    "malignant_epithelial": ["KRT19", "KRT8", "KRT18", "EPCAM", "SOX9"],
    "acinar_normal": ["PRSS1", "PRSS2", "CPA1", "CPA2", "CELA3A", "AMY2A"],
    "ductal_normal": ["KRT19", "CFTR", "SLC4A4"],
    "caf_fibroblast": ["COL1A1", "COL1A2", "DCN", "PDGFRB", "ACTA2"],
    "endothelial": ["PECAM1", "VWF", "CDH5"],
    "myeloid": ["CD68", "CSF1R", "LYZ"],
    "tcell_nk": ["CD3D", "CD3E", "NKG7", "GNLY"],
    "bcell_plasma": ["MS4A1", "CD79A", "MZB1"],
}

# Mapping from GSE202051 score columns to standard vocabulary (legacy path,
# used only if the richer 'new_celltypes' column below is unavailable)
GSE202051_TYPE_MAP = {
    "MALIGNANT CELLS": "malignant_epithelial",
    "ACINAR": "acinar_normal",
    "DUCTAL": "ductal_normal",
    "FIBROBLASTS": "caf_fibroblast",
    "CAF": "caf_fibroblast",
    "ENDOTHELIAL": "endothelial",
    "IMMUNE": "myeloid",
}


def map_gse202051_new_celltype(label):
    """Map the full GSE202051 object's fine-grained 'new_celltypes' labels
    (e.g. 'Epithelial-Malignant', 'Fibroblast-myCAF', 'Immune-CD8pos_Tcells')
    to this pipeline's standard cell-type vocabulary (CELL_TYPES_ORDERED
    below). Cell types with no clean standard-vocabulary equivalent
    (pericytes, vascular smooth muscle, endocrine, Schwann, adipocyte)
    map to 'unknown', consistent with how unmapped types are handled
    elsewhere in this pipeline."""
    if not isinstance(label, str):
        return "unknown"
    if label.startswith("Epithelial-Malignant"):
        return "malignant_epithelial"
    if label.startswith(("Epithelial-Ductal", "Epithelial-Atypical_Ductal", "Epithelial-ADM")):
        return "ductal_normal"
    if label.startswith("Epithelial-Acinar"):
        return "acinar_normal"
    if label.startswith(("Fibroblast-CAF", "Fibroblast-myCAF", "Fibroblast-iCAF")):
        return "caf_fibroblast"
    if label.startswith("Fibroblast-"):
        return "unknown"  # Pericyte, vSMC — not in standard vocabulary
    if label.startswith("Endothelial-"):
        return "endothelial"
    if label == "Immune-Macrophage" or label.startswith(("Immune-cDC", "Immune-pDC", "Immune-DC_activated",
                                                           "Immune-Mast", "Immune-Neutrophil")):
        return "myeloid"
    if label in ("Immune-B", "Immune-Plasma"):
        return "bcell_plasma"
    if label.startswith("Immune-"):
        return "tcell_nk"  # CD4/CD8/NK/Treg/dysfunctional-T subsets
    return "unknown"  # Endocrine-*, Schwann, Adipocyte, etc.

CELL_TYPES_ORDERED = [
    "malignant_epithelial", "caf_fibroblast", "endothelial",
    "myeloid", "tcell_nk", "bcell_plasma", "acinar_normal", "ductal_normal"
]


def load_configs():
    with open(os.path.join(CONFIG_DIR, "singlecell_cohorts.yml")) as f:
        cohort_cfg = yaml.safe_load(f)
    with open(os.path.join(CONFIG_DIR, "gene_sets.yml")) as f:
        gene_sets = yaml.safe_load(f)
    return cohort_cfg, gene_sets


# ---------------------------------------------------------------------------
# Full annotation pipeline (for GSE154778 and Peng_et_al)
# ---------------------------------------------------------------------------

def run_full_annotation_pipeline(adata, name):
    """Run PCA → neighbors → UMAP → Leiden → score → assign cell type."""
    np.random.seed(1234)

    # Highly variable genes
    n_pcs = min(50, adata.n_vars - 1, adata.n_obs - 1)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    if adata.var["highly_variable"].sum() < 50:
        adata.var["highly_variable"] = True

    sc.tl.pca(adata, svd_solver="arpack", n_comps=n_pcs, random_state=1234)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(n_pcs, 30), random_state=1234)
    sc.tl.umap(adata, random_state=1234)
    sc.tl.leiden(adata, resolution=0.5, random_state=1234)

    # Score each canonical cell type
    score_keys = []
    for ct in CELL_TYPES_ORDERED:
        markers = CANONICAL_MARKERS.get(ct, [])
        available = [g for g in markers if g in adata.var_names]
        if len(available) >= 1:
            key = f"score_{ct}"
            sc.tl.score_genes(adata, available, score_name=key, random_state=1234)
            score_keys.append((ct, key))

    # Assign cell type to each cluster based on highest average marker score
    cluster_cell_type = {}
    for cluster in adata.obs["leiden"].unique():
        cluster_mask = adata.obs["leiden"] == cluster
        best_ct = "unknown"
        best_score = -np.inf
        for ct, key in score_keys:
            if key in adata.obs.columns:
                avg = adata.obs.loc[cluster_mask, key].mean()
                if avg > best_score:
                    best_score = avg
                    best_ct = ct
        cluster_cell_type[cluster] = best_ct

    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_cell_type)
    return adata, score_keys


def annotate_with_simulation_labels(adata):
    """For simulated data: use true_cell_type directly."""
    adata.obs["cell_type"] = adata.obs["true_cell_type"].values
    cluster_cell_type = {}
    for cluster in adata.obs["leiden"].unique():
        cluster_mask = adata.obs["leiden"] == cluster
        ct_counts = adata.obs.loc[cluster_mask, "true_cell_type"].value_counts()
        cluster_cell_type[cluster] = ct_counts.index[0] if len(ct_counts) > 0 else "unknown"
    adata.obs["leiden_cell_type"] = adata.obs["leiden"].map(cluster_cell_type)
    return adata


# ---------------------------------------------------------------------------
# GSE202051 special: use existing annotations
# ---------------------------------------------------------------------------

def annotate_gse202051(adata, name):
    """Assign cell types using existing annotations; reuse existing UMAP/Leiden."""
    type_candidates = {}
    if "new_celltypes" in adata.obs.columns:
        # Full 43-patient object: use the real fine-grained annotation column
        # directly rather than the old argmax-over-score-columns heuristic.
        n_fine = adata.obs["new_celltypes"].nunique()
        print(f"  GSE202051: mapping {n_fine} real 'new_celltypes' labels to standard vocabulary.")
        adata.obs["cell_type"] = adata.obs["new_celltypes"].map(map_gse202051_new_celltype)
        unmapped = int((adata.obs["cell_type"] == "unknown").sum())
        print(f"    {unmapped}/{adata.n_obs} cells map to 'unknown' (pericyte/endocrine/Schwann/etc.)")
    else:
        print(f"  GSE202051: using existing UMAP, Leiden, and cell-type scores.")

        # Build priority list (use only columns that exist)
        type_candidates = {}
        for col, std_type in GSE202051_TYPE_MAP.items():
            if col in adata.obs.columns:
                if std_type not in type_candidates:
                    type_candidates[std_type] = col
                # Prefer 'MALIGNANT CELLS' over 'CAF' for malignant, etc.

        # For each cell: pick standard type with highest obs score
        score_df = pd.DataFrame(index=adata.obs.index)
        for std_type, col in type_candidates.items():
            score_df[std_type] = adata.obs[col].values.astype(float)

        if score_df.shape[1] > 0:
            adata.obs["cell_type"] = score_df.idxmax(axis=1).values
        else:
            adata.obs["cell_type"] = "unknown"

    # Ensure leiden column exists
    if "leiden" not in adata.obs.columns:
        print("  WARNING: leiden not found; assigning placeholder")
        adata.obs["leiden"] = "0"

    # Ensure UMAP exists
    if "X_umap" not in adata.obsm:
        print("  WARNING: X_umap not in obsm; running UMAP now")
        n_pcs = min(30, adata.n_vars - 1, adata.n_obs - 1)
        if "X_pca" not in adata.obsm:
            sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            sc.tl.pca(adata, n_comps=n_pcs, random_state=1234)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(n_pcs, 30), random_state=1234)
        sc.tl.umap(adata, random_state=1234)

    score_keys = [(std, f"gse202051_{std}") for std in type_candidates]
    return adata, score_keys


# ---------------------------------------------------------------------------
# UMAP figure
# ---------------------------------------------------------------------------

def save_umap_figure(adata, name, score_keys):
    """Save UMAP colored by cell type and Leiden clusters."""
    os.makedirs(FIGURES_DIR, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    umap_coords = adata.obsm["X_umap"]
    cell_types_unique = adata.obs["cell_type"].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, max(len(cell_types_unique), 1)))
    ct_color_map = dict(zip(cell_types_unique, colors))

    ax = axes[0]
    for ct in cell_types_unique:
        mask = adata.obs["cell_type"] == ct
        ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                   c=[ct_color_map[ct]], label=ct, s=2, alpha=0.6)
    ax.set_title(f"{name}: Cell Type Annotation")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(markerscale=4, bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=7)

    ax2 = axes[1]
    try:
        leiden_vals = adata.obs["leiden"].astype(int).values
    except (ValueError, TypeError):
        leiden_vals = pd.Categorical(adata.obs["leiden"]).codes
    scatter = ax2.scatter(umap_coords[:, 0], umap_coords[:, 1],
                          c=leiden_vals, cmap="tab20", s=2, alpha=0.6)
    ax2.set_title(f"{name}: Leiden Clusters")
    ax2.set_xlabel("UMAP1")
    ax2.set_ylabel("UMAP2")
    plt.colorbar(scatter, ax=ax2, label="Cluster")

    if adata.uns.get("simulated", False):
        fig.suptitle(f"[SIMULATED DATA] {name}", color="red", fontsize=12)
    elif name == "GSE202051":
        fig.suptitle(f"{name} (snRNA-seq, pre-annotated)", fontsize=12)

    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, f"FigureS1_{name}_celltype_annotation_umap.pdf")
    plt.savefig(fig_path, bbox_inches="tight")
    plt.close()
    print(f"  UMAP figure saved: {fig_path}")


# ---------------------------------------------------------------------------
# Per-cohort annotation dispatch
# ---------------------------------------------------------------------------

def annotate_cohort(cohort_cfg, gene_sets):
    name = cohort_cfg["name"]
    processed_file = os.path.join(BASE_DIR, cohort_cfg["processed_file"])
    annotations_file = os.path.join(BASE_DIR, cohort_cfg["annotations_file"])

    if not os.path.exists(processed_file):
        print(f"  Preprocessed file not found: {processed_file}. Skipping.")
        return

    print(f"  Loading {name}...")
    adata = sc.read_h5ad(processed_file)
    is_simulated = adata.uns.get("simulated", True)
    is_real = not is_simulated

    # Choose annotation strategy
    if name == "GSE202051" and is_real:
        # Rich pre-existing annotations
        adata, score_keys = annotate_gse202051(adata, name)
    elif is_simulated and "true_cell_type" in adata.obs.columns:
        # Simulated: need PCA/UMAP/Leiden first, then use ground truth labels
        np.random.seed(1234)
        n_pcs = min(50, adata.n_vars - 1, adata.n_obs - 1)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        if adata.var["highly_variable"].sum() < 50:
            adata.var["highly_variable"] = True
        sc.tl.pca(adata, svd_solver="arpack", n_comps=n_pcs, random_state=1234)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(n_pcs, 30), random_state=1234)
        sc.tl.umap(adata, random_state=1234)
        sc.tl.leiden(adata, resolution=0.5, random_state=1234)

        score_keys = []
        for ct in CELL_TYPES_ORDERED:
            markers = CANONICAL_MARKERS.get(ct, [])
            available = [g for g in markers if g in adata.var_names]
            if len(available) >= 1:
                key = f"score_{ct}"
                sc.tl.score_genes(adata, available, score_name=key, random_state=1234)
                score_keys.append((ct, key))

        print(f"  Simulated data: using true_cell_type labels directly.")
        adata = annotate_with_simulation_labels(adata)
    else:
        # Real data without pre-existing annotations: full pipeline
        adata, score_keys = run_full_annotation_pipeline(adata, name)

    # Save UMAP figure
    save_umap_figure(adata, name, score_keys)

    # Build annotation TSV
    ann_df = pd.DataFrame({
        "cell_id": adata.obs.index,
        "patient_id": adata.obs.get("patient_id", pd.Series(
            adata.obs.index, index=adata.obs.index)),
        "cohort": name,
        "leiden_cluster": adata.obs.get("leiden", pd.Series(
            ["0"] * adata.n_obs, index=adata.obs.index)).values,
        "cell_type": adata.obs["cell_type"].values,
        "UMAP1": adata.obsm["X_umap"][:, 0],
        "UMAP2": adata.obsm["X_umap"][:, 1],
        "is_simulated": adata.uns.get("simulated", False),
    })
    for ct, key in score_keys:
        if key in adata.obs.columns:
            ann_df[f"marker_score_{ct}"] = adata.obs[key].values

    os.makedirs(os.path.dirname(annotations_file), exist_ok=True)
    ann_df.to_csv(annotations_file, sep="\t", index=False)
    print(f"  Annotations saved: {annotations_file}")

    # Save updated h5ad with annotations
    adata.write_h5ad(processed_file)
    return adata


def main():
    print("=== Cell Type Annotation ===\n")
    cohort_cfg, gene_sets = load_configs()

    for cohort in cohort_cfg["cohorts"]:
        print(f"--- Annotating: {cohort['name']} ---")
        try:
            annotate_cohort(cohort, gene_sets)
        except Exception as e:
            import traceback
            print(f"  ERROR in {cohort['name']}: {e}")
            traceback.print_exc()
        print()

    print("=== Annotation Complete ===")


if __name__ == "__main__":
    main()
