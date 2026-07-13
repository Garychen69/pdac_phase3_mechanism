"""
sensitivity_treatment_status_GSE202051.py

Sensitivity analysis: re-run the GSE202051 hypoxia/acinar co-occurrence,
lipid cell-of-origin, EMT cell-of-origin, and CAF subtype comparisons
restricted to the 18 untreated patients only, to check whether pooling
untreated (18) and neoadjuvant-treated (25) patients in the main analysis
(Limitation 8 in PHASE3_MECHANISM_REPORT.md) materially changes any
conclusion. Reuses the exact same pseudobulk/testing logic as the
production scripts (imported, not reimplemented) so results are directly
comparable to the "full" (pooled) numbers already in the report.
"""

import os
import sys
import numpy as np
import pandas as pd
import yaml
import warnings
warnings.filterwarnings("ignore")

np.random.seed(1234)

from scipy import stats
from statsmodels.stats.multitest import multipletests
import anndata as ad
import scanpy as sc
sc.settings.verbosity = 0

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from analyze_lipid_cell_of_origin import (
    LIPID_SCORES, CELL_TYPES_TO_ANALYZE, MIN_CELLS_PER_PATIENT, MIN_PATIENTS_PER_ARM,
    pseudobulk_by_patient,
)
from analyze_caf_emt_cell_of_origin import get_patient_state, pseudobulk_median, score_caf_subtypes

TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
SCORES_FILE = os.path.join(BASE_DIR, "data", "processed", "singlecell", "GSE202051_cell_scores.tsv")
# Read from a local (non-OneDrive-synced) copy, not the project folder directly:
# scanpy's many small HDF5 reads against this ~5GB object are catastrophically
# slow through OneDrive's cloud-sync filter driver (same issue documented in
# config/singlecell_cohorts.yml from the 2026-07-12 GSE202051 upgrade).
PROCESSED_FILE = "C:/Users/erica/AppData/Local/Temp/pdac_phase3_local/GSE202051_preprocessed.h5ad"
GENE_SETS_FILE = os.path.join(BASE_DIR, "config", "gene_sets.yml")


def get_untreated_patients():
    """Patient IDs with treatment_status == 'Untreated' in the raw object
    (18 of 43 patients; the other 25 received some form of neoadjuvant
    chemo/radiotherapy per the treatment_status field)."""
    a = ad.read_h5ad(PROCESSED_FILE, backed="r")
    pt = a.obs.groupby("patient_id", observed=True)["treatment_status"].first()
    untreated = set(pt[pt == "Untreated"].index.astype(str))
    return untreated


def hypoxia_acinar_correlation(df, label, results):
    mal = df[df["cell_type"] == "malignant_epithelial"]
    r, p = stats.pearsonr(mal["hypoxia_score"], mal["acinar_identity_score"])
    frac_quadrant = ((mal["hypoxia_score"] >= mal["hypoxia_score"].median()) &
                      (mal["acinar_identity_score"] < mal["acinar_identity_score"].median())).mean()
    results.append({
        "analysis": "hypoxia_acinar_correlation", "subset": label,
        "n_malignant_cells": len(mal), "n_patients": mal["patient_id"].nunique(),
        "pearson_r": round(r, 4), "p_value": p, "frac_hypoxia_high_acinar_low": round(frac_quadrant, 4),
    })


def lipid_cell_of_origin(df, patient_state, label, results):
    for cell_type in CELL_TYPES_TO_ANALYZE:
        ct_df = df[df["cell_type"] == cell_type]
        pb = pseudobulk_by_patient(ct_df, [c + "_score" for c in LIPID_SCORES])
        pb["state"] = pb.index.map(patient_state)
        agg = pb[pb["state"] == "aggressive"]
        ref = pb[pb["state"] == "reference"]
        if len(agg) < MIN_PATIENTS_PER_ARM or len(ref) < MIN_PATIENTS_PER_ARM:
            continue
        pvals, rows = [], []
        for score_col, expected_dir in LIPID_SCORES.items():
            col = score_col + "_score"
            a_vals, r_vals = agg[col].values, ref[col].values
            stat, pval = stats.ranksums(a_vals, r_vals)
            observed_dir = "up" if np.median(a_vals) > np.median(r_vals) else "down"
            pvals.append(pval)
            rows.append((score_col, expected_dir, observed_dir, np.median(a_vals) - np.median(r_vals)))
        _, p_adj, _, _ = multipletests(pvals, method="fdr_bh")
        for j, (score_col, exp_dir, obs_dir, median_diff) in enumerate(rows):
            results.append({
                "analysis": "lipid_cell_of_origin", "subset": label, "cell_type": cell_type,
                "score": score_col, "n_aggressive": len(agg), "n_reference": len(ref),
                "expected_direction": exp_dir, "observed_direction": obs_dir,
                "median_diff_agg_minus_ref": round(median_diff, 4),
                "p_raw": round(pvals[j], 4), "p_adj_BH": round(p_adj[j], 4),
                "significant_FDR05": p_adj[j] < 0.05, "direction_reproduced": obs_dir == exp_dir,
            })


def emt_cell_of_origin(df, patient_state, label, results):
    for ct in ["malignant_epithelial", "caf_fibroblast", "myeloid"]:
        ct_df = df[df["cell_type"] == ct]
        pb, cell_counts = pseudobulk_median(ct_df, "emt_score")
        state_map = pb.index.map(patient_state)
        agg_emt, ref_emt = pb[state_map == "aggressive"].values, pb[state_map == "reference"].values
        n_agg, n_ref = int((state_map == "aggressive").sum()), int((state_map == "reference").sum())
        if n_agg < MIN_PATIENTS_PER_ARM or n_ref < MIN_PATIENTS_PER_ARM:
            continue
        stat, pval = stats.ranksums(agg_emt, ref_emt)
        median_diff = np.median(agg_emt) - np.median(ref_emt)
        results.append({
            "analysis": "emt_cell_of_origin", "subset": label, "cell_type": ct,
            "n_aggressive": n_agg, "n_reference": n_ref,
            "median_diff_agg_minus_ref": round(median_diff, 4), "p_raw": round(pval, 4),
            "observed_direction": "up" if median_diff > 0 else "down",
        })


def caf_subtype_proportions(scores_df, patient_state, untreated_patients, label, results):
    """Re-score CAF subtypes on the AnnData, matching the production script's
    order of operations exactly: score_genes is run on the full population
    of cells present (all cell types) BEFORE restricting to CAF cells, since
    scanpy's score_genes samples control genes from the whole dataset's
    expression distribution — scoring after subsetting to CAF-only cells
    changes that background distribution and produces non-comparable
    scores/assignments (verified: doing it in the wrong order reproduced
    implausible ~50% apCAF fractions that did not match the published
    Section 5 numbers). For 'untreated_only' the patient restriction is
    applied first (across all cell types), then scoring, then CAF subsetting
    — i.e. as if the cohort had only ever had these 18 patients."""
    with open(GENE_SETS_FILE) as f:
        gene_sets = yaml.safe_load(f)

    adata = sc.read_h5ad(PROCESSED_FILE)
    if label == "untreated_only":
        adata = adata[adata.obs["patient_id"].astype(str).isin(untreated_patients)].copy()

    score_keys, adata = score_caf_subtypes(adata, gene_sets)
    caf_adata = adata[adata.obs["cell_type"] == "caf_fibroblast"].copy()
    del adata
    subtype_df = pd.DataFrame(
        {st: caf_adata.obs[key].values for st, key in score_keys.items()},
        index=caf_adata.obs.index,
    )
    caf_adata.obs["caf_subtype"] = subtype_df.idxmax(axis=1).values
    caf_adata.obs["patient_state"] = caf_adata.obs["patient_id"].astype(str).map(patient_state).fillna("intermediate")

    for state in ["aggressive", "reference"]:
        sub = caf_adata.obs[caf_adata.obs["patient_state"] == state]
        if len(sub) == 0:
            continue
        counts = sub.groupby("patient_id", observed=True).size()
        valid = counts[counts >= MIN_CELLS_PER_PATIENT].index
        if len(valid) == 0:
            continue
        props = (sub[sub["patient_id"].isin(valid)]
                 .groupby("patient_id", observed=True)["caf_subtype"]
                 .value_counts(normalize=True).unstack(fill_value=0.0))
        mean_props = props.mean(axis=0)
        for st in score_keys:
            results.append({
                "analysis": "caf_subtype_proportion", "subset": label, "patient_group": state,
                "caf_subtype": st, "n_patients": len(valid), "n_caf_cells": len(sub),
                "mean_fraction": round(float(mean_props.get(st, 0.0)), 4),
            })


def main():
    print("=== GSE202051 Treatment-Status Sensitivity Analysis ===\n")
    untreated_patients = get_untreated_patients()
    print(f"Untreated patients: {len(untreated_patients)} of 43\n")

    df_full = pd.read_csv(SCORES_FILE, sep="\t")
    df_full["patient_id"] = df_full["patient_id"].astype(str)
    df_untreated = df_full[df_full["patient_id"].isin(untreated_patients)].copy()

    print(f"Full (pooled):     {df_full['patient_id'].nunique()} patients, {len(df_full)} cells")
    print(f"Untreated-only:    {df_untreated['patient_id'].nunique()} patients, {len(df_untreated)} cells\n")

    results = []
    for label, df in [("full", df_full), ("untreated_only", df_untreated)]:
        print(f"--- Subset: {label} ---")
        hypoxia_acinar_correlation(df, label, results)
        patient_state = get_patient_state(df)
        state_counts = pd.Series(patient_state).value_counts().to_dict()
        print(f"  Patient states: {state_counts}")
        lipid_cell_of_origin(df, patient_state, label, results)
        emt_cell_of_origin(df, patient_state, label, results)
        print(f"  Scoring CAF subtypes ({label})...")
        caf_subtype_proportions(df, patient_state, untreated_patients, label, results)
        print()

    results_df = pd.DataFrame(results)
    out_path = os.path.join(TABLES_DIR, "sensitivity_GSE202051_treatment_status.tsv")
    results_df.to_csv(out_path, sep="\t", index=False)
    print(f"Results saved: {out_path}\n")

    print("=== Summary: full vs untreated-only ===")
    for analysis in results_df["analysis"].unique():
        sub = results_df[results_df["analysis"] == analysis]
        print(f"\n[{analysis}]")
        print(sub.drop(columns=["analysis"]).to_string(index=False))

    print("\n=== Sensitivity Analysis Complete ===")


if __name__ == "__main__":
    main()
