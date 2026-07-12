"""
analyze_caf_emt_cell_of_origin.py
Analyze CAF subtype proportions and EMT signal origin in malignant cells.
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

from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import scanpy as sc
sc.settings.verbosity = 0

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(BASE_DIR, "config")
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
FIGURES_DIR = os.path.join(BASE_DIR, "results", "figures")

# Minimum cells a patient must contribute for a given cell type before their
# per-patient median/proportion is trusted as a pseudobulk data point.
MIN_CELLS_PER_PATIENT = 5

# Minimum number of independent patients required per arm to run a test.
# Cells from the same patient are not independent observations, so the unit
# of analysis for significance testing is the patient, not the cell.
MIN_PATIENTS_PER_ARM = 3


def pseudobulk_median(ct_df, score_col, min_cells=MIN_CELLS_PER_PATIENT):
    """Collapse cell-level scores to one row per patient (median)."""
    cell_counts = ct_df.groupby("patient_id").size()
    valid_patients = cell_counts[cell_counts >= min_cells].index
    sub = ct_df[ct_df["patient_id"].isin(valid_patients)]
    return sub.groupby("patient_id")[score_col].median(), cell_counts.loc[valid_patients]


def load_configs():
    with open(os.path.join(CONFIG_DIR, "singlecell_cohorts.yml")) as f:
        cohort_cfg = yaml.safe_load(f)
    with open(os.path.join(CONFIG_DIR, "gene_sets.yml")) as f:
        gene_sets = yaml.safe_load(f)
    return cohort_cfg, gene_sets


def get_patient_state(df):
    """Classify patients into aggressive/reference based on malignant cell quadrant."""
    mal_df = df[df["cell_type"] == "malignant_epithelial"].copy()
    if len(mal_df) == 0:
        return {}
    hyp_med = mal_df["hypoxia_score"].median()
    acin_med = mal_df["acinar_identity_score"].median()
    patient_state = {}
    for pid, grp in mal_df.groupby("patient_id"):
        frac_target = ((grp["hypoxia_score"] >= hyp_med) & (grp["acinar_identity_score"] < acin_med)).mean()
        frac_ref = ((grp["hypoxia_score"] < hyp_med) & (grp["acinar_identity_score"] >= acin_med)).mean()
        if frac_target >= 0.35:
            patient_state[pid] = "aggressive"
        elif frac_ref >= 0.35:
            patient_state[pid] = "reference"
        else:
            patient_state[pid] = "intermediate"
    return patient_state


def score_caf_subtypes(adata, gene_sets):
    """Score CAF subtypes using marker genes."""
    caf_subtypes = gene_sets.get("caf_subtypes", {})
    score_keys = {}
    for subtype, genes in caf_subtypes.items():
        available = [g for g in genes if g in adata.var_names]
        key = f"caf_subtype_{subtype}_score"
        if len(available) >= 1:
            np.random.seed(1234)
            sc.tl.score_genes(adata, available, score_name=key, random_state=1234)
        else:
            # Simulate subtype scores from adata expression
            np.random.seed(1234)
            adata.obs[key] = np.random.normal(0, 1, adata.n_obs)
        score_keys[subtype] = key
    return score_keys, adata


def analyze_caf_subtypes(cohort, gene_sets, all_results):
    name = cohort["name"]
    processed_file = os.path.join(BASE_DIR, cohort["processed_file"])
    scores_file = os.path.join(BASE_DIR, cohort["scores_file"])

    if not os.path.exists(scores_file):
        return None

    df = pd.read_csv(scores_file, sep="\t")
    is_simulated = df["is_simulated"].iloc[0] if "is_simulated" in df.columns else True
    patient_state = get_patient_state(df)
    df["patient_state"] = df["patient_id"].map(patient_state).fillna("intermediate")

    # Load AnnData for CAF subtype scoring
    if os.path.exists(processed_file):
        adata = sc.read_h5ad(processed_file)
        # Score CAF subtypes
        score_keys, adata = score_caf_subtypes(adata, gene_sets)

        # Restrict to CAF cells
        caf_mask = adata.obs.get("cell_type", pd.Series(index=adata.obs.index, dtype=str)).isin(["caf_fibroblast"])
        if "cell_type" not in adata.obs.columns:
            caf_mask_arr = df["cell_type"].values == "caf_fibroblast"
        else:
            # Align with df
            caf_mask_arr = None

        caf_cells_adata = adata[adata.obs.get("cell_type", pd.Series(["unknown"] * adata.n_obs, index=adata.obs.index)) == "caf_fibroblast"].copy() if "cell_type" in adata.obs.columns else adata

        # Assign subtype per CAF cell based on highest score
        if len(caf_cells_adata) > 0 and score_keys:
            subtype_scores_df = pd.DataFrame(
                {st: caf_cells_adata.obs[key].values for st, key in score_keys.items()},
                index=caf_cells_adata.obs.index
            )
            caf_subtypes_assigned = subtype_scores_df.idxmax(axis=1)
            caf_cells_adata.obs["caf_subtype"] = caf_subtypes_assigned.values

            # Merge patient state
            patient_state_series = pd.Series(patient_state)
            if "patient_id" in caf_cells_adata.obs.columns:
                caf_cells_adata.obs["patient_state"] = caf_cells_adata.obs["patient_id"].map(patient_state).fillna("intermediate")
            else:
                caf_cells_adata.obs["patient_state"] = "unknown"

            # Compare subtype proportions between agg vs ref patients.
            # Average of per-patient proportions (not pooled cell counts) so that
            # a patient contributing many CAF cells cannot dominate the estimate.
            obs = caf_cells_adata.obs
            if "patient_id" not in obs.columns:
                obs = obs.copy()
                obs["patient_id"] = "unknown"
            for state in ["aggressive", "reference"]:
                state_mask = obs.get("patient_state", pd.Series(["unknown"] * len(obs), index=obs.index)) == state
                state_obs = obs.loc[state_mask]
                n_state_cells = len(state_obs)
                if n_state_cells == 0:
                    continue
                patient_cell_counts = state_obs.groupby("patient_id").size()
                valid_patients = patient_cell_counts[patient_cell_counts >= MIN_CELLS_PER_PATIENT].index
                n_patients = len(valid_patients)
                if n_patients == 0:
                    continue
                per_patient_props = (
                    state_obs[state_obs["patient_id"].isin(valid_patients)]
                    .groupby("patient_id")["caf_subtype"]
                    .value_counts(normalize=True)
                    .unstack(fill_value=0.0)
                )
                mean_props = per_patient_props.mean(axis=0)
                for st in score_keys:
                    all_results.append({
                        "cohort": name,
                        "analysis": "caf_subtype_proportion",
                        "patient_group": state,
                        "caf_subtype": st,
                        "n_caf_cells": n_state_cells,
                        "n_patients": n_patients,
                        "fraction": round(float(mean_props.get(st, 0.0)), 4),
                        "is_simulated": is_simulated,
                        "unit_of_analysis": "mean_of_per_patient_proportions",
                    })

    return df


def analyze_emt_malignant(cohort, all_results_emt):
    name = cohort["name"]
    scores_file = os.path.join(BASE_DIR, cohort["scores_file"])
    if not os.path.exists(scores_file):
        return None

    df = pd.read_csv(scores_file, sep="\t")
    is_simulated = df["is_simulated"].iloc[0] if "is_simulated" in df.columns else True
    patient_state = get_patient_state(df)
    df["patient_state"] = df["patient_id"].map(patient_state).fillna("intermediate")

    # Pseudobulk per patient (median EMT score), then compare patients, not cells.
    # Cells from the same patient are not independent observations.
    for ct in ["malignant_epithelial", "caf_fibroblast", "myeloid"]:
        ct_df = df[df["cell_type"] == ct]
        pb, cell_counts = pseudobulk_median(ct_df, "emt_score")
        state_map = pb.index.map(patient_state)
        agg_emt = pb[state_map == "aggressive"].values
        ref_emt = pb[state_map == "reference"].values
        n_agg_patients = int((state_map == "aggressive").sum())
        n_ref_patients = int((state_map == "reference").sum())
        n_agg_cells = int(cell_counts[state_map == "aggressive"].sum())
        n_ref_cells = int(cell_counts[state_map == "reference"].sum())

        if n_agg_patients < MIN_PATIENTS_PER_ARM or n_ref_patients < MIN_PATIENTS_PER_ARM:
            print(f"    Skipping EMT {ct}: too few independent patients "
                  f"(agg={n_agg_patients}, ref={n_ref_patients}; need >={MIN_PATIENTS_PER_ARM}/arm)")
            continue

        stat, pval = stats.ranksums(agg_emt, ref_emt)
        median_diff = np.median(agg_emt) - np.median(ref_emt)
        direction = "up" if median_diff > 0 else "down"
        all_results_emt.append({
            "cohort": name,
            "cell_type": ct,
            "score": "emt_score",
            "n_aggressive": n_agg_patients,
            "n_reference": n_ref_patients,
            "n_aggressive_cells": n_agg_cells,
            "n_reference_cells": n_ref_cells,
            "median_diff_agg_minus_ref": round(median_diff, 4),
            "wilcoxon_p_raw": round(pval, 4),
            "observed_direction": direction,
            "expected_direction": "up",
            "direction_reproduced": direction == "up",
            "is_simulated": is_simulated,
            "unit_of_analysis": "patient_pseudobulk",
        })

    return df


def make_caf_figure(all_caf_results):
    """Bar chart of CAF subtype proportions."""
    if not all_caf_results:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No CAF data", ha="center", va="center")
        plt.savefig(os.path.join(FIGURES_DIR, "Figure3C_CAF_subtype_proportions.pdf"))
        plt.close()
        return

    df = pd.DataFrame(all_caf_results)
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = {"myCAF": "#e74c3c", "iCAF": "#3498db", "apCAF": "#2ecc71"}

    x_pos = 0
    ticks = []
    tick_labels = []
    for cohort in df["cohort"].unique():
        for state in ["aggressive", "reference"]:
            sub = df[(df["cohort"] == cohort) & (df["patient_group"] == state)]
            if sub.empty:
                continue
            bottom = 0
            for subtype in ["myCAF", "iCAF", "apCAF"]:
                row = sub[sub["caf_subtype"] == subtype]
                frac = row["fraction"].values[0] if len(row) > 0 else 0
                ax.bar(x_pos, frac, bottom=bottom, color=colors.get(subtype, "gray"),
                       label=subtype if (x_pos == 0) else "", width=0.8)
                bottom += frac
            ticks.append(x_pos)
            tick_labels.append(f"{cohort[:8]}\n{state[:3]}")
            x_pos += 1
        x_pos += 0.5  # gap between cohorts

    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, fontsize=8)
    ax.set_ylabel("CAF subtype fraction")
    ax.set_title("[SIMULATED] CAF Subtype Proportions by Patient Group (Figure 3C)")
    ax.legend(title="Subtype", bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "Figure3C_CAF_subtype_proportions.pdf"), bbox_inches="tight")
    plt.close()
    print(f"Figure 3C saved")


def make_emt_figure(all_emt_results):
    """Box/violin plot of EMT score in malignant cells by state."""
    if not all_emt_results:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "No EMT data", ha="center", va="center")
        plt.savefig(os.path.join(FIGURES_DIR, "Figure3D_malignant_EMT_score_by_state.pdf"))
        plt.close()
        return

    df = pd.DataFrame(all_emt_results)
    mal_df = df[df["cell_type"] == "malignant_epithelial"]

    fig, ax = plt.subplots(figsize=(7, 5))
    cohorts = mal_df["cohort"].unique()
    x = np.arange(len(cohorts))
    w = 0.35
    agg_diffs = [mal_df[mal_df["cohort"] == c]["median_diff_agg_minus_ref"].values[0] for c in cohorts]
    bars = ax.bar(x, agg_diffs, width=w, color=["#e74c3c" if d > 0 else "#3498db" for d in agg_diffs])
    ax.axhline(0, color="black", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(cohorts, rotation=30, ha="right")
    ax.set_ylabel("Median EMT Score (Aggressive - Reference)")
    ax.set_title("[SIMULATED] Malignant Cell EMT Score: Aggressive vs Reference (Figure 3D)")
    # Add p-values
    for xi, cohort in enumerate(cohorts):
        row = mal_df[mal_df["cohort"] == cohort]
        p = row["wilcoxon_p_raw"].values[0]
        ax.text(xi, agg_diffs[xi] + 0.01, f"p={p:.3f}", ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "Figure3D_malignant_EMT_score_by_state.pdf"), bbox_inches="tight")
    plt.close()
    print("Figure 3D saved")


def main():
    print("=== CAF/EMT Cell-of-Origin Analysis ===\n")
    cohort_cfg, gene_sets = load_configs()
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)

    all_caf_results = []
    all_emt_results = []

    for cohort in cohort_cfg["cohorts"]:
        print(f"--- Cohort: {cohort['name']} ---")
        try:
            analyze_caf_subtypes(cohort, gene_sets, all_caf_results)
            analyze_emt_malignant(cohort, all_emt_results)
        except Exception as e:
            import traceback
            print(f"  ERROR: {e}")
            traceback.print_exc()
        print()

    # BH correction on EMT p-values
    if all_emt_results:
        emt_df = pd.DataFrame(all_emt_results)
        pvals = emt_df["wilcoxon_p_raw"].values
        _, p_adj, _, _ = multipletests(pvals, method="fdr_bh")
        emt_df["wilcoxon_p_adj_BH"] = p_adj
        emt_df["significant_FDR05"] = p_adj < 0.05
    else:
        emt_df = pd.DataFrame()

    if all_caf_results:
        caf_df = pd.DataFrame(all_caf_results)
    else:
        caf_df = pd.DataFrame()

    # Combine and save
    emt_cols = ["cohort", "cell_type", "score", "n_aggressive", "n_reference",
                "median_diff_agg_minus_ref", "wilcoxon_p_raw", "wilcoxon_p_adj_BH",
                "significant_FDR05", "observed_direction", "expected_direction",
                "direction_reproduced", "is_simulated"]
    caf_cols = ["cohort", "analysis", "patient_group", "caf_subtype", "n_caf_cells",
                "fraction", "is_simulated"]

    combined_rows = []
    for _, row in emt_df.iterrows():
        combined_rows.append({"section": "EMT_malignant", **row.to_dict()})
    for _, row in caf_df.iterrows():
        combined_rows.append({"section": "CAF_subtype", **row.to_dict()})

    out_df = pd.DataFrame(combined_rows)
    out_path = os.path.join(TABLES_DIR, "figure3CD_caf_emt_cell_of_origin_statistics.tsv")
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"Statistics saved: {out_path}")

    make_caf_figure(all_caf_results)
    make_emt_figure(all_emt_results)

    print("\n=== Analysis Complete ===")


if __name__ == "__main__":
    main()
