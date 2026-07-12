"""
analyze_lipid_cell_of_origin.py
Determine which cell types drive lipid rewiring signals.
Compares lipid scores between aggressive (hypoxia_high/acinar_low) and reference cells.
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

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(BASE_DIR, "config")
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
FIGURES_DIR = os.path.join(BASE_DIR, "results", "figures")

LIPID_SCORES = {
    "lipid_synthesis_srebp": "up",
    "desaturation_elongation": "up",
    "fatty_acid_uptake_oxidation": "down",
}

CELL_TYPES_TO_ANALYZE = ["malignant_epithelial", "caf_fibroblast", "myeloid", "endothelial"]

# Minimum cells a patient must contribute for a given cell type before their
# per-patient median is trusted as a pseudobulk data point.
MIN_CELLS_PER_PATIENT = 5

# Minimum number of independent patients required per arm to run a test.
# Cells from the same patient are not independent observations, so the unit
# of analysis for significance testing is the patient, not the cell.
MIN_PATIENTS_PER_ARM = 3


def load_configs():
    with open(os.path.join(CONFIG_DIR, "singlecell_cohorts.yml")) as f:
        return yaml.safe_load(f)


def pseudobulk_by_patient(ct_df, score_cols, min_cells=MIN_CELLS_PER_PATIENT):
    """Collapse cell-level scores to one row per patient (median), so that each
    patient — not each cell — is a single independent observation in the test."""
    cell_counts = ct_df.groupby("patient_id").size()
    valid_patients = cell_counts[cell_counts >= min_cells].index
    pb = ct_df[ct_df["patient_id"].isin(valid_patients)].groupby("patient_id")[score_cols].median()
    pb["n_cells"] = cell_counts.loc[pb.index]
    return pb


def assign_state(df):
    """Assign aggressive/reference state based on quadrant."""
    mal_df = df[df["cell_type"] == "malignant_epithelial"].copy()
    if len(mal_df) == 0:
        return df, {}
    hyp_med = mal_df["hypoxia_score"].median()
    acin_med = mal_df["acinar_identity_score"].median()

    # State is assigned per patient based on malignant cell majority
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

    df = df.copy()
    df["patient_state"] = df["patient_id"].map(patient_state).fillna("intermediate")
    return df, patient_state


def analyze_lipid_by_celltype(cohort, all_results):
    name = cohort["name"]
    scores_file = os.path.join(BASE_DIR, cohort["scores_file"])
    if not os.path.exists(scores_file):
        return

    df = pd.read_csv(scores_file, sep="\t")
    is_simulated = df["is_simulated"].iloc[0] if "is_simulated" in df.columns else True

    df, patient_state = assign_state(df)
    state_counts = pd.Series(patient_state).value_counts().to_dict()
    print(f"  {name}: patient states = {state_counts}")

    score_cols = [c + "_score" for c in LIPID_SCORES]

    for cell_type in CELL_TYPES_TO_ANALYZE:
        ct_df = df[df["cell_type"] == cell_type]
        pb = pseudobulk_by_patient(ct_df, score_cols)
        pb["state"] = pb.index.map(patient_state)

        agg_patients = pb[pb["state"] == "aggressive"]
        ref_patients = pb[pb["state"] == "reference"]
        n_agg_cells = int(agg_patients["n_cells"].sum())
        n_ref_cells = int(ref_patients["n_cells"].sum())

        if len(agg_patients) < MIN_PATIENTS_PER_ARM or len(ref_patients) < MIN_PATIENTS_PER_ARM:
            print(f"    Skipping {cell_type}: too few independent patients "
                  f"(agg={len(agg_patients)}, ref={len(ref_patients)}; "
                  f"need >={MIN_PATIENTS_PER_ARM}/arm)")
            continue

        p_values_raw = []
        test_names = []
        for score_col, expected_dir in LIPID_SCORES.items():
            col = score_col + "_score"
            a_vals = agg_patients[col].values
            r_vals = ref_patients[col].values
            stat, pval = stats.ranksums(a_vals, r_vals)
            observed_dir = "up" if np.median(a_vals) > np.median(r_vals) else "down"
            p_values_raw.append(pval)
            test_names.append((score_col, expected_dir, observed_dir, np.median(a_vals) - np.median(r_vals)))

        # BH-FDR correction
        _, p_adj, _, _ = multipletests(p_values_raw, method="fdr_bh")

        for j, (score_col, exp_dir, obs_dir, median_diff) in enumerate(test_names):
            pval_raw = p_values_raw[j]
            pval_adj = p_adj[j]
            reproduced = (obs_dir == exp_dir) and (pval_adj < 0.05)
            all_results.append({
                "cohort": name,
                "cell_type": cell_type,
                "score": score_col,
                "n_aggressive": len(agg_patients),
                "n_reference": len(ref_patients),
                "n_aggressive_cells": n_agg_cells,
                "n_reference_cells": n_ref_cells,
                "expected_direction": exp_dir,
                "observed_direction": obs_dir,
                "median_diff_agg_minus_ref": round(median_diff, 4),
                "wilcoxon_p_raw": round(pval_raw, 4),
                "wilcoxon_p_adj_BH": round(pval_adj, 4),
                "significant_FDR05": pval_adj < 0.05,
                "direction_reproduced": obs_dir == exp_dir,
                "fully_reproduced": reproduced,
                "is_simulated": is_simulated,
                "unit_of_analysis": "patient_pseudobulk",
            })

    return df


def compute_cell_intrinsic_flag(results_df):
    """Flag if malignant signal is stronger and significant vs non-malignant."""
    results_df = results_df.copy()
    results_df["cell_intrinsic_supported"] = False

    for cohort in results_df["cohort"].unique():
        for score in LIPID_SCORES.keys():
            sub = results_df[(results_df["cohort"] == cohort) & (results_df["score"] == score)]
            mal_row = sub[sub["cell_type"] == "malignant_epithelial"]
            if len(mal_row) == 0:
                continue
            mal_diff = abs(mal_row["median_diff_agg_minus_ref"].values[0])
            mal_sig = mal_row["significant_FDR05"].values[0]
            mal_repro = mal_row["direction_reproduced"].values[0]

            non_mal = sub[sub["cell_type"] != "malignant_epithelial"]
            non_mal_max_diff = non_mal["median_diff_agg_minus_ref"].abs().max() if len(non_mal) > 0 else 0

            supported = mal_sig and mal_repro and (mal_diff > non_mal_max_diff)
            idx = results_df[(results_df["cohort"] == cohort) &
                             (results_df["score"] == score) &
                             (results_df["cell_type"] == "malignant_epithelial")].index
            results_df.loc[idx, "cell_intrinsic_supported"] = supported

    return results_df


def make_boxplot_figure(all_dfs, cohorts):
    """Create boxplot of lipid scores by cell type and state."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))

    score_cols = list(LIPID_SCORES.keys())
    score_labels = {
        "lipid_synthesis_srebp": "Lipid Synthesis (SREBP)",
        "desaturation_elongation": "Desaturation/Elongation",
        "fatty_acid_uptake_oxidation": "FA Uptake/Oxidation",
    }

    ct_colors = {
        "malignant_epithelial": "#e74c3c",
        "caf_fibroblast": "#3498db",
        "myeloid": "#2ecc71",
        "endothelial": "#f39c12",
    }

    # Combine data from all cohorts
    combined = pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()

    if combined.empty:
        for ax in axes:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, "Figure3B_lipid_scores_by_celltype_and_state.pdf"))
        plt.close()
        return

    for ax_idx, score_col in enumerate(score_cols):
        ax = axes[ax_idx]
        col = score_col + "_score"
        ct_data = {}
        for ct in CELL_TYPES_TO_ANALYZE:
            for state in ["aggressive", "reference"]:
                key = f"{ct[:6]}\n{state[:3]}"
                mask = (combined["cell_type"] == ct) & (combined["patient_state"] == state)
                vals = combined.loc[mask, col].dropna().values
                ct_data[key] = vals

        keys = list(ct_data.keys())
        vals_list = [ct_data[k] for k in keys]
        # Filter empty
        pairs = [(k, v) for k, v in zip(keys, vals_list) if len(v) > 0]
        if not pairs:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            continue
        keys_filt, vals_filt = zip(*pairs)
        bp = ax.boxplot(vals_filt, labels=keys_filt, patch_artist=True, widths=0.6)
        ct_color_list = []
        for ct in CELL_TYPES_TO_ANALYZE:
            for _ in ["aggressive", "reference"]:
                ct_color_list.append(ct_colors[ct])
        for patch, color in zip(bp["boxes"], ct_color_list[:len(bp["boxes"])]):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_title(score_labels[score_col], fontsize=10)
        ax.set_ylabel("Score", fontsize=9)
        ax.tick_params(axis="x", labelsize=7)
        ax.grid(axis="y", alpha=0.3)

    # Note simulated
    suptitle = "Lipid Scores by Cell Type and State [SIMULATED DATA]" if (
        combined["is_simulated"].all() if "is_simulated" in combined.columns else True
    ) else "Lipid Scores by Cell Type and State"
    plt.suptitle(suptitle, fontsize=12)
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "Figure3B_lipid_scores_by_celltype_and_state.pdf")
    plt.savefig(fig_path, bbox_inches="tight")
    plt.close()
    print(f"Figure saved: {fig_path}")


def main():
    print("=== Lipid Cell-of-Origin Analysis ===\n")
    config = load_configs()
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)

    all_results = []
    all_dfs = []

    for cohort in config["cohorts"]:
        print(f"--- Cohort: {cohort['name']} ---")
        try:
            df = analyze_lipid_by_celltype(cohort, all_results)
            if df is not None:
                all_dfs.append(df)
        except Exception as e:
            import traceback
            print(f"  ERROR: {e}")
            traceback.print_exc()
        print()

    if all_results:
        results_df = pd.DataFrame(all_results)
        results_df = compute_cell_intrinsic_flag(results_df)
        out_path = os.path.join(TABLES_DIR, "figure3B_lipid_cell_of_origin_statistics.tsv")
        results_df.to_csv(out_path, sep="\t", index=False)
        print(f"Statistics saved: {out_path}")
        make_boxplot_figure(all_dfs, config["cohorts"])
    else:
        cols = ["cohort", "cell_type", "score", "n_aggressive", "n_reference",
                "n_aggressive_cells", "n_reference_cells",
                "expected_direction", "observed_direction", "median_diff_agg_minus_ref",
                "wilcoxon_p_raw", "wilcoxon_p_adj_BH", "significant_FDR05",
                "direction_reproduced", "fully_reproduced", "cell_intrinsic_supported",
                "is_simulated", "unit_of_analysis"]
        pd.DataFrame(columns=cols).to_csv(
            os.path.join(TABLES_DIR, "figure3B_lipid_cell_of_origin_statistics.tsv"), sep="\t", index=False)
        make_boxplot_figure([], config["cohorts"])

    print("\n=== Analysis Complete ===")


if __name__ == "__main__":
    main()
