"""
analyze_hypoxia_acinar_cooccurrence.py
Test whether hypoxia-high and acinar-low co-occur in the SAME malignant cells.
Quadrant analysis + correlation tests.
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(BASE_DIR, "config")
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
FIGURES_DIR = os.path.join(BASE_DIR, "results", "figures")


def load_configs():
    with open(os.path.join(CONFIG_DIR, "singlecell_cohorts.yml")) as f:
        return yaml.safe_load(f)


def load_scores(cohort):
    scores_file = os.path.join(BASE_DIR, cohort["scores_file"])
    if not os.path.exists(scores_file):
        return None
    df = pd.read_csv(scores_file, sep="\t")
    return df


def assign_quadrant(hypoxia_scores, acinar_scores):
    """Assign cells to quadrants based on median-split."""
    hyp_med = np.median(hypoxia_scores)
    acin_med = np.median(acinar_scores)
    quadrants = []
    for h, a in zip(hypoxia_scores, acinar_scores):
        if h >= hyp_med and a < acin_med:
            quadrants.append("hypoxia_high_acinar_low")
        elif h < hyp_med and a >= acin_med:
            quadrants.append("hypoxia_low_acinar_high")
        elif h >= hyp_med and a >= acin_med:
            quadrants.append("hypoxia_high_acinar_high")
        else:
            quadrants.append("hypoxia_low_acinar_low")
    return np.array(quadrants)


def analyze_cohort(cohort, all_results):
    name = cohort["name"]
    df = load_scores(cohort)
    if df is None:
        print(f"  No scores file for {name}. Skipping.")
        return None

    is_simulated = df["is_simulated"].iloc[0] if "is_simulated" in df.columns else True

    # Restrict to malignant epithelial cells
    mal_df = df[df["cell_type"] == "malignant_epithelial"].copy()
    if len(mal_df) < 20:
        print(f"  WARNING: Only {len(mal_df)} malignant cells in {name}.")

    # Correlations
    hyp_scores = mal_df["hypoxia_score"].values
    acin_scores = mal_df["acinar_identity_score"].values

    pearson_r, pearson_p = stats.pearsonr(hyp_scores, acin_scores) if len(mal_df) >= 3 else (np.nan, np.nan)
    spearman_r, spearman_p = stats.spearmanr(hyp_scores, acin_scores) if len(mal_df) >= 3 else (np.nan, np.nan)

    # Quadrant assignment
    mal_df = mal_df.copy()
    mal_df["quadrant"] = assign_quadrant(hyp_scores, acin_scores)

    # Fraction of cells per patient in each quadrant
    patient_quadrant = mal_df.groupby(["patient_id", "quadrant"]).size().unstack(fill_value=0)
    patient_quadrant = patient_quadrant.div(patient_quadrant.sum(axis=1), axis=0)

    # Overall quadrant fractions
    quadrant_counts = mal_df["quadrant"].value_counts()
    quadrant_fracs = (quadrant_counts / quadrant_counts.sum()).to_dict()

    # Fraction in hypoxia_high_acinar_low quadrant
    frac_target = quadrant_fracs.get("hypoxia_high_acinar_low", 0.0)

    result_row = {
        "cohort": name,
        "n_malignant_cells": len(mal_df),
        "n_patients": mal_df["patient_id"].nunique(),
        "pearson_r_hypoxia_acinar": round(pearson_r, 4),
        "pearson_p": round(pearson_p, 4) if not np.isnan(pearson_p) else np.nan,
        "spearman_r_hypoxia_acinar": round(spearman_r, 4),
        "spearman_p": round(spearman_p, 4) if not np.isnan(spearman_p) else np.nan,
        "frac_hypoxia_high_acinar_low": round(frac_target, 4),
        "frac_hypoxia_low_acinar_high": round(quadrant_fracs.get("hypoxia_low_acinar_high", 0.0), 4),
        "frac_hypoxia_high_acinar_high": round(quadrant_fracs.get("hypoxia_high_acinar_high", 0.0), 4),
        "frac_hypoxia_low_acinar_low": round(quadrant_fracs.get("hypoxia_low_acinar_low", 0.0), 4),
        "cooccurrence_interpretation": (
            "Co-occur in same cells (true composite state)" if frac_target > 0.30
            else "Largely separate populations (composite artifact)"
        ),
        "is_simulated": is_simulated,
    }
    all_results.append(result_row)

    # Plot
    os.makedirs(FIGURES_DIR, exist_ok=True)
    quadrant_colors = {
        "hypoxia_high_acinar_low": "#e74c3c",
        "hypoxia_low_acinar_high": "#2ecc71",
        "hypoxia_high_acinar_high": "#f39c12",
        "hypoxia_low_acinar_low": "#3498db",
    }
    fig, ax = plt.subplots(figsize=(7, 6))
    for qname, color in quadrant_colors.items():
        mask = mal_df["quadrant"] == qname
        ax.scatter(mal_df.loc[mask, "hypoxia_score"],
                   mal_df.loc[mask, "acinar_identity_score"],
                   c=color, label=f"{qname}\n(n={mask.sum()}, {quadrant_fracs.get(qname, 0):.1%})",
                   s=10, alpha=0.5)
    # Median lines
    ax.axvline(np.median(hyp_scores), color="gray", linestyle="--", linewidth=1, alpha=0.7)
    ax.axhline(np.median(acin_scores), color="gray", linestyle="--", linewidth=1, alpha=0.7)
    ax.set_xlabel("Hypoxia Score", fontsize=12)
    ax.set_ylabel("Acinar Identity Score", fontsize=12)
    title_prefix = "[SIMULATED] " if is_simulated else ""
    ax.set_title(f"{title_prefix}{name}: Malignant Cell Quadrants\n"
                 f"Pearson r={pearson_r:.3f}, p={pearson_p:.3f}", fontsize=11)
    ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, f"Figure3A_malignant_hypoxia_acinar_quadrants_{name}.pdf")
    plt.savefig(fig_path, bbox_inches="tight")
    plt.close()
    print(f"  Figure saved: {fig_path}")
    print(f"  Pearson r={pearson_r:.3f}, p={pearson_p:.4f}; frac target quadrant={frac_target:.2%}")

    return mal_df


def main():
    print("=== Hypoxia/Acinar Co-occurrence Analysis ===\n")
    config = load_configs()
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)

    all_results = []
    for cohort in config["cohorts"]:
        print(f"--- Cohort: {cohort['name']} ---")
        try:
            analyze_cohort(cohort, all_results)
        except Exception as e:
            import traceback
            print(f"  ERROR: {e}")
            traceback.print_exc()
        print()

    if all_results:
        results_df = pd.DataFrame(all_results)
        out_path = os.path.join(TABLES_DIR, "figure3A_hypoxia_acinar_cooccurrence_by_cohort.tsv")
        results_df.to_csv(out_path, sep="\t", index=False)
        print(f"Results saved: {out_path}")
        print(results_df[["cohort", "pearson_r_hypoxia_acinar", "frac_hypoxia_high_acinar_low",
                           "cooccurrence_interpretation"]].to_string(index=False))
    else:
        # Create empty result file
        pd.DataFrame(columns=[
            "cohort", "n_malignant_cells", "n_patients",
            "pearson_r_hypoxia_acinar", "pearson_p",
            "spearman_r_hypoxia_acinar", "spearman_p",
            "frac_hypoxia_high_acinar_low", "frac_hypoxia_low_acinar_high",
            "frac_hypoxia_high_acinar_high", "frac_hypoxia_low_acinar_low",
            "cooccurrence_interpretation", "is_simulated"
        ]).to_csv(os.path.join(TABLES_DIR, "figure3A_hypoxia_acinar_cooccurrence_by_cohort.tsv"), sep="\t", index=False)

    print("\n=== Analysis Complete ===")


if __name__ == "__main__":
    main()
