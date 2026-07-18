"""
analyze_acadl_moffitt_subtype.py
Links ACADL to the established Moffitt et al. 2015 (Nat Genet, "Virtual
Microdissection") basal-like/classical PDAC subtype call, rather than only
this pipeline's own composite hypoxia-high/acinar-low signature.

GSE71729 IS the Moffitt et al. cohort -- its clinical annotation already
carries the original published tumor subtype call
(tumor_subtype_0na_1classical_2basal: 1=classical, 2=basal) and the paired
stroma subtype call (stroma_subtype_0na_1low_2normal_3activated), for all
145 samples (0/NA does not occur in this cohort). No new data acquisition or
subtype classifier of our own is needed -- this is a real, first-party,
citable call already sitting unused in the processed clinical table.

Rationale: the 2026-07-13 novelty check flagged that this pipeline's own
composite signature overlaps conceptually with, and should be positioned
relative to, established PDAC subtype frameworks (Bailey ADEX, Moffitt
basal/classical) rather than presented as a new axis. Testing ACADL directly
against Moffitt's own subtype call answers the reviewer question "is this
just re-discovering basal-like PDAC?" head-on, and if so, upgrades "ACADL is
reduced in our own study-specific aggressive group" to the more citable
"ACADL is reduced in basal-like PDAC" (Moffitt et al.'s own classification,
CPTAC/TCGA-familiar terminology).

Two tests:
  1. ACADL expression: basal-like vs classical (direct, primary).
  2. Crosstab: does this pipeline's own hypoxia-high/acinar-low "aggressive"
     call actually enrich for basal-like tumors, or is it an orthogonal axis?
     (context for the Discussion, not required for claim #1 to hold either way.)
"""

import os
import numpy as np
import pandas as pd
import random
import warnings
warnings.filterwarnings("ignore")

np.random.seed(1234)
random.seed(1234)

from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESEARCH_DIR = os.path.dirname(BASE_DIR)
PHASE2_DIR = os.path.join(RESEARCH_DIR, "pdac_phase2_validation", "data", "processed", "validation")
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
FIGURES_DIR = os.path.join(BASE_DIR, "results", "figures")
os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

COHORT = "GSE71729"
SUBTYPE_MAP = {1: "classical", 2: "basal"}


def main():
    print("=== ACADL vs Moffitt Basal/Classical Subtype (GSE71729) ===\n")

    clin_path = os.path.join(PHASE2_DIR, f"{COHORT}_clinical_clean.tsv")
    expr_path = os.path.join(PHASE2_DIR, f"{COHORT}_expression_gene_symbol.tsv.gz")
    scores_path = os.path.join(PHASE2_DIR, f"{COHORT}_sample_scores.tsv")

    clin = pd.read_csv(clin_path, sep="\t", index_col=0)
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    scores = pd.read_csv(scores_path, sep="\t", index_col=0)

    if "ACADL" not in expr.index:
        print("ACADL not found in expression matrix. Exiting.")
        return

    df = pd.DataFrame(index=clin.index)
    df["tumor_subtype_code"] = clin["tumor_subtype_0na_1classical_2basal"]
    df["stroma_subtype_code"] = clin["stroma_subtype_0na_1low_2normal_3activated"]
    df["ACADL"] = expr.loc["ACADL"].reindex(df.index)
    df["aggressive"] = scores["aggressive"].reindex(df.index)
    df["group"] = scores["group"].reindex(df.index)

    df = df.dropna(subset=["tumor_subtype_code", "ACADL"])
    df = df[df["tumor_subtype_code"].isin([1, 2])]
    df["subtype"] = df["tumor_subtype_code"].map(SUBTYPE_MAP)
    print(f"Samples with tumor subtype + ACADL: {len(df)}")
    print(df["subtype"].value_counts().to_string())

    # --- Test 1: ACADL basal vs classical ---
    basal = df.loc[df["subtype"] == "basal", "ACADL"]
    classical = df.loc[df["subtype"] == "classical", "ACADL"]
    stat, pval = stats.ranksums(basal, classical)
    direction = "lower_in_basal" if basal.median() < classical.median() else "higher_in_basal"
    effect = float(basal.median() - classical.median())

    print(f"\n--- ACADL: basal (n={len(basal)}) vs classical (n={len(classical)}) ---")
    print(f"  Median ACADL: basal={basal.median():.4f}, classical={classical.median():.4f}")
    print(f"  Direction: {direction}, effect size (median diff)={effect:.4f}")
    print(f"  Wilcoxon rank-sum p={pval:.6g}")

    result_row = {
        "cohort": COHORT,
        "n_basal": len(basal),
        "n_classical": len(classical),
        "median_ACADL_basal": round(float(basal.median()), 4),
        "median_ACADL_classical": round(float(classical.median()), 4),
        "effect_size": round(effect, 4),
        "direction": direction,
        "p_value": pval,
    }

    # --- Test 2: does the pipeline's own "aggressive" call enrich for basal? ---
    df_agg = df.dropna(subset=["aggressive"])
    df_agg = df_agg[df_agg["group"].isin(["Hypoxia-high / Acinar-low", "Hypoxia-low / Acinar-high",
                                            "Hypoxia-high / Acinar-high", "Hypoxia-low / Acinar-low"])] \
             if "group" in df_agg.columns else df_agg
    crosstab = pd.crosstab(df_agg["aggressive"].astype(int), df_agg["subtype"])
    print(f"\n--- Crosstab: this pipeline's 'aggressive' call vs Moffitt subtype (n={len(df_agg)}) ---")
    print(crosstab.to_string())
    odds_ratio, fisher_p = stats.fisher_exact(crosstab.values) if crosstab.shape == (2, 2) else (np.nan, np.nan)
    print(f"  Fisher's exact test: OR={odds_ratio:.3f}, p={fisher_p:.4g}" if not np.isnan(fisher_p) else "  Fisher's exact test: not applicable (table shape)")

    crosstab_out = crosstab.copy()
    crosstab_out.to_csv(os.path.join(TABLES_DIR, "acadl_aggressive_vs_moffitt_subtype_crosstab.tsv"), sep="\t")

    # --- Save primary result ---
    result_df = pd.DataFrame([result_row])
    out_path = os.path.join(TABLES_DIR, "acadl_moffitt_subtype.tsv")
    result_df.to_csv(out_path, sep="\t", index=False)
    print(f"\nPrimary result saved: {out_path}")

    # --- Figure: boxplot ACADL by subtype ---
    fig, ax = plt.subplots(figsize=(5, 4.5))
    palette = {"basal": "#e74c3c", "classical": "#3498db"}
    sns.boxplot(data=df, x="subtype", y="ACADL", order=["basal", "classical"], palette=palette, ax=ax, width=0.5, fliersize=2)
    sns.stripplot(data=df, x="subtype", y="ACADL", order=["basal", "classical"], color="black", alpha=0.3, size=3, ax=ax)
    ax.set_title(f"GSE71729 (Moffitt et al. 2015 cohort)\nACADL by published subtype, p={pval:.3g}", fontsize=10)
    ax.set_xlabel("")
    ax.set_ylabel("ACADL expression")
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "Figure_ACADL_by_Moffitt_subtype.pdf")
    plt.savefig(fig_path, bbox_inches="tight")
    plt.close()
    print(f"Figure saved: {fig_path}")

    print("\n=== ACADL vs Moffitt Subtype Complete ===")


if __name__ == "__main__":
    main()
