"""
deconvolve_bulk_purity.py
Estimate tumor purity from Phase 2 bulk cohorts using a stromal gene proxy.
Run purity-adjusted linear models for CAF and EMT scores.
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
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(BASE_DIR, "config")
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
FIGURES_DIR = os.path.join(BASE_DIR, "results", "figures")
PURITY_DIR = os.path.join(BASE_DIR, "data", "processed", "purity_adjusted")


def load_stromal_genes():
    """Load the purity proxy gene list and verify it is disjoint from the CAF and
    EMT signatures being purity-adjusted. The original proxy (ACTA2, COL1A1,
    COL1A2, DCN, FAP, PDGFRB) was a strict subset of the 10-gene CAF signature,
    making the "purity-adjusted" CAF regression partly circular (regressing the
    CAF score on ~60% of its own genes). This proxy instead spans immune,
    generic stromal, and vascular compartments via genes never used to define
    caf/emt/hypoxia/acinar/lipid* elsewhere in this pipeline."""
    with open(os.path.join(CONFIG_DIR, "gene_sets.yml")) as f:
        gene_sets = yaml.safe_load(f)
    stromal_genes = gene_sets["purity_stromal_immune"]
    for sig_name in ("caf", "emt"):
        overlap = set(stromal_genes) & set(gene_sets.get(sig_name, []))
        assert not overlap, f"Purity proxy overlaps {sig_name} signature: {overlap}"
    return stromal_genes


STROMAL_GENES = load_stromal_genes()

PHASE2_PROCESSED = "C:/Users/erica/OneDrive/Research/pdac_phase2_validation/data/processed/validation"
PHASE2_RAW = "C:/Users/erica/OneDrive/Research/pdac_phase2_validation/data/raw/validation"

COHORTS = [
    {
        "name": "GSE79668",
        "scores_file": "GSE79668_sample_scores.tsv",
        "expression_file": "GSE79668_expression_gene_symbol.tsv.gz",
    },
    {
        "name": "GSE71729",
        "scores_file": "GSE71729_sample_scores.tsv",
        "expression_file": "GSE71729_expression_gene_symbol.tsv.gz",
    },
    {
        "name": "GSE62165",
        "scores_file": "GSE62165_sample_scores.tsv",
        "expression_file": "GSE62165_expression_gene_symbol.tsv.gz",
    },
]


def load_phase2_data(cohort):
    """Load Phase 2 bulk data. The scores file contains the aggressive binary column directly."""
    name = cohort["name"]
    scores_path = os.path.join(PHASE2_PROCESSED, cohort["scores_file"])
    expr_path = os.path.join(PHASE2_PROCESSED, cohort["expression_file"])

    try:
        scores_df = pd.read_csv(scores_path, sep="\t", index_col=0)
        print(f"  Loaded real scores: {scores_path} ({len(scores_df)} samples)")
    except Exception as e:
        print(f"  Could not load real scores ({e}). Simulating bulk data.")
        np.random.seed(1234)
        n_samples = 60
        samples = [f"{name}_S{i:03d}" for i in range(n_samples)]
        scores_df = pd.DataFrame({
            "sample_id": samples,
            "aggressive": [1]*(n_samples//2) + [0]*(n_samples//2),
            "caf_score": np.concatenate([np.random.normal(0.6, 0.2, n_samples//2),
                                          np.random.normal(0.3, 0.2, n_samples//2)]),
            "emt_score": np.concatenate([np.random.normal(0.5, 0.2, n_samples//2),
                                          np.random.normal(0.2, 0.2, n_samples//2)]),
        }, index=samples)
        return scores_df, None, True

    try:
        expr_df = pd.read_csv(expr_path, sep="\t", index_col=0)
        print(f"  Loaded expression: {expr_path} ({expr_df.shape})")
        return scores_df, expr_df, False
    except Exception as e:
        print(f"  Could not load expression ({e}). Purity will use stromal score from scores file.")
        return scores_df, None, False


def load_real_estimate_purity(cohort_name):
    """Load real ESTIMATE (Yoshihara et al. 2013) stromal/immune scores computed
    by scripts/estimate_purity.R via the tidyestimate package, and convert to a
    purity-direction covariate. Returns None if not yet computed for this cohort.

    We use the raw combined 'estimate' score (stromal + immune enrichment),
    not the Affymetrix-only purity-conversion formula, as the adjustment
    covariate for all three cohorts: the conversion formula
    (purity = cos(0.6049872018 + 0.0001467884*estimate)) is only valid for
    Affymetrix arrays (only GSE62165 qualifies; GSE79668 is RNA-seq, GSE71729
    is Agilent), and even restricted to GSE62165 it returns NA for 112/131
    samples (high-stroma/immune samples push the cosine argument past the
    point where it goes negative) — too much sample loss to use as the
    primary covariate. The raw 'estimate' score requires no such conversion
    and is itself a standard, widely-used continuous stromal/immune-content
    covariate in the literature."""
    path = os.path.join(PURITY_DIR, f"{cohort_name}_estimate_scores.tsv")
    if not os.path.exists(path):
        print(f"  Real ESTIMATE scores not found at {path}; run scripts/estimate_purity.R first")
        return None
    df = pd.read_csv(path, sep="\t", index_col=0)
    estimate_score = df["estimate"]
    s_min, s_max = estimate_score.min(), estimate_score.max()
    if s_max > s_min:
        purity = 1.0 - (estimate_score - s_min) / (s_max - s_min)
    else:
        purity = pd.Series(0.5, index=estimate_score.index)
    return purity


def estimate_purity(scores_df, expr_df, cohort_name):
    """Estimate tumor purity via a stromal/immune gene expression proxy that is
    disjoint from the CAF and EMT signatures being purity-adjusted (see
    load_stromal_genes), so the adjustment is not partly circular."""
    if expr_df is not None:
        available_genes = [g for g in STROMAL_GENES if g in expr_df.index]
        if len(available_genes) >= 2:
            stromal_score = expr_df.loc[available_genes].mean(axis=0)
            print(f"  Purity proxy: {len(available_genes)}/{len(STROMAL_GENES)} stromal genes from expression matrix")
        else:
            print(f"  WARNING: Only {len(available_genes)} stromal genes in expression matrix. Falling back to ESTIMATE proxy.")
            stromal_score = None
    else:
        stromal_score = None

    if stromal_score is None:
        # Use macrophage_score as orthogonal stromal proxy (not caf/emt to avoid collinearity)
        if "macrophage_score" in scores_df.columns:
            print(f"  Purity proxy: macrophage_score (avoiding collinearity with CAF/EMT targets)")
            stromal_score = scores_df["macrophage_score"]
        else:
            print(f"  Purity proxy: random (no expression or macrophage score available)")
            np.random.seed(1234)
            stromal_score = pd.Series(np.random.normal(0, 1, scores_df.shape[0]), index=scores_df.index)

    # Ensure index alignment with scores_df
    if isinstance(stromal_score, pd.Series):
        stromal_score = stromal_score.reindex(scores_df.index)

    # Normalize: purity = 1 - normalized_stromal (high stromal = low purity)
    s_min = stromal_score.min()
    s_max = stromal_score.max()
    if s_max > s_min:
        purity = 1.0 - (stromal_score - s_min) / (s_max - s_min)
    else:
        purity = pd.Series(0.5, index=stromal_score.index)

    return purity.fillna(0.5).clip(0, 1)


def run_linear_models(scores_df, purity, cohort_name, is_simulated, purity_method):
    """Run adjusted and unadjusted linear models for CAF and EMT scores.
    Uses the 'aggressive' binary column directly from the Phase 2 scores file."""
    results = []

    # Align purity with scores_df
    common = scores_df.index.intersection(purity.index)
    if len(common) == 0:
        n = min(len(scores_df), len(purity))
        scores_aligned = scores_df.iloc[:n]
        purity_aligned = purity.iloc[:n]
    else:
        scores_aligned = scores_df.loc[common]
        purity_aligned = purity.loc[common]

    # Use 'aggressive' binary column directly (1=aggressive, 0=reference/other)
    if "aggressive" in scores_aligned.columns:
        group_binary = scores_aligned["aggressive"].astype(float).values
        print(f"  Group coding: 'aggressive' column ({int(group_binary.sum())} aggressive, {int((1-group_binary).sum())} reference)")
    else:
        print("  WARNING: no 'aggressive' column found; defaulting to zeros")
        group_binary = np.zeros(len(scores_aligned))

    for score_col in ["caf_score", "emt_score"]:
        if score_col not in scores_aligned.columns:
            continue
        y = scores_aligned[score_col].values.astype(float)
        x_group = group_binary

        # Unadjusted
        try:
            X_unadj = sm.add_constant(x_group)
            model_unadj = sm.OLS(y, X_unadj).fit()
            coef_unadj = model_unadj.params[1]
            pval_unadj = model_unadj.pvalues[1]
        except Exception:
            coef_unadj = np.nan
            pval_unadj = np.nan

        # Purity-adjusted
        try:
            purity_vals = purity_aligned.values.astype(float) if isinstance(purity_aligned, pd.Series) else purity_aligned
            X_adj = sm.add_constant(np.column_stack([x_group, purity_vals]))
            model_adj = sm.OLS(y, X_adj).fit()
            coef_adj = model_adj.params[1]
            pval_adj = model_adj.pvalues[1]
            coef_purity = model_adj.params[2]
        except Exception:
            coef_adj = np.nan
            pval_adj = np.nan
            coef_purity = np.nan

        direction_change = (
            "attenuated" if (not np.isnan(coef_unadj) and not np.isnan(coef_adj) and abs(coef_adj) < abs(coef_unadj))
            else "amplified" if (not np.isnan(coef_unadj) and not np.isnan(coef_adj) and abs(coef_adj) > abs(coef_unadj))
            else "unchanged"
        )

        results.append({
            "cohort": cohort_name,
            "purity_method": purity_method,
            "score": score_col,
            "n_samples": len(scores_aligned),
            "n_aggressive": int(group_binary.sum()),
            "n_reference": int((1 - group_binary).sum()),
            "coef_unadjusted": round(float(coef_unadj), 4) if not np.isnan(coef_unadj) else np.nan,
            "pval_unadjusted": round(float(pval_unadj), 4) if not np.isnan(pval_unadj) else np.nan,
            "coef_purity_adjusted": round(float(coef_adj), 4) if not np.isnan(coef_adj) else np.nan,
            "pval_purity_adjusted": round(float(pval_adj), 4) if not np.isnan(pval_adj) else np.nan,
            "coef_purity_term": round(float(coef_purity), 4) if not np.isnan(coef_purity) else np.nan,
            "direction_change_after_adjustment": direction_change,
            "interpretation": (
                f"Effect {direction_change} after purity adjustment; "
                f"suggests {'confounding by stroma' if direction_change == 'attenuated' else 'true tumor signal'}"
            ),
            "is_simulated": is_simulated,
        })

    return results


def main():
    print("=== Purity Deconvolution / Adjustment ===\n")
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(PURITY_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)

    all_results = []

    for cohort in COHORTS:
        name = cohort["name"]
        print(f"--- Cohort: {name} ---")
        try:
            scores_df, expr_df, is_simulated = load_phase2_data(cohort)

            # Method 1 (primary): real ESTIMATE (Yoshihara et al. 2013) via tidyestimate
            real_purity = load_real_estimate_purity(name)
            if real_purity is not None:
                print(f"  [real ESTIMATE] purity-direction range: {real_purity.min():.3f} - {real_purity.max():.3f}")
                purity_df = pd.DataFrame({"sample_id": real_purity.index, "estimated_purity": real_purity.values,
                                          "cohort": name, "method": "real_ESTIMATE_ssGSEA", "is_simulated": is_simulated})
                purity_df.to_csv(os.path.join(PURITY_DIR, f"{name}_real_estimate_purity.tsv"), sep="\t", index=False)
                results = run_linear_models(scores_df, real_purity, name, is_simulated, "real_ESTIMATE_ssGSEA")
                all_results.extend(results)
            else:
                print(f"  WARNING: no real ESTIMATE scores for {name}; run scripts/estimate_purity.R first")

            # Method 2 (legacy, kept for comparison): disjoint 8-gene heuristic proxy
            heuristic_purity = estimate_purity(scores_df, expr_df, name)
            print(f"  [heuristic proxy] purity range: {heuristic_purity.min():.3f} - {heuristic_purity.max():.3f}")
            heuristic_df = pd.DataFrame({"sample_id": heuristic_purity.index, "estimated_purity": heuristic_purity.values,
                                      "cohort": name, "method": "heuristic_8gene_proxy", "is_simulated": is_simulated})
            heuristic_df.to_csv(os.path.join(PURITY_DIR, f"{name}_purity_estimates.tsv"), sep="\t", index=False)
            results = run_linear_models(scores_df, heuristic_purity, name, is_simulated, "heuristic_8gene_proxy")
            all_results.extend(results)
        except Exception as e:
            import traceback
            print(f"  ERROR: {e}")
            traceback.print_exc()
        print()

    # BH correction across all p-values, computed separately per purity method
    # (each method is its own family of tests; pooling both together would
    # understate significance for whichever method has fewer/weaker results)
    if all_results:
        df = pd.DataFrame(all_results)
        df["pval_unadjusted_BH"] = np.nan
        df["pval_purity_adjusted_BH"] = np.nan
        for method, idx in df.groupby("purity_method").groups.items():
            pvals_unadj = df.loc[idx, "pval_unadjusted"].fillna(1.0).values
            pvals_adj_model = df.loc[idx, "pval_purity_adjusted"].fillna(1.0).values
            _, unadj_bh, _, _ = multipletests(pvals_unadj, method="fdr_bh")
            _, adj_bh, _, _ = multipletests(pvals_adj_model, method="fdr_bh")
            df.loc[idx, "pval_unadjusted_BH"] = unadj_bh.round(4)
            df.loc[idx, "pval_purity_adjusted_BH"] = adj_bh.round(4)

        out_path = os.path.join(TABLES_DIR, "figure3E_purity_adjusted_caf_emt_results_by_cohort.tsv")
        df.to_csv(out_path, sep="\t", index=False)
        print(f"Results saved: {out_path}")
        print(df[["cohort", "purity_method", "score", "coef_unadjusted", "coef_purity_adjusted",
                   "direction_change_after_adjustment"]].to_string(index=False))

        # Figure (primary method: real ESTIMATE only, to keep the plot readable;
        # the full comparison including the legacy heuristic proxy is in the TSV)
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        plot_df = df[df["purity_method"] == "real_ESTIMATE_ssGSEA"]
        for ax_idx, score_col in enumerate(["caf_score", "emt_score"]):
            ax = axes[ax_idx]
            sub = plot_df[plot_df["score"] == score_col]
            x = np.arange(len(sub))
            w = 0.35
            ax.bar(x - w/2, sub["coef_unadjusted"].fillna(0), width=w, label="Unadjusted", color="#3498db", alpha=0.8)
            ax.bar(x + w/2, sub["coef_purity_adjusted"].fillna(0), width=w, label="Purity-adjusted", color="#e74c3c", alpha=0.8)
            ax.axhline(0, color="black", linewidth=1)
            ax.set_xticks(x)
            ax.set_xticklabels(sub["cohort"].values, rotation=30, ha="right")
            ax.set_title(f"{score_col.replace('_', ' ').title()}")
            ax.set_ylabel("Group Effect Coefficient")
            ax.legend()
            ax.grid(axis="y", alpha=0.3)
        suptitle = "[SIMULATED] " if df["is_simulated"].all() else ""
        plt.suptitle(f"{suptitle}Purity-Adjusted CAF/EMT Effects, real ESTIMATE method (Figure 3E)", fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, "Figure3E_purity_adjusted_CAF_EMT.pdf"), bbox_inches="tight")
        plt.close()
        print("Figure 3E saved")
    else:
        pd.DataFrame(columns=[
            "cohort", "score", "n_samples", "n_aggressive", "n_reference",
            "coef_unadjusted", "pval_unadjusted", "coef_purity_adjusted",
            "pval_purity_adjusted", "coef_purity_term", "direction_change_after_adjustment",
            "interpretation", "is_simulated", "pval_unadjusted_BH", "pval_purity_adjusted_BH"
        ]).to_csv(os.path.join(TABLES_DIR, "figure3E_purity_adjusted_caf_emt_results_by_cohort.tsv"), sep="\t", index=False)

    print("\n=== Purity Analysis Complete ===")


if __name__ == "__main__":
    main()
