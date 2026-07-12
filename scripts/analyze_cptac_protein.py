"""
analyze_cptac_protein.py
Protein-level validation of lipid rewiring genes using real CPTAC-PDA data.
Uses the cptac Python package (source: umich proteomics + washu transcriptomics).
Group assignment via hypoxia/acinar scores computed from matched RNA-seq.
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
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
FIGURES_DIR = os.path.join(BASE_DIR, "results", "figures")
PROT_PROCESSED = os.path.join(BASE_DIR, "data", "processed", "proteomics")
os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(PROT_PROCESSED, exist_ok=True)

# Target proteins from Phase 2 gene sets
LIPID_PROTEINS = {
    "lipid_synthesis_srebp": {
        "genes": ["SREBF1", "FASN", "ACACA", "ACLY", "HMGCR", "SQLE"],
        "expected_direction": "up",
    },
    "desaturation_elongation": {
        "genes": ["SCD", "FADS1", "FADS2", "ELOVL6"],
        "expected_direction": "up",
    },
    "fatty_acid_uptake_oxidation": {
        "genes": ["CD36", "FABP4", "FABP5", "CPT1A", "CPT1B", "ACADL", "HADHA"],
        "expected_direction": "down",
    },
}

# Phase 2 gene sets for group assignment via transcriptomics
HYPOXIA_GENES = ["EPAS1", "VEGFA", "CA9", "ADM", "EGLN3", "LOX", "SLC2A1", "BNIP3", "ANGPT2"]
ACINAR_GENES = ["PTF1A", "BHLHA15", "RBPJL", "CPA1", "CPA2", "PRSS1", "PRSS2", "CEL", "CELA3A", "AMY2A", "REG1A"]


def load_cptac_dataset():
    """Load the shared CPTAC PDAC object plus transcriptomics/purity once."""
    import cptac
    print("  Loading CPTAC PDAC dataset...")
    pdac = cptac.Pdac()

    # WashU transcriptomics (shared group assignment for both sources)
    print("  Fetching washu transcriptomics...")
    trans_raw = pdac.get_transcriptomics(source="washu")
    trans = trans_raw[~trans_raw.index.str.endswith(".N")].copy()
    print(f"  Transcriptomics: {trans.shape[0]} tumor samples x {trans.shape[1]} genes")

    # ESTIMATE tumor purity
    print("  Fetching tumor purity (ESTIMATE)...")
    try:
        purity_raw = pdac.get_tumor_purity(source="washu")
        purity_df = purity_raw[~purity_raw.index.str.endswith(".N")].copy()
        print(f"  Purity: {purity_df.shape[0]} samples")
    except Exception as e:
        print(f"  Purity load failed: {e}")
        purity_df = None

    return pdac, trans, purity_df


def load_proteomics(pdac, source):
    """Load and clean proteomics for one source ('umich' or 'bcm')."""
    print(f"  Fetching {source} proteomics...")
    prot_raw = pdac.get_proteomics(source=source)
    prot = prot_raw[~prot_raw.index.str.endswith(".N")].copy()

    # Flatten MultiIndex columns (gene, site) → gene name; average duplicates
    if isinstance(prot.columns, pd.MultiIndex):
        prot.columns = prot.columns.get_level_values(0)
        prot = prot.T.groupby(level=0).mean().T

    print(f"  {source} proteomics: {prot.shape[0]} tumor samples x {prot.shape[1]} proteins")
    return prot


def assign_groups(trans):
    """Assign hypoxia-high/acinar-low groups using Phase 2 gene set scoring."""
    print("  Computing hypoxia and acinar scores from transcriptomics...")

    # Score hypoxia
    hyp_avail = [g for g in HYPOXIA_GENES if g in trans.columns]
    acin_avail = [g for g in ACINAR_GENES if g in trans.columns]
    print(f"  Hypoxia genes found: {len(hyp_avail)}/{len(HYPOXIA_GENES)}: {hyp_avail}")
    print(f"  Acinar genes found: {len(acin_avail)}/{len(ACINAR_GENES)}: {acin_avail}")

    hypoxia_score = trans[hyp_avail].mean(axis=1)
    acinar_score = trans[acin_avail].mean(axis=1)

    # Median split (same as Phase 2)
    hyp_median = hypoxia_score.median()
    acin_median = acinar_score.median()
    hyp_high = hypoxia_score >= hyp_median
    acin_low = acinar_score < acin_median

    groups = pd.Series("other", index=trans.index)
    groups[hyp_high & acin_low] = "aggressive"
    groups[~hyp_high & ~acin_low] = "reference"

    group_df = pd.DataFrame({
        "hypoxia_score": hypoxia_score,
        "acinar_score": acinar_score,
        "group": groups,
        "aggressive": (groups == "aggressive").astype(int),
    })

    print(f"  Groups: aggressive={int((groups=='aggressive').sum())}, "
          f"reference={int((groups=='reference').sum())}, "
          f"other={int((groups=='other').sum())}")
    return group_df


def run_protein_comparison(prot, group_df):
    """Compare protein abundance between aggressive and reference groups."""
    # Align samples
    common = prot.index.intersection(group_df.index)
    print(f"  Samples with both proteomics + group assignment: {len(common)}")

    prot_aligned = prot.loc[common]
    grp_aligned = group_df.loc[common]

    agg_idx = grp_aligned["group"] == "aggressive"
    ref_idx = grp_aligned["group"] == "reference"
    print(f"  Aggressive: {agg_idx.sum()}, Reference: {ref_idx.sum()}")

    results = []
    for gene_set, info in LIPID_PROTEINS.items():
        expected = info["expected_direction"]
        for gene in info["genes"]:
            if gene not in prot_aligned.columns:
                print(f"    MISSING: {gene}")
                continue

            agg_clean = prot_aligned.loc[agg_idx, gene].dropna()
            ref_clean = prot_aligned.loc[ref_idx, gene].dropna()

            if len(agg_clean) < 3 or len(ref_clean) < 3:
                continue

            stat, pval = stats.ranksums(agg_clean, ref_clean)
            median_agg = float(agg_clean.median())
            median_ref = float(ref_clean.median())
            effect_size = (median_agg - median_ref)
            observed = "up" if median_agg > median_ref else "down"
            reproduced = (observed == expected)

            results.append({
                "gene_set": gene_set,
                "protein": gene,
                "comparison": "aggressive_vs_reference",
                "n_aggressive": len(agg_clean),
                "n_reference": len(ref_clean),
                "median_aggressive": round(median_agg, 4),
                "median_reference": round(median_ref, 4),
                "effect_size": round(effect_size, 4),
                "p_value": round(pval, 6),
                "fdr": np.nan,
                "expected_direction": expected,
                "observed_direction": observed,
                "direction_reproduced_at_protein_level": reproduced,
                "is_simulated": False,
            })

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(results)
    _, fdr, _, _ = multipletests(df["p_value"].values, method="fdr_bh")
    df["fdr"] = fdr.round(4)
    return df


def make_figure(prot, group_df, results_df):
    """Generate Figure 3F: protein abundance by group for lipid genes."""
    common = prot.index.intersection(group_df.index)
    prot_a = prot.loc[common]
    grp_a = group_df.loc[common]
    plot_data = prot_a.join(grp_a[["group"]]).query("group in ['aggressive','reference']")

    proteins_to_plot = results_df["protein"].tolist()
    n_proteins = len(proteins_to_plot)
    if n_proteins == 0:
        return

    cols = min(5, n_proteins)
    rows = (n_proteins + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3, rows * 3.5))
    axes = np.array(axes).flatten() if n_proteins > 1 else [axes]

    palette = {"aggressive": "#e74c3c", "reference": "#3498db"}

    for i, row in results_df.iterrows():
        if i >= len(axes):
            break
        gene = row["protein"]
        ax = axes[i]
        sub = plot_data[["group", gene]].dropna()
        sns.boxplot(data=sub, x="group", y=gene, palette=palette, ax=ax,
                    order=["aggressive", "reference"], width=0.5, fliersize=2)
        direction = "↑" if row["observed_direction"] == "up" else "↓"
        expected_symbol = "✓" if row["direction_reproduced_at_protein_level"] else "✗"
        fdr_label = f"FDR={row['fdr']:.3f}" if row["fdr"] < 0.1 else f"FDR={row['fdr']:.2f}"
        color = "#27ae60" if row["direction_reproduced_at_protein_level"] else "#c0392b"
        ax.set_title(f"{gene} {direction}{expected_symbol}\n{fdr_label}", fontsize=9, color=color)
        ax.set_xlabel("")
        ax.set_ylabel("Protein abundance\n(log2 ratio)", fontsize=8)
        ax.tick_params(axis="x", labelsize=8)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle("CPTAC-PDA Lipid Enzyme Protein Abundance\nAggressive vs Reference (H-hi/A-lo vs H-lo/A-hi)\n"
                 "[REAL DATA — umich proteomics + washu transcriptomics]",
                 fontsize=10, y=1.01)
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "Figure3F_CPTAC_lipid_protein_by_group.pdf")
    plt.savefig(fig_path, bbox_inches="tight")
    plt.close()
    print(f"  Figure saved: {fig_path}")


def print_summary(results_df, source):
    print(f"\n--- {source} Summary ---")
    for gs in LIPID_PROTEINS:
        sub = results_df[results_df["gene_set"] == gs]
        n_repro = sub["direction_reproduced_at_protein_level"].sum()
        n_sig_repro = ((sub["fdr"] < 0.05) & sub["direction_reproduced_at_protein_level"]).sum()
        exp = LIPID_PROTEINS[gs]["expected_direction"]
        print(f"  {gs} (expected {exp}): {n_repro}/{len(sub)} correct direction, "
              f"{n_sig_repro}/{len(sub)} FDR<0.05+correct")
    total_repro = results_df["direction_reproduced_at_protein_level"].sum()
    total = len(results_df)
    sig_repro = ((results_df["fdr"] < 0.05) & results_df["direction_reproduced_at_protein_level"]).sum()
    print(f"  OVERALL: {total_repro}/{total} directionally concordant, {sig_repro}/{total} FDR<0.05+concordant")


def make_replication_figure(umich_df, bcm_df):
    """Figure 3F-BCM: side-by-side effect sizes for both sources."""
    proteins = sorted(set(umich_df["protein"]) | set(bcm_df["protein"]))
    umich_map = umich_df.set_index("protein")
    bcm_map = bcm_df.set_index("protein")

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(proteins))
    w = 0.35

    def _bar_vals(df_map, proteins):
        effects, colors, alphas = [], [], []
        for p in proteins:
            if p in df_map.index:
                row = df_map.loc[p]
                effects.append(row["effect_size"])
                repro = row["direction_reproduced_at_protein_level"]
                sig = row["fdr"] < 0.05
                colors.append("#27ae60" if (repro and sig) else
                              "#f39c12" if repro else "#c0392b")
                alphas.append(1.0 if sig else 0.5)
            else:
                effects.append(0)
                colors.append("grey")
                alphas.append(0.3)
        return effects, colors, alphas

    u_eff, u_col, u_alp = _bar_vals(umich_map, proteins)
    b_eff, b_col, b_alp = _bar_vals(bcm_map, proteins)

    for i, (ue, uc, ua) in enumerate(zip(u_eff, u_col, u_alp)):
        ax.bar(x[i] - w/2, ue, width=w, color=uc, alpha=ua, label="umich" if i == 0 else "")
    for i, (be, bc, ba) in enumerate(zip(b_eff, b_col, b_alp)):
        ax.bar(x[i] + w/2, be, width=w, color=bc, alpha=ba, hatch="//",
               label="bcm" if i == 0 else "")

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(proteins, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Effect size (median aggressive − reference)", fontsize=10)
    ax.set_title("CPTAC-PDA Lipid Enzyme Protein Abundance: umich vs BCM Replication\n"
                 "Green=FDR<0.05+concordant  Orange=concordant  Red=discordant  Striped=BCM",
                 fontsize=10)

    # Shade gene set regions
    boundaries = []
    pos = 0
    for gs, info in LIPID_PROTEINS.items():
        found = [p for p in info["genes"] if p in proteins]
        boundaries.append((pos, pos + len(found), gs, info["expected_direction"]))
        pos += len(found)
    shades = ["#eaf4fb", "#fef9e7", "#fdf2f8"]
    for idx, (start, end, gs, exp) in enumerate(boundaries):
        ax.axvspan(start - 0.5, end - 0.5, alpha=0.15, color=shades[idx % 3])
        ax.text((start + end - 1) / 2, ax.get_ylim()[1] * 0.95,
                f"{gs}\n(expected {exp})", ha="center", va="top", fontsize=7)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#27ae60", label="FDR<0.05 + concordant"),
        Patch(facecolor="#f39c12", label="Concordant, ns"),
        Patch(facecolor="#c0392b", label="Discordant"),
        Patch(facecolor="grey", label="Not detected"),
        Patch(facecolor="white", edgecolor="black", label="Solid=umich / Striped=BCM"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8)
    plt.tight_layout()
    path = os.path.join(FIGURES_DIR, "Figure3F_CPTAC_lipid_protein_by_group.pdf")
    plt.savefig(path, bbox_inches="tight")
    plt.close()
    print(f"  Figure saved: {path}")


def make_concordance_table(umich_df, bcm_df):
    """Merge both sources into a single replication table."""
    u = umich_df[["gene_set","protein","expected_direction",
                  "median_aggressive","median_reference","effect_size",
                  "p_value","fdr","observed_direction",
                  "direction_reproduced_at_protein_level"]].copy()
    u.columns = ["gene_set","protein","expected_direction",
                 "umich_median_agg","umich_median_ref","umich_effect",
                 "umich_pval","umich_fdr","umich_direction","umich_reproduced"]

    b = bcm_df[["protein","median_aggressive","median_reference","effect_size",
                "p_value","fdr","observed_direction",
                "direction_reproduced_at_protein_level"]].copy()
    b.columns = ["protein","bcm_median_agg","bcm_median_ref","bcm_effect",
                 "bcm_pval","bcm_fdr","bcm_direction","bcm_reproduced"]

    merged = u.merge(b, on="protein", how="outer")
    merged["replicated_both"] = (
        merged["umich_reproduced"].fillna(False) &
        merged["bcm_reproduced"].fillna(False)
    )
    merged["significant_both"] = (
        (merged["umich_fdr"].fillna(1) < 0.05) &
        (merged["bcm_fdr"].fillna(1) < 0.05)
    )
    return merged


def main():
    print("=== CPTAC-PDA Protein Validation — umich + BCM Replication ===\n")

    # Load shared resources once
    pdac, trans, purity_df = load_cptac_dataset()

    # Compute group assignments from transcriptomics (shared)
    print("\nAssigning groups from transcriptomics...")
    group_df = assign_groups(trans)
    group_df.to_csv(os.path.join(PROT_PROCESSED, "cptac_pda_group_assignments.tsv"), sep="\t")

    all_results = {}

    for source in ["umich", "bcm"]:
        print(f"\n{'='*50}")
        print(f"  Source: {source}")
        print(f"{'='*50}")
        prot = load_proteomics(pdac, source)

        # Save lipid protein matrix
        all_target = [g for gs in LIPID_PROTEINS.values() for g in gs["genes"]]
        found_cols = [g for g in all_target if g in prot.columns]
        prot[found_cols].to_csv(
            os.path.join(PROT_PROCESSED, f"cptac_pda_{source}_lipid_protein_matrix.tsv"), sep="\t"
        )

        print(f"\nRunning protein-level comparisons ({source})...")
        results_df = run_protein_comparison(prot, group_df)
        results_df["source"] = source

        if results_df.empty:
            print(f"  No results for {source}.")
            continue

        # Save per-source table
        results_df.to_csv(
            os.path.join(TABLES_DIR, f"cptac_{source}_lipid_protein_statistics.tsv"),
            sep="\t", index=False
        )
        print_summary(results_df, source)
        all_results[source] = results_df

        # Per-source figure (simple)
        make_figure(prot, group_df, results_df)

    # Combined figure + replication table
    if "umich" in all_results and "bcm" in all_results:
        print("\nBuilding cross-source replication table...")
        concordance = make_concordance_table(all_results["umich"], all_results["bcm"])
        out_path = os.path.join(TABLES_DIR, "figure3F_cptac_lipid_protein_statistics.tsv")
        concordance.to_csv(out_path, sep="\t", index=False)
        print(f"Replication table saved: {out_path}")

        print("\nGenerating combined replication figure...")
        make_replication_figure(all_results["umich"], all_results["bcm"])

        # Print replication summary
        print("\n--- Cross-Source Replication Summary ---")
        n_repro_both = concordance["replicated_both"].sum()
        n_sig_both = concordance["significant_both"].sum()
        n_total = len(concordance)
        print(f"  Proteins in both sources: {n_total}")
        print(f"  Concordant direction in BOTH: {n_repro_both}/{n_total}")
        print(f"  FDR<0.05 + concordant in BOTH: {n_sig_both}/{n_total}")
        print("\n  Protein-level replication:")
        for _, row in concordance.iterrows():
            u_sig = "FDR<0.05" if row.get("umich_fdr", 1) < 0.05 else "ns"
            b_sig = "FDR<0.05" if row.get("bcm_fdr", 1) < 0.05 else "ns"
            rep = "REPLICATED" if row["replicated_both"] else ""
            print(f"    {row['protein']:8s}  umich: {str(row.get('umich_direction','?')):4s} {u_sig:10s}  "
                  f"bcm: {str(row.get('bcm_direction','?')):4s} {b_sig:10s}  {rep}")

    # Save ESTIMATE purity
    if purity_df is not None:
        purity_df.to_csv(os.path.join(PROT_PROCESSED, "cptac_pda_estimate_purity.tsv"), sep="\t")

    print("\n=== CPTAC Protein Analysis Complete ===")


if __name__ == "__main__":
    main()
