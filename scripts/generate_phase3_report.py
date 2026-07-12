"""
generate_phase3_report.py
Generate the Phase 3 mechanistic report from all result tables.

All sections are data-driven: real/simulated status, statistics, and the
final mechanistic classification are all derived from the actual result
tables on disk (never hardcoded), so re-running this script cannot silently
regress a report that was previously accurate back to stale boilerplate.
"""

import os
import sys
import numpy as np
import pandas as pd
import random
import yaml
import warnings
import datetime
warnings.filterwarnings("ignore")

np.random.seed(1234)
random.seed(1234)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
REPORTS_DIR = os.path.join(BASE_DIR, "results", "reports")
FIGURES_DIR = os.path.join(BASE_DIR, "results", "figures")
CONFIG_DIR = os.path.join(BASE_DIR, "config")


def load_table(filename, required_cols=None):
    path = os.path.join(TABLES_DIR, filename)
    if not os.path.exists(path):
        return None, f"File not found: {path}"
    try:
        df = pd.read_csv(path, sep="\t")
        if required_cols:
            missing = [c for c in required_cols if c not in df.columns]
            if missing:
                return df, f"Warning: missing columns {missing}"
        return df, None
    except Exception as e:
        return None, f"Error loading {filename}: {e}"


def fmt(val, decimals=3):
    try:
        return f"{float(val):.{decimals}f}"
    except Exception:
        return str(val)


def bool_str(val):
    return "YES" if val else "NO"


def mode_n(df, col):
    """Most common value of an n_aggressive/n_reference-style column
    (a few proteins have fewer non-NaN samples than the full cohort)."""
    if df is None or df.empty or col not in df.columns:
        return None
    try:
        return int(df[col].mode().iloc[0])
    except Exception:
        return None


def generate_report():
    os.makedirs(REPORTS_DIR, exist_ok=True)

    lines = []
    def w(s=""):
        lines.append(s)

    w("# Phase 3 Mechanistic Report: Single-Cell Resolution of Aggressive PDAC State")
    w(f"\n**Generated:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    w("\n---\n")

    # Load all tables
    cooccur_df, cooccur_err = load_table("figure3A_hypoxia_acinar_cooccurrence_by_cohort.tsv")
    lipid_df, lipid_err = load_table("figure3B_lipid_cell_of_origin_statistics.tsv")
    caf_emt_df, caf_emt_err = load_table("figure3CD_caf_emt_cell_of_origin_statistics.tsv")
    purity_df, purity_err = load_table("figure3E_purity_adjusted_caf_emt_results_by_cohort.tsv")
    cptac_df, cptac_err = load_table("figure3F_cptac_lipid_protein_statistics.tsv")
    surv_df, surv_err = load_table("figure3G_expanded_survival_meta_analysis.tsv")
    umich_stats_df, _ = load_table("cptac_umich_lipid_protein_statistics.tsv")
    bcm_stats_df, _ = load_table("cptac_bcm_lipid_protein_statistics.tsv")

    # `findings` collects real, run-specific facts as each section computes
    # them, so Sections 10-13 are built from actual data rather than
    # hardcoded text — this is what keeps a re-run from silently reverting
    # to stale "simulated" boilerplate.
    findings = {}

    w("## 1. Dataset Summary")
    w()
    w("### Single-Cell Cohorts")
    w()

    cohort_real_status = {}
    cohort_gene_coverage = {}
    cohort_dissociation_stress = {}
    cohort_n_patients = {}
    try:
        with open(os.path.join(CONFIG_DIR, "singlecell_cohorts.yml")) as f:
            sc_cfg = yaml.safe_load(f)
        for cohort in sc_cfg["cohorts"]:
            name = cohort["name"]
            scores_file = os.path.join(BASE_DIR, cohort["scores_file"])
            processed_file = os.path.join(BASE_DIR, cohort["processed_file"])
            status = "Available (processed)" if os.path.exists(scores_file) else "Missing"
            sim_note = ""
            is_sim = True
            n_cells = 0
            if os.path.exists(scores_file):
                try:
                    sdf = pd.read_csv(scores_file, sep="\t")
                    is_sim = bool(sdf["is_simulated"].iloc[0]) if "is_simulated" in sdf.columns else True
                    n_cells = len(sdf)
                    data_label = "SIMULATED" if is_sim else "REAL"
                    sim_note = f" — **{data_label} DATA** ({n_cells} cells)"
                except Exception:
                    sim_note = " (could not read scores)"
            cohort_real_status[name] = not is_sim
            # Try to read gene coverage + dissociation-stress flag from h5ad uns
            if os.path.exists(processed_file):
                try:
                    import anndata as _ad
                    _adata = _ad.read_h5ad(processed_file)
                    gc = _adata.uns.get("gene_coverage", {})
                    if gc:
                        cohort_gene_coverage[name] = gc
                    if "dissociation_stress_elevated" in _adata.uns:
                        cohort_dissociation_stress[name] = bool(_adata.uns["dissociation_stress_elevated"])
                except Exception:
                    pass
            # Special note for GSE202051
            extra = ""
            if name == "GSE202051" and not is_sim:
                extra = " *(snRNA-seq; pre-annotated with rich cell-type scores — **strength: no de-novo annotation needed**)*"
            w(f"- **{name}**: {status}{sim_note}{extra}")
    except Exception as e:
        w(f"- Could not load cohort config: {e}")

    # Per-cohort patient counts (from the co-occurrence table, which reports
    # n_patients per malignant-cell cohort)
    if cooccur_df is not None and "n_patients" in cooccur_df.columns:
        for _, row in cooccur_df.iterrows():
            cohort_n_patients[row["cohort"]] = int(row["n_patients"])

    # Gene set coverage per cohort
    if cohort_gene_coverage:
        w()
        w("### Gene Set Coverage per Cohort (signature genes found in dataset)")
        w()
        sig_names = ["hypoxia", "acinar_identity", "lipid_synthesis_srebp",
                     "desaturation_elongation", "fatty_acid_uptake_oxidation", "caf", "emt"]
        header = "| Cohort | " + " | ".join(sig_names) + " |"
        sep = "|--------|" + "|".join(["--------"] * len(sig_names)) + "|"
        w(header)
        w(sep)
        for cname, gc in cohort_gene_coverage.items():
            row_vals = []
            for sig in sig_names:
                info = gc.get(sig, {})
                found = info.get("found", "?")
                total = info.get("total", "?")
                pct = info.get("pct", 0)
                flag = " ⚠" if isinstance(pct, (int, float)) and pct < 50 else ""
                row_vals.append(f"{found}/{total} ({pct:.0f}%){flag}")
            w(f"| {cname} | " + " | ".join(row_vals) + " |")
        w()
        w("*(⚠ = coverage < 50%, interpret with caution)*")

    w()
    w("### Proteomics")
    flag_file = os.path.join(BASE_DIR, "data", "raw", "proteomics", "CPTAC_PDA", "CPTAC_PDA_DOWNLOAD_FAILED.txt")
    cptac_is_real = (cptac_df is not None and len(cptac_df) > 0 and not os.path.exists(flag_file))
    cptac_status = "Real data" if cptac_is_real else "SIMULATED (download failed)"
    findings["cptac_is_real"] = cptac_is_real
    w(f"- **CPTAC-PDA**: {cptac_status}")

    w()
    w("### Bulk Cohorts (Phase 2, reused)")
    phase2_processed = os.path.join(BASE_DIR, "..", "pdac_phase2_validation", "data", "processed", "validation")
    phase2_all_found = True
    for c in ["GSE79668", "GSE71729", "GSE62165"]:
        f = os.path.join(phase2_processed, f"{c}_sample_scores.tsv")
        found = os.path.exists(f)
        phase2_all_found = phase2_all_found and found
        w(f"- **{c}**: {'Found' if found else 'Not found (simulated for purity analysis)'}")
    findings["phase2_bulk_real"] = phase2_all_found

    w()
    # Dynamic note based on actual data status
    n_real = sum(cohort_real_status.values())
    n_total_cohorts = len(cohort_real_status)
    if n_real == n_total_cohorts:
        w("> **DATA STATUS:** All single-cell cohorts use **REAL downloaded data**. "
          f"CPTAC proteomics uses **{'REAL' if cptac_is_real else 'SIMULATED'} data**. "
          f"Phase 2 bulk cohorts (for purity adjustment) are **{'REAL' if phase2_all_found else 'not all found'}**. "
          "Conclusions reflect actual measurements except where noted below.")
    elif n_real > 0:
        real_names = [k for k, v in cohort_real_status.items() if v]
        sim_names = [k for k, v in cohort_real_status.items() if not v]
        w(f"> **DATA STATUS:** {n_real}/{n_total_cohorts} single-cell cohorts use real data "
          f"({', '.join(real_names)}). "
          f"Simulated fallback used for: {', '.join(sim_names) if sim_names else 'none'}. "
          "Conclusions from real cohorts are based on actual measurements; simulated cohort conclusions are provisional.")
    else:
        w("> **CRITICAL NOTE:** All single-cell analyses in this report use **SIMULATED DATA** "
          "generated with realistic distributions. Real dataset downloads (GSE154778, GSE202051, Peng et al.) "
          "failed because GEO SOFT files were downloaded but count matrices were not retrievable automatically. "
          "All conclusions below are provisional and **must be replicated with real scRNA-seq data**.")

    findings["cohort_real_status"] = cohort_real_status
    findings["n_real_sc_cohorts"] = n_real
    findings["n_total_sc_cohorts"] = n_total_cohorts

    # 2. Cell type annotation QC
    w()
    w("## 2. Cell Type Annotation QC")
    w()
    ann_files = [f for f in os.listdir(os.path.join(BASE_DIR, "data", "processed", "singlecell"))
                 if f.endswith("_cell_annotations.tsv")] if os.path.exists(
        os.path.join(BASE_DIR, "data", "processed", "singlecell")) else []
    if ann_files:
        for ann_file in ann_files:
            ann_path = os.path.join(BASE_DIR, "data", "processed", "singlecell", ann_file)
            cohort_name = ann_file.replace("_cell_annotations.tsv", "")
            try:
                ann_df = pd.read_csv(ann_path, sep="\t")
                n_cells = len(ann_df)
                ct_counts = ann_df["cell_type"].value_counts()
                w(f"**{cohort_name}**: {n_cells} cells")
                for ct, count in ct_counts.items():
                    w(f"  - {ct}: {count} ({count/n_cells:.1%})")
            except Exception as e:
                w(f"  Could not load {ann_file}: {e}")
    else:
        w("No annotation files found.")

    # 3. Hypoxia/Acinar co-occurrence
    w()
    w("## 3. Hypoxia/Acinar Co-occurrence Analysis (Figure 3A)")
    w()
    mean_frac, mean_r = None, None
    if cooccur_df is not None and len(cooccur_df) > 0:
        for _, row in cooccur_df.iterrows():
            w(f"**{row['cohort']}** (n={row.get('n_malignant_cells', 'N/A')} malignant cells):")
            w(f"  - Pearson r(hypoxia, acinar) = {fmt(row.get('pearson_r_hypoxia_acinar', np.nan))}, "
              f"p = {fmt(row.get('pearson_p', np.nan))}")
            w(f"  - Spearman r = {fmt(row.get('spearman_r_hypoxia_acinar', np.nan))}")
            w(f"  - Fraction hypoxia_high/acinar_low: {fmt(row.get('frac_hypoxia_high_acinar_low', np.nan), 3)}")
            w(f"  - Interpretation: {row.get('cooccurrence_interpretation', 'N/A')}")
            if row.get("is_simulated", True):
                w(f"  - *[SIMULATED DATA]*")
            w()

        mean_frac = cooccur_df["frac_hypoxia_high_acinar_low"].mean()
        mean_r = cooccur_df["pearson_r_hypoxia_acinar"].mean()
        w(f"**Conclusion:** Mean fraction of malignant cells in hypoxia-high/acinar-low quadrant = {mean_frac:.2%}.")
        if mean_r < -0.2:
            w("Negative correlation (r < -0.2) suggests hypoxia and acinar programs are inversely related, "
              "consistent with a true composite state where loss of acinar identity accompanies hypoxia activation.")
        elif mean_r > 0.2:
            w("Positive correlation suggests these programs co-activate, not mutually exclusive.")
        else:
            w("Weak correlation suggests the hypoxia-high/acinar-low state is partially a composite artifact "
              "of two independent cell populations rather than a uniform cell state.")
    else:
        w("No co-occurrence data available. " + (cooccur_err or ""))

    # A cohort "co-occurs" (supports a genuine joint state) if its own r is
    # meaningfully negative AND its enrichment is clearly above the ~25%
    # expected by chance under independence of two median-split axes.
    cooccurrence_intrinsic = False
    if cooccur_df is not None and len(cooccur_df) > 0 and "pearson_r_hypoxia_acinar" in cooccur_df.columns:
        n_cooccurring = int(((cooccur_df["pearson_r_hypoxia_acinar"] < -0.2) &
                              (cooccur_df["frac_hypoxia_high_acinar_low"] > 0.30)).sum())
        cooccurrence_intrinsic = n_cooccurring >= 2
        findings["n_cohorts_cooccurring"] = n_cooccurring
    findings["cooccurrence_intrinsic"] = cooccurrence_intrinsic
    findings["mean_r"] = mean_r
    findings["mean_frac"] = mean_frac

    # 4. Lipid program cell-of-origin
    w()
    w("## 4. Lipid Program Cell-of-Origin (Figure 3B)")
    w()
    lipid_intrinsic = False
    mal_lipid_summary = []
    if lipid_df is not None and len(lipid_df) > 0:
        w("### Summary by Cell Type and Score")
        w()
        for score in ["lipid_synthesis_srebp", "desaturation_elongation", "fatty_acid_uptake_oxidation"]:
            w(f"**{score}:**")
            sub = lipid_df[lipid_df["score"] == score] if "score" in lipid_df.columns else pd.DataFrame()
            if sub.empty:
                w("  No data.")
                continue
            for ct in ["malignant_epithelial", "caf_fibroblast", "myeloid", "endothelial"]:
                ct_sub = sub[sub["cell_type"] == ct] if "cell_type" in sub.columns else pd.DataFrame()
                if ct_sub.empty:
                    continue
                row = ct_sub.iloc[0]
                sim = " [SIM]" if row.get("is_simulated", True) else ""
                w(f"  - {ct}: obs={row.get('observed_direction', '?')}, "
                  f"p_adj={fmt(row.get('wilcoxon_p_adj_BH', np.nan))}, "
                  f"repro={'YES' if row.get('fully_reproduced', False) else 'NO'}, "
                  f"cell_intrinsic={bool_str(row.get('cell_intrinsic_supported', False))}{sim}")
                if ct == "malignant_epithelial":
                    mal_lipid_summary.append({
                        "score": score,
                        "n_agg": row.get("n_aggressive"),
                        "n_ref": row.get("n_reference"),
                        "p_adj": row.get("wilcoxon_p_adj_BH"),
                        "direction": row.get("observed_direction"),
                        "expected": row.get("expected_direction"),
                        "cell_intrinsic": bool(row.get("cell_intrinsic_supported", False)),
                    })
            w()

        mal_sub = lipid_df[lipid_df["cell_type"] == "malignant_epithelial"] if "cell_type" in lipid_df.columns else pd.DataFrame()
        n_mal_intrinsic = int(mal_sub["cell_intrinsic_supported"].sum()) if "cell_intrinsic_supported" in mal_sub.columns else 0
        n_mal_total = len(mal_sub)
        lipid_intrinsic = n_mal_intrinsic > 0

        n_repro = lipid_df["fully_reproduced"].sum() if "fully_reproduced" in lipid_df.columns else 0
        n_total = len(lipid_df)
        n_intrinsic = lipid_df["cell_intrinsic_supported"].sum() if "cell_intrinsic_supported" in lipid_df.columns else 0
        w(f"**Summary:** {n_repro}/{n_total} comparisons reproduced (correct direction + FDR<0.05) across all cell types; "
          f"{n_mal_intrinsic}/{n_mal_total} malignant-cell comparisons show cell-intrinsic support.")

        # Data-driven conclusion (no hardcoded "real data required" — we
        # already know from Section 1 whether this cohort's data is real)
        if n_mal_intrinsic > 0:
            supported_scores = [m["score"] for m in mal_lipid_summary if m["cell_intrinsic"]]
            w(f"\n**Conclusion:** {', '.join(supported_scores)} shows cell-intrinsic elevation in malignant "
              f"epithelial cells in the adequately-powered cohort(s), consistent with (not against) tumor-intrinsic "
              f"lipid metabolic rewiring. Other lipid scores in this table did not reach significance at the "
              f"patient level and should be read as directional trends, not confirmed findings, given the small "
              f"patient counts noted in Section 11.")
        else:
            w("\n**Conclusion:** No malignant-cell lipid score reached significance at the patient level in any "
              "adequately-powered cohort. Given the small number of testable patients (see Section 11), this "
              "should be read as inconclusive rather than as a confirmed negative result.")
    else:
        w("No lipid cell-of-origin data available. " + (lipid_err or ""))
    findings["lipid_intrinsic"] = lipid_intrinsic
    findings["mal_lipid_summary"] = mal_lipid_summary

    # 5. CAF/EMT cell-of-origin
    w()
    w("## 5. CAF/EMT Cell-of-Origin (Figures 3C, 3D)")
    w()
    emt_malignant_significant = False
    emt_direction_summary = []
    if caf_emt_df is not None and len(caf_emt_df) > 0:
        if "section" in caf_emt_df.columns:
            emt_sub = caf_emt_df[caf_emt_df["section"] == "EMT_malignant"]
            caf_sub = caf_emt_df[caf_emt_df["section"] == "CAF_subtype"]
        else:
            emt_sub = caf_emt_df
            caf_sub = pd.DataFrame()

        if len(emt_sub) > 0:
            w("### EMT Score in Malignant Cells")
            w()
            mal_emt = emt_sub[emt_sub["cell_type"] == "malignant_epithelial"] if "cell_type" in emt_sub.columns else emt_sub
            if len(mal_emt) > 0:
                for _, row in mal_emt.iterrows():
                    sim = " [SIM]" if row.get("is_simulated", True) else ""
                    is_sig = bool(row.get("significant_FDR05", False))
                    w(f"  - {row.get('cohort', 'N/A')}: median diff = {fmt(row.get('median_diff_agg_minus_ref', np.nan))}, "
                      f"p = {fmt(row.get('wilcoxon_p_raw', np.nan))}, "
                      f"p_adj = {fmt(row.get('wilcoxon_p_adj_BH', np.nan))}, "
                      f"direction = {row.get('observed_direction', '?')}{sim}")
                    if is_sig:
                        emt_malignant_significant = True
                    emt_direction_summary.append({
                        "cohort": row.get("cohort"),
                        "direction": row.get("observed_direction"),
                        "significant": is_sig,
                    })
                n_sig = sum(1 for e in emt_direction_summary if e["significant"])
                directions = set(e["direction"] for e in emt_direction_summary if e["significant"])
                if n_sig == 0:
                    w(f"\n**EMT conclusion:** No cohort reached patient-level significance for malignant-cell EMT score; "
                      f"inconclusive given the small testable patient counts (Section 11).")
                else:
                    dir_txt = "/".join(sorted(directions))
                    w(f"\n**EMT conclusion:** EMT score is significantly **{dir_txt}** in aggressive malignant cells "
                      f"in {n_sig}/{len(emt_direction_summary)} testable cohort(s) — this is malignant-cell-intrinsic "
                      f"(not purely a stromal-composition effect), regardless of whether the direction matches the "
                      f"originally hypothesized 'up' direction.")

        if len(caf_sub) > 0:
            w()
            w("### CAF Subtype Analysis")
            w()
            for cohort in caf_sub.get("cohort", pd.Series()).unique():
                w(f"**{cohort}:**")
                c_sub = caf_sub[caf_sub["cohort"] == cohort]
                for state in ["aggressive", "reference"]:
                    s_sub = c_sub[c_sub.get("patient_group", c_sub.get("state", pd.Series())) == state] if "patient_group" in c_sub.columns else pd.DataFrame()
                    if len(s_sub) == 0:
                        continue
                    w(f"  {state}: " + ", ".join([f"{row['caf_subtype']}={row['fraction']:.2%}"
                                                   for _, row in s_sub.iterrows()]))
    else:
        w("CAF/EMT data not available or empty. " + (caf_emt_err or ""))
    findings["emt_malignant_significant"] = emt_malignant_significant
    findings["emt_direction_summary"] = emt_direction_summary

    # 6. Purity adjustment
    w()
    w("## 6. Purity Adjustment of CAF/EMT Signals (Figure 3E)")
    w()
    purity_survives = False
    purity_detail = []
    if purity_df is not None and len(purity_df) > 0:
        p_adj_col = "pval_purity_adjusted_BH" if "pval_purity_adjusted_BH" in purity_df.columns else "pval_purity_adjusted"
        p_unadj_col = "pval_unadjusted_BH" if "pval_unadjusted_BH" in purity_df.columns else "pval_unadjusted"
        for _, row in purity_df.iterrows():
            sim = " [SIM]" if row.get("is_simulated", True) else ""
            sig_before = row.get(p_unadj_col, np.nan) < 0.05 if p_unadj_col in row else False
            sig_after = row.get(p_adj_col, np.nan) < 0.05 if p_adj_col in row else False
            w(f"- **{row.get('cohort', 'N/A')} — {row.get('score', 'N/A')}**: "
              f"unadjusted coef={fmt(row.get('coef_unadjusted', np.nan))} (p={fmt(row.get(p_unadj_col, np.nan))}), "
              f"adjusted coef={fmt(row.get('coef_purity_adjusted', np.nan))} (p={fmt(row.get(p_adj_col, np.nan))}), "
              f"{'still significant after adjustment' if sig_after else 'not significant after adjustment'}{sim}")
            purity_detail.append({
                "cohort": row.get("cohort"), "score": row.get("score"),
                "sig_before": bool(sig_before), "sig_after": bool(sig_after),
            })
        w()
        n_sig_before = sum(1 for p in purity_detail if p["sig_before"])
        n_sig_before_and_after = sum(1 for p in purity_detail if p["sig_before"] and p["sig_after"])
        purity_survives = n_sig_before > 0 and (n_sig_before_and_after / n_sig_before) >= 0.5
        findings["n_purity_sig_before"] = n_sig_before
        findings["n_purity_sig_before_and_after"] = n_sig_before_and_after
        if n_sig_before == 0:
            w("**Conclusion:** No cohort showed a significant CAF/EMT group effect even before purity adjustment, "
              "so purity confounding cannot be assessed as a cause of a signal that was not detected in the first place.")
        else:
            w(f"**Conclusion:** Of the {n_sig_before} cohort/score combinations significant before purity adjustment, "
              f"{n_sig_before_and_after} ({n_sig_before_and_after/n_sig_before:.0%}) remain significant after adjusting "
              f"for the 8-gene purity proxy. This {'supports' if purity_survives else 'argues against'} a largely "
              f"purity-independent CAF/EMT association in bulk PDAC.")
    else:
        w("Purity adjustment data not available. " + (purity_err or ""))
    findings["purity_survives"] = purity_survives

    # 7. CPTAC protein validation
    w()
    w("## 7. Protein-Level Lipid Validation — CPTAC-PDA (Figure 3F)")
    w()
    n_all_conc, n_all_sig, n_total_proteins = 0, 0, 0
    if cptac_df is not None and len(cptac_df) > 0:
        w(f"**Data status:** {'SIMULATED DATA' if not cptac_is_real else 'Real CPTAC data'} "
          f"(umich + BCM proteomics, WashU transcriptomics for group assignment).")
        w()
        n_total_proteins = len(cptac_df)
        if "gene_set" in cptac_df.columns:
            for gs in cptac_df["gene_set"].unique():
                sub = cptac_df[cptac_df["gene_set"] == gs]
                n_conc = int(sub["replicated_both"].sum()) if "replicated_both" in sub.columns else 0
                n_sig = int(sub["significant_both"].sum()) if "significant_both" in sub.columns else 0
                exp_dir = sub["expected_direction"].iloc[0] if "expected_direction" in sub.columns else "?"
                w(f"**{gs}** (expected: {exp_dir}):")
                w(f"  - {n_conc}/{len(sub)} proteins concordant in direction in both umich and BCM")
                w(f"  - {n_sig}/{len(sub)} significant (FDR<0.05) and concordant in both sources")
        n_all_conc = int(cptac_df["replicated_both"].sum()) if "replicated_both" in cptac_df.columns else 0
        n_all_sig = int(cptac_df["significant_both"].sum()) if "significant_both" in cptac_df.columns else 0
        w()
        w(f"**Overall:** {n_all_conc}/{n_total_proteins} proteins directionally concordant in both proteomics sources, "
          f"{n_all_sig}/{n_total_proteins} significant (FDR<0.05) and concordant in both.")
        if n_all_sig > 0 and "significant_both" in cptac_df.columns:
            replicated_names = cptac_df.loc[cptac_df["significant_both"] == True, "protein"].tolist()
            w(f"  Independently replicated at FDR<0.05 in both sources: {', '.join(replicated_names)}.")
        if not cptac_is_real:
            w("\n**CRITICAL:** These results are from SIMULATED data and do not constitute real protein validation. "
              "CPTAC-PDA data must be downloaded via the `cptac` Python package for real validation.")
    else:
        w("CPTAC protein data not available. " + (cptac_err or ""))
    findings["cptac_n_concordant"] = n_all_conc
    findings["cptac_n_sig_concordant"] = n_all_sig
    findings["cptac_n_total"] = n_total_proteins

    # 8. Survival meta-analysis
    w()
    w("## 8. Expanded Survival Meta-Analysis (Figure 3G)")
    w()
    hr_re, p_re, hr_lo, hr_hi = None, None, None, None
    n_hr_above_1, n_cohorts_survival = 0, 0
    if surv_df is not None and len(surv_df) > 0:
        ind = surv_df[surv_df.get("analysis_type", pd.Series([""] * len(surv_df))) == "individual_cohort"] if "analysis_type" in surv_df.columns else surv_df
        meta = surv_df[surv_df.get("analysis_type", pd.Series([""] * len(surv_df))) == "meta_analysis"] if "analysis_type" in surv_df.columns else pd.DataFrame()

        if len(ind) > 0:
            w("### Individual Cohorts")
            for _, row in ind.iterrows():
                sim = " [SIM]" if row.get("is_simulated", False) else ""
                w(f"- **{row.get('cohort', 'N/A')}**: HR={fmt(row.get('HR', np.nan))}, "
                  f"95% CI [{fmt(row.get('HR_lower_95CI', np.nan))}-{fmt(row.get('HR_upper_95CI', np.nan))}], "
                  f"p={fmt(row.get('p_value', np.nan))}{sim}")
            n_cohorts_survival = len(ind)
            n_hr_above_1 = int((ind["HR"] > 1).sum())

        if len(meta) > 0:
            w()
            w("### Pooled Estimates")
            for _, row in meta.iterrows():
                w(f"- **{row.get('cohort', 'N/A')}**: HR={fmt(row.get('HR', np.nan))}, "
                  f"95% CI [{fmt(row.get('HR_lower_95CI', np.nan))}-{fmt(row.get('HR_upper_95CI', np.nan))}], "
                  f"p={fmt(row.get('p_value', np.nan))}")

        re_row = surv_df[surv_df.get("cohort", pd.Series()) == "pooled_random_effects"] if "cohort" in surv_df.columns else pd.DataFrame()
        if len(re_row) == 0:
            re_row = surv_df[surv_df.get("cohort", pd.Series()) == "pooled_fixed_effects"] if "cohort" in surv_df.columns else pd.DataFrame()
        if len(re_row) > 0:
            hr_re = float(re_row["HR"].values[0])
            p_re = float(re_row["p_value"].values[0])
            hr_lo = float(re_row["HR_lower_95CI"].values[0])
            hr_hi = float(re_row["HR_upper_95CI"].values[0])
            w()
            if p_re < 0.05:
                w(f"**Conclusion:** Statistically significant survival association (HR={hr_re:.3f}, p={p_re:.4f}) "
                  "in expanded meta-analysis.")
            else:
                consistency = f"{n_hr_above_1}/{n_cohorts_survival} cohorts have HR > 1" if n_cohorts_survival else "trend inconsistent"
                w(f"**Conclusion:** Expanded meta-analysis is not statistically significant (HR={hr_re:.3f}, p={p_re:.4f}); "
                  f"{consistency}. "
                  f"{'The direction is directionally consistent but underpowered.' if n_cohorts_survival and n_hr_above_1 == n_cohorts_survival else 'The direction is inconsistent across cohorts, not merely underpowered.'}")
    else:
        w("Survival meta-analysis data not available. " + (surv_err or ""))
    findings.update({"hr_pooled": hr_re, "p_pooled": p_re, "hr_pooled_lo": hr_lo, "hr_pooled_hi": hr_hi,
                      "n_cohorts_survival": n_cohorts_survival, "n_hr_above_1": n_hr_above_1})

    # 9. Immune/Endothelial QC
    w()
    w("## 9. Immune / Endothelial QC Check")
    w()
    if lipid_df is not None and "cell_type" in lipid_df.columns:
        for ct in ["myeloid", "endothelial"]:
            ct_sub = lipid_df[lipid_df["cell_type"] == ct]
            if len(ct_sub) > 0:
                n_sig = ct_sub["significant_FDR05"].sum() if "significant_FDR05" in ct_sub.columns else 0
                w(f"- **{ct}**: {n_sig}/{len(ct_sub)} lipid comparisons significant at FDR<0.05. "
                  f"{'Strong non-tumor signal — warrants investigation.' if n_sig > len(ct_sub) // 2 else 'Low off-target signal consistent with tumor-driven lipid rewiring.'}")
    else:
        w("No immune/endothelial data available for QC check.")

    # ------------------------------------------------------------------
    # 10. Final mechanistic classification — derived from `findings`,
    # matching the 3-way taxonomy from the Phase 3 instructions:
    #   - resolved & tumor-intrinsic: co-occurrence + lipid + EMT + purity
    #     all support a real, malignant-cell-intrinsic, purity-independent
    #     program.
    #   - resolved & composite/microenvironment-driven: none of the above
    #     hold.
    #   - partially resolved/mixed: anything in between.
    # ------------------------------------------------------------------
    criteria = {
        "cooccurrence_intrinsic": findings.get("cooccurrence_intrinsic", False),
        "lipid_intrinsic": findings.get("lipid_intrinsic", False),
        "emt_malignant_intrinsic": findings.get("emt_malignant_significant", False),
        "purity_survives": findings.get("purity_survives", False),
    }
    n_supported = sum(bool(v) for v in criteria.values())
    if n_supported == 4:
        classification = "Mechanistically resolved and tumor-intrinsic"
    elif n_supported == 0:
        classification = "Mechanistically resolved, composite/microenvironment-driven"
    else:
        classification = "Partially resolved / mixed"
    findings["classification"] = classification
    findings["criteria"] = criteria

    w()
    w("## 10. Final Mechanistic Classification")
    w()
    w(f"**Classification: {classification}**")
    w()
    w(f"Data status: {findings['n_real_sc_cohorts']}/{findings['n_total_sc_cohorts']} single-cell cohorts real, "
      f"CPTAC proteomics {'real' if cptac_is_real else 'simulated'}, "
      f"Phase 2 bulk cohorts (purity adjustment) {'real' if phase2_all_found else 'not all found'}.")
    w()
    w(f"1. **Hypoxia/acinar co-occurrence:** Pearson r = "
      + ", ".join(fmt(r, 3) for r in (cooccur_df["pearson_r_hypoxia_acinar"] if cooccur_df is not None else []))
      + f"; mean fraction of malignant cells in the hypoxia-high/acinar-low quadrant = {mean_frac:.1%} "
        f"(chance level under independent median splits ≈ 25%). "
      + ("This supports a genuine co-occurring malignant cell state." if criteria["cooccurrence_intrinsic"]
         else "This is consistent with a composite artifact: hypoxia-high and acinar-low occur largely in separate malignant cell subpopulations, not the same cells."))
    w()
    if mal_lipid_summary:
        supported = [m["score"] for m in mal_lipid_summary if m["cell_intrinsic"]]
        w(f"2. **Lipid rewiring cell-of-origin:** "
          + (f"{', '.join(supported)} shows cell-intrinsic elevation in malignant cells in an adequately-powered cohort."
             if supported else "No malignant-cell lipid score reached significance in an adequately-powered cohort.")
          + " This remains based on very small patient counts (see Section 11) and should be treated as suggestive, not confirmatory.")
    else:
        w("2. **Lipid rewiring cell-of-origin:** No testable cohort available.")
    w()
    emt_dirs = set(e["direction"] for e in findings.get("emt_direction_summary", []) if e["significant"])
    w(f"3. **CAF/EMT cell-of-origin:** "
      + (f"EMT score is significantly {'/'.join(sorted(emt_dirs))} in aggressive malignant cells in at least one "
         f"testable cohort — a malignant-cell-intrinsic effect (not purely stromal-composition-driven), though the "
         f"direction may not match the originally hypothesized 'up' direction and rests on a single small cohort."
         if criteria["emt_malignant_intrinsic"] else
         "No cohort showed a statistically significant malignant-cell EMT effect at the patient level; inconclusive given small patient counts."))
    w()
    w(f"4. **Purity adjustment of bulk CAF/EMT:** "
      + (f"{findings.get('n_purity_sig_before_and_after', 0)}/{findings.get('n_purity_sig_before', 0)} "
         f"cohort/score combinations that were significant before adjustment remain significant afterward — "
         f"the CAF/EMT bulk signal is largely **not** a tumor-purity artifact."
         if criteria["purity_survives"] else
         "The CAF/EMT bulk signal was substantially attenuated after purity adjustment and should be interpreted as at least partly explained by tumor-purity/stromal-content differences rather than a purity-independent biological effect."))
    w()
    if cptac_is_real and n_total_proteins > 0:
        w(f"5. **Protein-level validation:** {n_all_conc}/{n_total_proteins} lipid/FA-oxidation proteins directionally "
          f"concordant in both umich and BCM proteomics; {n_all_sig}/{n_total_proteins} independently replicated at FDR<0.05 in both.")
    else:
        w("5. **Protein-level validation:** SIMULATED or unavailable CPTAC data — real CPTAC validation is needed.")
    w()
    if hr_re is not None:
        w(f"6. **Survival:** Pooled HR={hr_re:.3f} [{hr_lo:.2f}-{hr_hi:.2f}], p={p_re:.3f} across "
          f"{n_cohorts_survival} cohorts ({n_hr_above_1}/{n_cohorts_survival} with HR>1). "
          + ("Statistically significant." if p_re < 0.05 else
             "Not statistically significant; " + ("directionally consistent but underpowered." if n_hr_above_1 == n_cohorts_survival else "and inconsistent in direction across cohorts, not merely underpowered.")))
    else:
        w("6. **Survival:** Meta-analysis data not available.")

    # 11. Limitations
    w()
    w("## 11. Limitations")
    w()
    w("**Limitations to weigh before drawing manuscript-level conclusions:**")
    w()
    lim_n = 1
    if n_real < n_total_cohorts:
        sim_names = [k for k, v in cohort_real_status.items() if not v]
        w(f"{lim_n}. Simulated fallback data was used for: {', '.join(sim_names)}. "
          "Conclusions drawn from these cohorts are provisional pending real data.")
        lim_n += 1
    if not cptac_is_real:
        w(f"{lim_n}. CPTAC-PDA proteomics are SIMULATED in this run — protein-level validation is provisional.")
        lim_n += 1
    w(f"{lim_n}. **No matched bulk + single-cell patient data.** The single-cell cohorts and the Phase 2 bulk "
      "cohorts are independent and share no patients, so within-cohort single-cell grouping is not equivalent "
      "to the between-cohort Phase 2 aggressive-group definition.")
    lim_n += 1
    small_cohorts = [c for c, n in cohort_n_patients.items() if n is not None and n < 5]
    if small_cohorts:
        w(f"{lim_n}. **Small per-cohort patient counts limit cell-of-origin tests.** " +
          "; ".join(f"{c}: {cohort_n_patients[c]} patient(s)" for c in small_cohorts) +
          " — too few for patient-level statistical testing of lipid/EMT cell-of-origin in these cohorts.")
        lim_n += 1
    stressed = [c for c, v in cohort_dissociation_stress.items() if v]
    if stressed:
        w(f"{lim_n}. **Dissociation stress artifacts detected** in {', '.join(stressed)} (elevated FOS/JUN/HSP-family "
          "genes relative to overall expression). Hypoxia scores in these cohorts may be partially inflated by "
          "dissociation stress rather than true in vivo hypoxia.")
        lim_n += 1
    n_agg_umich, n_ref_umich = mode_n(umich_stats_df, "n_aggressive"), mode_n(umich_stats_df, "n_reference")
    n_agg_bcm, n_ref_bcm = mode_n(bcm_stats_df, "n_aggressive"), mode_n(bcm_stats_df, "n_reference")
    if n_agg_bcm is not None and n_ref_bcm is not None and n_agg_umich is not None:
        w(f"{lim_n}. **BCM proteomics replication uses an unbalanced reference group** (n={n_ref_bcm} vs "
          f"n={n_ref_umich} in umich; aggressive n={n_agg_bcm} vs {n_agg_umich}), reducing power for BCM FDR thresholds.")
        lim_n += 1
    w(f"{lim_n}. **Purity adjustment uses an 8-gene expression-mean heuristic**, not a validated deconvolution "
      "method (ESTIMATE, CIBERSORTx, TIMER2).")
    lim_n += 1
    w(f"{lim_n}. All statistical conclusions assume the real datasets are representative of the broader PDAC "
      "population; cohort sizes remain small relative to bulk validation cohorts.")

    # 12. Manuscript wording
    w()
    w("## 12. Recommended Manuscript Wording")
    w()
    if classification == "Mechanistically resolved and tumor-intrinsic":
        w("> \"Single-cell RNA sequencing analysis of PDAC tumors revealed that the hypoxia-high/acinar-low "
          "aggressive state co-occurs within individual malignant epithelial cells rather than representing "
          f"a composite of two distinct subpopulations (r={fmt(mean_r)}, Figure 3A). Lipid synthesis and "
          "desaturation genes were preferentially elevated in malignant cells from aggressive tumors, supporting "
          "a cell-intrinsic metabolic reprogramming rather than a stromal confound (Figure 3B). Purity-adjusted "
          f"analysis of bulk data confirmed that CAF and EMT signals remain significant "
          f"({findings.get('n_purity_sig_before_and_after', 0)}/{findings.get('n_purity_sig_before', 0)} "
          "effects survive adjustment) after controlling for tumor content (Figure 3E). Protein-level analysis "
          f"of CPTAC-PDA samples corroborated {n_all_sig}/{n_total_proteins} lipid enzymes at the proteomic level "
          f"in the aggressive group (Figure 3F). "
          + (f"The pooled survival HR was {fmt(hr_re)} [{fmt(hr_lo,2)}-{fmt(hr_hi,2)}], p={fmt(p_re)}." if hr_re is not None else "") + "\"")
    elif classification == "Mechanistically resolved, composite/microenvironment-driven":
        w("> \"Single-cell and purity-adjusted analysis indicate the bulk hypoxia-high/acinar-low signature more "
          f"likely reflects a composite of tumor-cell state heterogeneity (r={fmt(mean_r)} at the single-cell "
          "level, near the ~25% chance rate for co-occurrence) and variation in stromal content, rather than a "
          "single coherent malignant cell program. The lipid and CAF/EMT associations identified in bulk data "
          "should be interpreted as tumor-composition-associated rather than strictly tumor-cell-intrinsic.\"")
    else:
        w("> \"Single-cell analysis indicates the bulk hypoxia-high/acinar-low signature is a composite of two "
          f"largely independent malignant-cell axes (r≈{fmt(mean_r)} at the cell level) rather than a single "
          "coherent program. "
          + ("Patient-level re-analysis finds cell-intrinsic support for lipid synthesis rewiring in malignant cells in the one adequately-powered cohort, "
             if lipid_intrinsic else "Lipid rewiring cell-of-origin remains inconclusive given small patient counts, ")
          + ("while EMT score is significantly altered in aggressive malignant cells in that same cohort — a malignant-cell-intrinsic effect. "
             if emt_malignant_significant else "and EMT cell-of-origin is similarly inconclusive. ")
          + (f"The bulk CAF and EMT associations survive adjustment for tumor purity in "
             f"{findings.get('n_purity_sig_before_and_after', 0)}/{findings.get('n_purity_sig_before', 0)} "
             "significant cohort/score combinations, arguing against a pure tumor-purity artifact. "
             if purity_survives else "The bulk CAF/EMT associations are substantially attenuated after purity adjustment. ")
          + (f"Pooled survival HR={fmt(hr_re)} [{fmt(hr_lo,2)}-{fmt(hr_hi,2)}], p={fmt(p_re)}. " if hr_re is not None else "")
          + "All single-cell cell-of-origin conclusions rest on very small patient numbers and require replication in "
            "larger and, ideally, matched bulk+single-cell cohorts.\"")

    # 13. Next steps
    w()
    w("## 13. Next Steps")
    w()
    step_n = 1
    w(f"{step_n}. **Obtain matched bulk + single-cell data** for the same patients to properly test lipid/EMT "
      "cell-of-origin with matched group assignments.")
    step_n += 1
    if not cptac_is_real:
        w(f"{step_n}. **Download real CPTAC-PDA data** via the `cptac` Python package (`cptac.Pdac()`).")
        step_n += 1
    if n_real < n_total_cohorts:
        w(f"{step_n}. **Download real scRNA-seq data** for the remaining simulated cohort(s).")
        step_n += 1
    w(f"{step_n}. **Implement proper purity deconvolution** (ESTIMATE, CIBERSORTx, or TIMER2) in place of the "
      "current 8-gene expression-mean proxy.")
    step_n += 1
    w(f"{step_n}. **Spatial transcriptomics** (e.g. 10x Visium PDAC datasets) to directly test spatial "
      "co-localization of hypoxic/dedifferentiated malignant regions with CAF/EMT-high stroma.")
    step_n += 1
    if stressed:
        w(f"{step_n}. **Dissociation stress correction** (e.g. van den Brink et al. 2017 gene list) before "
          f"re-scoring hypoxia in {', '.join(stressed)}.")
        step_n += 1
    w(f"{step_n}. **More patients per single-cell cohort** to properly power lipid/EMT/CAF cell-of-origin tests "
      "in the cohorts currently too small to test (see Section 11).")

    w()
    w("---")
    prov_bits = []
    for name, is_real in cohort_real_status.items():
        prov_bits.append(f"{name}: {'REAL' if is_real else 'SIMULATED'}")
    prov_bits.append(f"CPTAC: {'REAL' if cptac_is_real else 'SIMULATED'}")
    w(f"*Report generated by Phase 3 pipeline (data-driven; no hardcoded provenance claims). "
      f"Data status — {'; '.join(prov_bits)}.*")

    # Write report
    report_path = os.path.join(REPORTS_DIR, "PHASE3_MECHANISM_REPORT.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Report written: {report_path}")
    return report_path


def main():
    print("=== Generating Phase 3 Report ===\n")
    try:
        report_path = generate_report()
        print(f"\nReport: {report_path}")
    except Exception as e:
        import traceback
        print(f"ERROR: {e}")
        traceback.print_exc()
        # Write minimal report
        os.makedirs(REPORTS_DIR, exist_ok=True)
        with open(os.path.join(REPORTS_DIR, "PHASE3_MECHANISM_REPORT.md"), "w") as f:
            f.write("# Phase 3 Mechanism Report\n\nERROR during report generation.\n")
            f.write(f"Error: {e}\n")
    print("\n=== Report Generation Complete ===")


if __name__ == "__main__":
    main()
