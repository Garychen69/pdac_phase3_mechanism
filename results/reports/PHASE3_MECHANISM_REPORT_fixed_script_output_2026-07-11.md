# Phase 3 Mechanistic Report: Single-Cell Resolution of Aggressive PDAC State

**Generated:** 2026-07-11 21:46:13

---

## 1. Dataset Summary

### Single-Cell Cohorts

- **GSE154778**: Available (processed) — **REAL DATA** (8000 cells)
- **GSE202051**: Available (processed) — **REAL DATA** (2607 cells) *(snRNA-seq; pre-annotated with rich cell-type scores — **strength: no de-novo annotation needed**)*
- **Peng_et_al**: Available (processed) — **REAL DATA** (43888 cells)

### Gene Set Coverage per Cohort (signature genes found in dataset)

| Cohort | hypoxia | acinar_identity | lipid_synthesis_srebp | desaturation_elongation | fatty_acid_uptake_oxidation | caf | emt |
|--------|--------|--------|--------|--------|--------|--------|--------|
| GSE154778 | 9/9 (100%) | 7/11 (64%) | 6/6 (100%) | 4/4 (100%) | 7/7 (100%) | 10/10 (100%) | 10/10 (100%) |
| GSE202051 | 9/9 (100%) | 7/11 (64%) | 6/6 (100%) | 4/4 (100%) | 6/7 (86%) | 8/10 (80%) | 7/10 (70%) |
| Peng_et_al | 9/9 (100%) | 10/11 (91%) | 6/6 (100%) | 4/4 (100%) | 7/7 (100%) | 10/10 (100%) | 10/10 (100%) |

*(⚠ = coverage < 50%, interpret with caution)*

### Proteomics
- **CPTAC-PDA**: Real data

### Bulk Cohorts (Phase 2, reused)
- **GSE79668**: Found
- **GSE71729**: Found
- **GSE62165**: Found

> **DATA STATUS:** All single-cell cohorts use **REAL downloaded data**. CPTAC proteomics uses **REAL data**. Phase 2 bulk cohorts (for purity adjustment) are **REAL**. Conclusions reflect actual measurements except where noted below.

## 2. Cell Type Annotation QC

**GSE154778**: 8000 cells
  - malignant_epithelial: 3710 (46.4%)
  - caf_fibroblast: 2101 (26.3%)
  - myeloid: 1257 (15.7%)
  - bcell_plasma: 367 (4.6%)
  - tcell_nk: 352 (4.4%)
  - endothelial: 213 (2.7%)
**GSE202051**: 2607 cells
  - malignant_epithelial: 2039 (78.2%)
  - ductal_normal: 475 (18.2%)
  - endothelial: 46 (1.8%)
  - acinar_normal: 25 (1.0%)
  - caf_fibroblast: 18 (0.7%)
  - myeloid: 4 (0.2%)
**Peng_et_al**: 43888 cells
  - myeloid: 12159 (27.7%)
  - malignant_epithelial: 11370 (25.9%)
  - tcell_nk: 10449 (23.8%)
  - endothelial: 3984 (9.1%)
  - bcell_plasma: 2587 (5.9%)
  - caf_fibroblast: 2478 (5.6%)
  - acinar_normal: 544 (1.2%)
  - ductal_normal: 317 (0.7%)

## 3. Hypoxia/Acinar Co-occurrence Analysis (Figure 3A)

**GSE154778** (n=3710 malignant cells):
  - Pearson r(hypoxia, acinar) = 0.060, p = 0.000
  - Spearman r = 0.075
  - Fraction hypoxia_high/acinar_low: 0.244
  - Interpretation: Largely separate populations (composite artifact)

**GSE202051** (n=2039 malignant cells):
  - Pearson r(hypoxia, acinar) = -0.030, p = 0.168
  - Spearman r = -0.008
  - Fraction hypoxia_high/acinar_low: 0.261
  - Interpretation: Largely separate populations (composite artifact)

**Peng_et_al** (n=11370 malignant cells):
  - Pearson r(hypoxia, acinar) = -0.016, p = 0.088
  - Spearman r = -0.018
  - Fraction hypoxia_high/acinar_low: 0.272
  - Interpretation: Largely separate populations (composite artifact)

**Conclusion:** Mean fraction of malignant cells in hypoxia-high/acinar-low quadrant = 25.91%.
Weak correlation suggests the hypoxia-high/acinar-low state is partially a composite artifact of two independent cell populations rather than a uniform cell state.

## 4. Lipid Program Cell-of-Origin (Figure 3B)

### Summary by Cell Type and Score

**lipid_synthesis_srebp:**
  - malignant_epithelial: obs=up, p_adj=0.018, repro=YES, cell_intrinsic=YES
  - caf_fibroblast: obs=up, p_adj=0.273, repro=NO, cell_intrinsic=NO
  - myeloid: obs=down, p_adj=0.855, repro=NO, cell_intrinsic=NO
  - endothelial: obs=up, p_adj=0.855, repro=NO, cell_intrinsic=NO

**desaturation_elongation:**
  - malignant_epithelial: obs=up, p_adj=0.584, repro=NO, cell_intrinsic=NO
  - caf_fibroblast: obs=up, p_adj=0.273, repro=NO, cell_intrinsic=NO
  - myeloid: obs=up, p_adj=0.855, repro=NO, cell_intrinsic=NO
  - endothelial: obs=up, p_adj=0.698, repro=NO, cell_intrinsic=NO

**fatty_acid_uptake_oxidation:**
  - malignant_epithelial: obs=up, p_adj=0.410, repro=NO, cell_intrinsic=NO
  - caf_fibroblast: obs=up, p_adj=0.273, repro=NO, cell_intrinsic=NO
  - myeloid: obs=up, p_adj=0.855, repro=NO, cell_intrinsic=NO
  - endothelial: obs=down, p_adj=0.698, repro=NO, cell_intrinsic=NO

**Summary:** 1/12 comparisons reproduced (correct direction + FDR<0.05) across all cell types; 1/3 malignant-cell comparisons show cell-intrinsic support.

**Conclusion:** lipid_synthesis_srebp shows cell-intrinsic elevation in malignant epithelial cells in the adequately-powered cohort(s), consistent with (not against) tumor-intrinsic lipid metabolic rewiring. Other lipid scores in this table did not reach significance at the patient level and should be read as directional trends, not confirmed findings, given the small patient counts noted in Section 11.

## 5. CAF/EMT Cell-of-Origin (Figures 3C, 3D)

### EMT Score in Malignant Cells

  - Peng_et_al: median diff = -0.135, p = 0.006, p_adj = 0.019, direction = down

**EMT conclusion:** EMT score is significantly **down** in aggressive malignant cells in 1/1 testable cohort(s) — this is malignant-cell-intrinsic (not purely a stromal-composition effect), regardless of whether the direction matches the originally hypothesized 'up' direction.

### CAF Subtype Analysis

**GSE154778:**
  aggressive: myCAF=6.55%, iCAF=2.92%, apCAF=0.54%
  reference: myCAF=13.53%, iCAF=5.84%, apCAF=0.63%
**Peng_et_al:**
  aggressive: myCAF=18.68%, iCAF=9.54%, apCAF=1.19%
  reference: myCAF=23.73%, iCAF=9.95%, apCAF=1.61%

## 6. Purity Adjustment of CAF/EMT Signals (Figure 3E)

- **GSE79668 — caf_score**: unadjusted coef=0.419 (p=0.152), adjusted coef=0.048 (p=0.802), not significant after adjustment
- **GSE79668 — emt_score**: unadjusted coef=0.395 (p=0.152), adjusted coef=0.047 (p=0.802), not significant after adjustment
- **GSE71729 — caf_score**: unadjusted coef=0.527 (p=0.001), adjusted coef=0.461 (p=0.000), still significant after adjustment
- **GSE71729 — emt_score**: unadjusted coef=0.493 (p=0.000), adjusted coef=0.446 (p=0.000), still significant after adjustment
- **GSE62165 — caf_score**: unadjusted coef=0.440 (p=0.011), adjusted coef=0.431 (p=0.000), still significant after adjustment
- **GSE62165 — emt_score**: unadjusted coef=0.560 (p=0.000), adjusted coef=0.552 (p=0.000), still significant after adjustment

**Conclusion:** Of the 4 cohort/score combinations significant before purity adjustment, 4 (100%) remain significant after adjusting for the 8-gene purity proxy. This supports a largely purity-independent CAF/EMT association in bulk PDAC.

## 7. Protein-Level Lipid Validation — CPTAC-PDA (Figure 3F)

**Data status:** Real CPTAC data (umich + BCM proteomics, WashU transcriptomics for group assignment).

**lipid_synthesis_srebp** (expected: up):
  - 4/6 proteins concordant in direction in both umich and BCM
  - 0/6 significant (FDR<0.05) and concordant in both sources
**fatty_acid_uptake_oxidation** (expected: down):
  - 1/6 proteins concordant in direction in both umich and BCM
  - 1/6 significant (FDR<0.05) and concordant in both sources
**desaturation_elongation** (expected: up):
  - 2/3 proteins concordant in direction in both umich and BCM
  - 0/3 significant (FDR<0.05) and concordant in both sources

**Overall:** 7/15 proteins directionally concordant in both proteomics sources, 1/15 significant (FDR<0.05) and concordant in both.
  Independently replicated at FDR<0.05 in both sources: ACADL.

## 8. Expanded Survival Meta-Analysis (Figure 3G)

### Individual Cohorts
- **GSE79668**: HR=1.232, 95% CI [0.629-2.413], p=0.543
- **GSE71729**: HR=1.079, 95% CI [0.661-1.761], p=0.761
- **GSE21501**: HR=0.927, 95% CI [0.519-1.655], p=0.798

### Pooled Estimates
- **pooled_fixed_effects**: HR=1.061, 95% CI [0.765-1.471], p=0.723
- **pooled_random_effects**: HR=1.061, 95% CI [0.765-1.471], p=0.723

**Conclusion:** Expanded meta-analysis is not statistically significant (HR=1.061, p=0.7231); 2/3 cohorts have HR > 1. The direction is inconsistent across cohorts, not merely underpowered.

## 9. Immune / Endothelial QC Check

- **myeloid**: 0/3 lipid comparisons significant at FDR<0.05. Low off-target signal consistent with tumor-driven lipid rewiring.
- **endothelial**: 0/3 lipid comparisons significant at FDR<0.05. Low off-target signal consistent with tumor-driven lipid rewiring.

## 10. Final Mechanistic Classification

**Classification: Partially resolved / mixed**

Data status: 3/3 single-cell cohorts real, CPTAC proteomics real, Phase 2 bulk cohorts (purity adjustment) real.

1. **Hypoxia/acinar co-occurrence:** Pearson r = 0.060, -0.030, -0.016; mean fraction of malignant cells in the hypoxia-high/acinar-low quadrant = 25.9% (chance level under independent median splits ≈ 25%). This is consistent with a composite artifact: hypoxia-high and acinar-low occur largely in separate malignant cell subpopulations, not the same cells.

2. **Lipid rewiring cell-of-origin:** lipid_synthesis_srebp shows cell-intrinsic elevation in malignant cells in an adequately-powered cohort. This remains based on very small patient counts (see Section 11) and should be treated as suggestive, not confirmatory.

3. **CAF/EMT cell-of-origin:** EMT score is significantly down in aggressive malignant cells in at least one testable cohort — a malignant-cell-intrinsic effect (not purely stromal-composition-driven), though the direction may not match the originally hypothesized 'up' direction and rests on a single small cohort.

4. **Purity adjustment of bulk CAF/EMT:** 4/4 cohort/score combinations that were significant before adjustment remain significant afterward — the CAF/EMT bulk signal is largely **not** a tumor-purity artifact.

5. **Protein-level validation:** 7/15 lipid/FA-oxidation proteins directionally concordant in both umich and BCM proteomics; 1/15 independently replicated at FDR<0.05 in both.

6. **Survival:** Pooled HR=1.061 [0.77-1.47], p=0.723 across 3 cohorts (2/3 with HR>1). Not statistically significant; and inconsistent in direction across cohorts, not merely underpowered.

## 11. Limitations

**Limitations to weigh before drawing manuscript-level conclusions:**

1. **No matched bulk + single-cell patient data.** The single-cell cohorts and the Phase 2 bulk cohorts are independent and share no patients, so within-cohort single-cell grouping is not equivalent to the between-cohort Phase 2 aggressive-group definition.
2. **Small per-cohort patient counts limit cell-of-origin tests.** GSE202051: 1 patient(s) — too few for patient-level statistical testing of lipid/EMT cell-of-origin in these cohorts.
3. **Dissociation stress artifacts detected** in GSE154778, Peng_et_al (elevated FOS/JUN/HSP-family genes relative to overall expression). Hypoxia scores in these cohorts may be partially inflated by dissociation stress rather than true in vivo hypoxia.
4. **BCM proteomics replication uses an unbalanced reference group** (n=23 vs n=46 in umich; aggressive n=41 vs 46), reducing power for BCM FDR thresholds.
5. **Purity adjustment uses an 8-gene expression-mean heuristic**, not a validated deconvolution method (ESTIMATE, CIBERSORTx, TIMER2).
6. All statistical conclusions assume the real datasets are representative of the broader PDAC population; cohort sizes remain small relative to bulk validation cohorts.

## 12. Recommended Manuscript Wording

> "Single-cell analysis indicates the bulk hypoxia-high/acinar-low signature is a composite of two largely independent malignant-cell axes (r≈0.004 at the cell level) rather than a single coherent program. Patient-level re-analysis finds cell-intrinsic support for lipid synthesis rewiring in malignant cells in the one adequately-powered cohort, while EMT score is significantly altered in aggressive malignant cells in that same cohort — a malignant-cell-intrinsic effect. The bulk CAF and EMT associations survive adjustment for tumor purity in 4/4 significant cohort/score combinations, arguing against a pure tumor-purity artifact. Pooled survival HR=1.061 [0.77-1.47], p=0.723. All single-cell cell-of-origin conclusions rest on very small patient numbers and require replication in larger and, ideally, matched bulk+single-cell cohorts."

## 13. Next Steps

1. **Obtain matched bulk + single-cell data** for the same patients to properly test lipid/EMT cell-of-origin with matched group assignments.
2. **Implement proper purity deconvolution** (ESTIMATE, CIBERSORTx, or TIMER2) in place of the current 8-gene expression-mean proxy.
3. **Spatial transcriptomics** (e.g. 10x Visium PDAC datasets) to directly test spatial co-localization of hypoxic/dedifferentiated malignant regions with CAF/EMT-high stroma.
4. **Dissociation stress correction** (e.g. van den Brink et al. 2017 gene list) before re-scoring hypoxia in GSE154778, Peng_et_al.
5. **More patients per single-cell cohort** to properly power lipid/EMT/CAF cell-of-origin tests in the cohorts currently too small to test (see Section 11).

---
*Report generated by Phase 3 pipeline (data-driven; no hardcoded provenance claims). Data status — GSE154778: REAL; GSE202051: REAL; Peng_et_al: REAL; CPTAC: REAL.*