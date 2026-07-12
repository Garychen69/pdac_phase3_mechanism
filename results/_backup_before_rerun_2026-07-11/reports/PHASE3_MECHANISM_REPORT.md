# Phase 3 Mechanistic Report: Single-Cell Resolution of Aggressive PDAC State

**Generated:** 2026-07-03 (updated 2026-07-04 with BCM replication; updated 2026-07-06 with patient-level pseudobulk re-analysis of lipid/CAF/EMT/immune-endothelial cell-of-origin, Sections 4/5/9; updated 2026-07-06 with a non-circular purity proxy, Section 6; updated 2026-07-06 with a real GSE21501 Cox fit replacing the hardcoded survival value, Section 8)

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
- **CPTAC-PDA**: **REAL DATA** — umich (145 tumor samples) + BCM (105 tumor samples), via `cptac.Pdac()`; WashU transcriptomics for group assignment

### Bulk Cohorts (Phase 2, reused)
- **GSE79668**: Found
- **GSE71729**: Found
- **GSE62165**: Found

> **DATA STATUS:** All single-cell cohorts use **REAL downloaded data**. CPTAC proteomics uses **REAL data** (umich + BCM). Purity adjustment uses **REAL Phase 2 bulk expression data** with a purity proxy verified disjoint from the CAF/EMT signatures. Survival meta-analysis uses a **real Cox fit for GSE21501** (parsed from the downloaded SOFT file, HR=0.93, p=0.80, n=102), alongside the two real Phase 2 cohort HRs. Conclusions reflect actual measurements except where noted.

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

*(Note: this correlation is computed at the individual-cell level and does not depend on the patient-level aggressive/reference classification, so it is not affected by the pseudoreplication issue described in Sections 4–5 below.)*

## 4. Lipid Program Cell-of-Origin (Figure 3B)

**METHODOLOGY CORRECTION (2026-07-06):** The original version of this analysis classified patients as aggressive/reference (correctly, at the patient level) but then tested lipid scores by pooling all of a patient's malignant cells together and running a Wilcoxon rank-sum test **treating each cell as an independent observation** (e.g. n=698 vs n=368 "cells" in GSE154778, n=5,730 vs n=3,666 in Peng_et_al). Since these cells come from only 10 (GSE154778) and 17 (Peng_et_al) actual patients, this is pseudoreplication: it inflates the effective sample size by 30–500x and can produce very small p-values from tiny, patient-driven (not necessarily biological) median differences. The analysis has been re-run using **patient-level pseudobulk** (median score per patient, per cell type, requiring ≥5 cells/patient and ≥3 patients per arm) as the unit of statistical testing.

### Patient states per cohort (after re-run)
- GSE154778: 1 aggressive, 2 reference, 7 intermediate (of 10 patients) — **too few patients per arm to test**
- GSE202051: 1 patient total — not testable (as before)
- Peng_et_al: 5 aggressive, 6 reference, 6 intermediate (of 17 patients) — testable

### Summary by Cell Type and Score (Peng_et_al only; GSE154778/GSE202051 underpowered at patient level)

**lipid_synthesis_srebp:**
  - malignant_epithelial: obs=up, p_adj=0.019, repro=**YES**, cell_intrinsic=**YES** (n=5 aggressive vs 6 reference patients)
  - caf_fibroblast: obs=up, p_adj=0.273, repro=NO (not significant)
  - myeloid: obs=down, p_adj=0.855, repro=NO
  - endothelial: obs=up, p_adj=0.855, repro=NO

**desaturation_elongation:**
  - malignant_epithelial: obs=up (correct direction), p_adj=0.584, repro=NO (not significant)
  - caf_fibroblast: obs=up, p_adj=0.273, repro=NO
  - myeloid: obs=up, p_adj=0.855, repro=NO
  - endothelial: obs=up, p_adj=0.698, repro=NO

**fatty_acid_uptake_oxidation:**
  - malignant_epithelial: obs=up (wrong direction, expected down), p_adj=0.410, repro=NO
  - caf_fibroblast: obs=up, p_adj=0.273, repro=NO
  - myeloid: obs=up, p_adj=0.855, repro=NO
  - endothelial: obs=down (correct direction), p_adj=0.698, repro=NO

**Summary:** 1/3 lipid scores reproduced (correct direction + FDR<0.05) in the one adequately powered cohort (Peng_et_al, malignant cells): lipid_synthesis_srebp is elevated in aggressive malignant cells (p_adj=0.019), consistent with a cell-intrinsic effect. GSE154778 previously reported all three lipid scores as significantly reversed in malignant cells (p<0.001), but that result was driven by only 1 aggressive vs 2 reference patients tested as if 698 vs 368 independent cells — it cannot be evaluated with real patient-level power and should not be treated as a finding.

**Revised conclusion:** The earlier claim that "lipid rewiring does not reproduce and reverses direction in malignant cells" was largely a pseudoreplication artifact. With patient-level testing, the only well-powered comparison (Peng_et_al, 5 vs 6 patients) shows lipid_synthesis_srebp elevated in aggressive malignant cells in the expected direction and significant after FDR correction — weak support *for* cell-intrinsic lipid synthesis rewiring, not against it. Desaturation/elongation trends in the expected direction but is not significant; FA-oxidation does not reproduce. Given n=5 vs 6 patients, this should be read as suggestive, not confirmatory — more patients (or matched bulk+single-cell data, as previously recommended) are needed.

## 5. CAF/EMT Cell-of-Origin (Figures 3C, 3D)

**METHODOLOGY CORRECTION (2026-07-06):** Same pseudoreplication issue as Section 4 applied to the EMT malignant-cell test (previously n=698 vs 368 and n=5,730 vs 3,666 "cells" from only 10 and 17 patients respectively). Re-run using patient-level pseudobulk. The CAF subtype proportions (below) have also been changed from pooled-cell proportions to the **mean of each patient's own subtype proportion**, so that patients contributing more CAF cells no longer dominate the estimate.

### EMT Score in Malignant Cells (patient-level pseudobulk)

  - GSE154778: **not testable** (1 aggressive vs 2 reference patients; below the 3-per-arm minimum) — the previously reported p<0.001 result for this cohort was a pseudoreplication artifact and should be disregarded.
  - Peng_et_al (5 aggressive vs 6 reference patients): median diff = -0.135, p_raw=0.0062, p_adj(BH)=0.019, direction = **down** (opposite of expected "up")

**EMT conclusion:** Unlike the lipid finding, the EMT reversal **survives** patient-level correction in the one testable cohort: EMT score is still significantly *lower* in aggressive malignant cells than reference (Peng_et_al, n=5 vs 6 patients, FDR-adjusted p=0.019). This does not support malignant-cell-intrinsic EMT as the driver of the bulk EMT signal — but note this is now based on a single cohort with only 11 patients total, not "2/2 cohorts" as previously claimed; GSE154778 cannot corroborate or refute it with real patient-level power.

### CAF Subtype Analysis (mean of per-patient proportions)

**GSE154778** (aggressive: 1 patient, reference: 2 patients — descriptive only, not a statistical comparison):
  aggressive: myCAF=6.55%, iCAF=2.92%, apCAF=0.54%
  reference: myCAF=13.53%, iCAF=5.84%, apCAF=0.63%

**Peng_et_al** (aggressive: 5 patients, reference: 6 patients):
  aggressive: myCAF=18.68%, iCAF=9.54%, apCAF=1.19%
  reference: myCAF=23.73%, iCAF=9.95%, apCAF=1.61%

*(Note: absolute proportions shifted from the previous pooled-cell version, e.g. GSE154778 aggressive myCAF was reported as 65.48% before — that number was dominated by whichever single aggressive patient had the most CAF cells. The corrected numbers here average across patients rather than across cells, but with 1–2 patients per arm in GSE154778 they remain descriptive, not inferential.)*

## 6. Purity Adjustment of CAF/EMT Signals (Figure 3E)

*(Unaffected by the pseudoreplication fix — this analysis already operated at one row per bulk sample/patient. Re-run 2026-07-06 with a corrected, non-circular purity proxy — see methodology note below.)*

**METHODOLOGY CORRECTION (2026-07-06):** The original 6-gene purity proxy (ACTA2, COL1A1, COL1A2, DCN, FAP, PDGFRB) was a strict subset of the 10-gene CAF signature (ACTA2, COL1A1, COL1A2, COL3A1, DCN, LUM, FAP, PDGFRB, TAGLN, POSTN) — regressing the CAF score on a covariate built from 60% of its own genes, which mechanically inflates apparent attenuation regardless of true confounding. The proxy has been rebuilt (`purity_stromal_immune` in `config/gene_sets.yml`) using 8 genes spanning immune, generic stromal, and vascular compartments (PTPRC, THY1, SPARC, CD2, CD14, CD34, LAPTM5, SERPING1) — verified disjoint from both the CAF and EMT signatures, and confirmed present in all three cohort expression matrices (8/8 coverage each).

| Cohort | Score | n | Unadjusted coef (p) | Purity-adjusted coef (p) | Change |
|---|---|---|---|---|---|
| GSE79668 | caf_score | 49 | 0.419 (p=0.152, ns) | 0.048 (p=0.802, ns) | Large drop, but neither was significant to begin with |
| GSE79668 | emt_score | 49 | 0.395 (p=0.128, ns) | 0.047 (p=0.760, ns) | Large drop, but neither was significant to begin with |
| GSE71729 | caf_score | 145 | 0.527 (p=0.0003) | 0.461 (p<0.0001) | **Trivial (13% drop), remains highly significant** |
| GSE71729 | emt_score | 145 | 0.493 (p<0.0001) | 0.446 (p<0.0001) | **Trivial (9% drop), remains highly significant** |
| GSE62165 | caf_score | 131 | 0.440 (p=0.0075) | 0.431 (p<0.0001) | **Trivial (2% drop), remains highly significant** |
| GSE62165 | emt_score | 131 | 0.560 (p=0.0001) | 0.553 (p<0.0001) | **Trivial (1% drop), remains highly significant** |

**Revised conclusion:** This reverses the previous headline finding. With the circularity removed, the CAF and EMT associations **survive purity adjustment essentially intact and remain highly significant** in the two adequately powered cohorts (GSE71729 n=145, GSE62165 n=131) — coefficients shrink by only 1–13%, and in GSE62165 the adjusted p-value is actually *more* significant than the unadjusted one (purity absorbs residual noise rather than explaining away the group effect). Only in the smallest cohort (GSE79668, n=49) does the coefficient drop substantially — but that cohort's CAF/EMT association was not statistically significant even before adjustment (p=0.13–0.15), so this is better read as "an already-underpowered, non-significant effect got smaller," not as confirmation of stromal confounding. **The previous claim that "the CAF bulk signal is largely explained by tumor purity variation" was itself substantially an artifact of the circular proxy — the corrected analysis instead supports a real, largely purity-independent CAF/EMT association in bulk PDAC tumors.**

## 7. Protein-Level Lipid Validation — CPTAC-PDA (Figure 3F)

**Data status:** REAL DATA — umich proteomics (145 tumor samples) and BCM proteomics (105 tumor samples with matched transcriptomics), WashU transcriptomics for group assignment (9/9 hypoxia genes, 10/11 acinar genes). Group assignment: 46 aggressive / 46 reference (umich); 41 aggressive / 23 reference (BCM — unequal due to smaller overlap with 140-sample transcriptomics set). Missing proteins in both sources: ELOVL6, CPT1B (not detected).

### umich results (n=46 aggressive, n=46 reference)

| Gene set | Protein | Expected | umich direction | umich FDR | Concordant |
|---|---|---|---|---|---|
| lipid_synthesis_srebp | FASN | up | up | 0.015 | ✓ |
| lipid_synthesis_srebp | ACACA | up | up | 0.001 | ✓ |
| lipid_synthesis_srebp | SQLE | up | up | 0.050 | ✓ |
| lipid_synthesis_srebp | HMGCR | up | up | 0.516 | — (ns) |
| lipid_synthesis_srebp | SREBF1 | up | down | 0.516 | ✗ |
| lipid_synthesis_srebp | ACLY | up | down | 0.169 | ✗ |
| desaturation_elongation | SCD | up | up | 0.862 | — (ns) |
| desaturation_elongation | FADS1 | up | up | 0.870 | — (ns) |
| desaturation_elongation | FADS2 | up | down | 0.617 | ✗ |
| fatty_acid_uptake_oxidation | ACADL | down | down | <0.001 | ✓ |
| fatty_acid_uptake_oxidation | CD36 | down | down | 0.705 | — (ns) |
| fatty_acid_uptake_oxidation | FABP5 | down | down | 0.078 | — (borderline) |
| fatty_acid_uptake_oxidation | FABP4 | down | up | 0.870 | ✗ |
| fatty_acid_uptake_oxidation | CPT1A | down | up | 0.705 | ✗ |
| fatty_acid_uptake_oxidation | HADHA | down | up | 0.870 | ✗ |

**umich overall:** 9/15 directionally concordant; 4/15 FDR<0.05 and concordant (FASN, ACACA, SQLE, ACADL).

### BCM replication (n=41 aggressive, n=23 reference)

| Protein | Expected | BCM direction | BCM FDR | Replicates umich |
|---|---|---|---|---|
| ACACA | up | up | ns | ✓ (direction) |
| FASN | up | up | ns | ✓ (direction) |
| SQLE | up | up | ns | ✓ (direction) |
| HMGCR | up | up | ns | ✓ (direction) |
| SREBF1 | up | down | ns | — (both discordant) |
| ACLY | up | up | ns | — (direction flipped vs umich) |
| SCD | up | up | ns | ✓ (direction) |
| FADS1 | up | up | ns | ✓ (direction) |
| FADS2 | up | up | ns | — (direction flipped vs umich) |
| ACADL | down | down | **FDR<0.05** | **✓ (replicated, FDR<0.05 both)** |
| CD36 | down | up | ns | ✗ |
| FABP5 | down | up | ns | ✗ |
| FABP4 | down | up | ns | — (both discordant in same direction) |
| CPT1A | down | down | ns | — (direction flipped vs umich) |
| HADHA | down | up | ns | — (both discordant) |

**BCM overall:** 10/15 directionally concordant with expected direction; 1/15 FDR<0.05+concordant.

**Cross-source replication summary:** 7/15 proteins concordant in both sources; **1 protein fully replicated with FDR<0.05 in both: ACADL** (fatty acid β-oxidation, consistently reduced in aggressive PDAC). Among the 4 umich-significant proteins, 3 (FASN, ACACA, SQLE) show directional concordance in BCM but do not reach FDR<0.05, likely reflecting the smaller and unbalanced BCM reference group (n=23 vs n=46 in umich). SREBF1 and HADHA are discordant in both sources.

**Interpretation:** The core lipogenic enzymes and the β-oxidation enzyme ACADL show the most consistent protein-level evidence. ACACA (FDR=0.001 umich) and FASN (FDR=0.015 umich) are elevated, and ACADL (FDR<0.001 in both umich and BCM) is reduced in aggressive tumors — the two most mechanistically anchored predictions from Phase 2. SQLE (cholesterol synthesis) is elevated in umich (FDR=0.050) with directional concordance in BCM. ACADL is the single protein achieving independent replication across both CPTAC proteomics centers. SREBF1 itself is not elevated at protein level in either source, consistent with post-translational rather than abundance-level regulation of SREBP1 activity. Overall, the protein data partially corroborates the lipid synthesis/FA-oxidation axis from bulk RNA-seq, with the strongest evidence at enzymatic nodes most directly linked to lipid flux.

## 8. Expanded Survival Meta-Analysis (Figure 3G)

*(Unaffected by the pseudoreplication/purity fixes. Re-run 2026-07-06: GSE21501 now uses a real Cox fit instead of a hardcoded literature value — see methodology note below.)*

**METHODOLOGY CORRECTION (2026-07-06):** Previously, `run_expanded_survival_meta.py` found the real `GSE21501_family.soft.gz` file on disk but never parsed it — it returned a hardcoded literature value (HR=1.40, p=0.21) labeled `is_simulated: False`, which was misleading. A new script, `scripts/parse_gse21501_survival.py`, streams the ~6M-line, 188MB family SOFT file once (132 samples on the GPL4133 Agilent platform) and extracts: (1) real clinical fields from `!Sample_characteristics_ch2` — `os_time`, `os_event` (102/132 samples have both; the remaining 30 lack clinical annotation per the series description); (2) expression values for the hypoxia (9/9 found) and acinar_identity (9/11 found; AMY2A and CELA3A are absent from this platform) marker genes, via the platform's probe→GENE_SYMBOL table. `run_expanded_survival_meta.py` then reproduces the exact Phase 2 methodology (z-score each marker gene across the cohort, average into hypoxia/acinar scores, median-split, define "aggressive" = hypoxia-high & acinar-low) and fits a univariable Cox model (lifelines `CoxPHFitter`, aggressive vs. everyone else) — the same approach used for GSE79668/GSE71729 in Phase 2.

### Individual Cohorts
- **GSE79668**: HR=1.232, 95% CI [0.629-2.413], p=0.543 (Phase 2, unchanged)
- **GSE71729**: HR=1.079, 95% CI [0.661-1.761], p=0.761 (Phase 2, unchanged)
- **GSE21501**: HR=0.927, 95% CI [0.519-1.655], p=0.798 (**real Cox fit**, n=102, 66 events, 24 aggressive — replaces the previous hardcoded HR=1.40, p=0.21)

### Pooled Estimates
- **pooled_fixed_effects**: HR=1.061, 95% CI [0.765-1.471], p=0.723
- **pooled_random_effects**: HR=1.061, 95% CI [0.765-1.471], p=0.723 (I²=0%)

**Revised conclusion:** The real GSE21501 result does **not** match the previously assumed "HR~1.3–1.4, consistent with the Phase 2 trend" — it comes out essentially null and, if anything, in the *opposite* direction (HR=0.93, i.e. very slightly protective, though nowhere near significant, p=0.80). This is no longer "3/3 cohorts with HR>1" as previously implied; it is 2/3, with the third real cohort centered almost exactly on the null. The pooled estimate consequently drops from the previous HR=1.22 [0.89–1.68] to **HR=1.06 [0.77–1.47], p=0.72** — moving further from, not closer to, significance. **The survival trend should now be described as inconsistent across the three available cohorts, not as a uniform directional trend that merely lacks power.** As a sanity check, the 24/102 (23.5%) fraction of GSE21501 samples falling in the hypoxia-high/acinar-low quadrant closely matches the ~25% expected under independence — consistent with the Phase 3 single-cell finding (Section 3) that hypoxia and acinar status behave as largely independent axes, which by construction puts any median-split "aggressive" group at roughly one-quarter of a cohort regardless of platform.

## 9. Immune / Endothelial QC Check

**METHODOLOGY CORRECTION (2026-07-06):** This section draws from the same `figure3B_lipid_cell_of_origin_statistics.tsv` table used in Section 4, which has already been re-computed at the patient-pseudobulk level (≥5 cells/patient, ≥3 patients/arm) — no separate script exists for this section. GSE154778 could not be tested for any cell type (only 1 aggressive patient total), so all figures below are from Peng_et_al (5 aggressive vs 6 reference patients) only.

- **myeloid** (Peng_et_al, n=5 vs 6 patients): 0/3 lipid comparisons significant at FDR<0.05 (lipid_synthesis_srebp p_adj=0.855; desaturation_elongation p_adj=0.855; fatty_acid_uptake_oxidation p_adj=0.855, wrong direction). Previously reported as "2/3 significant" from a pseudoreplicated cell-level test (n in the thousands); at the patient level there is no significant myeloid lipid signal. Low off-target signal consistent with tumor-driven lipid rewiring, though this could equally reflect the small patient count (n=11 total) rather than a true null.
- **endothelial** (Peng_et_al, n=5 vs 6 patients): 0/3 lipid comparisons significant at FDR<0.05 (lipid_synthesis_srebp p_adj=0.855; desaturation_elongation p_adj=0.698; fatty_acid_uptake_oxidation p_adj=0.698, correct direction but ns). Previously reported as "2/3 significant" from the same flawed cell-level test. No significant endothelial off-target lipid signal at the patient level.
- **caf_fibroblast** (Peng_et_al, n=5 vs 6 patients; see also Section 4): 0/3 lipid comparisons significant at FDR<0.05 (p_adj = 0.273 for all three). Previously reported as significant under the cell-level test. No significant CAF lipid signal at the patient level.

**Revised conclusion:** After correcting for pseudoreplication, none of the three non-malignant cell types (myeloid, endothelial, CAF) show a statistically significant lipid difference between aggressive and reference patients in the one testable cohort. This removes the previous "myeloid/endothelial/CAF lipid signal" as a confound to worry about, but the absence of significance here is at least partly a power issue (n=11 patients) rather than confirmed evidence of no effect — it should not be read as a clean negative any more than the malignant-cell result should be read as a clean positive.

## 10. Final Mechanistic Classification

**Classification: Partially resolved / mixed (revised 2026-07-06 after pseudoreplication fix, purity-proxy circularity fix, and real GSE21501 survival fix)**

Based on **real single-cell RNA-seq data** (GSE154778: 8,000 primary tumor cells, 10 patients; GSE202051: 2,607 snRNA-seq nuclei, 1 patient; Peng et al.: 43,888 cells, 17 patients), real Phase 2 bulk cohorts for purity adjustment, and **real CPTAC proteomics** (umich 145 + BCM 105 tumor samples):

1. **Hypoxia/acinar co-occurrence (REAL DATA — negative result, unaffected by the fix):** Pearson r(hypoxia, acinar) in malignant cells = +0.06, −0.03, −0.02 across three real scRNA-seq cohorts (all near-zero; p > 0.05 in 2/3 cohorts). The fraction of malignant cells in the hypoxia-high/acinar-low quadrant is 24–27%, indistinguishable from the 25% expected under independence. This still supports the composite-artifact interpretation: hypoxia-high and acinar-low programs largely occur in separate malignant cell subpopulations, not the same cells.

2. **Lipid rewiring cell-of-origin (REVISED — weak positive signal, not a clean negative):** The original claim that lipid scores reverse direction in malignant cells was driven almost entirely by pseudoreplication in GSE154778 (1 vs 2 patients tested as 698 vs 368 "independent" cells) — that result is now excluded as untestable. In the one adequately powered cohort (Peng_et_al, 5 vs 6 patients), lipid_synthesis_srebp is significantly elevated in aggressive malignant cells in the **expected** direction (FDR=0.019), consistent with (not against) cell-intrinsic lipid synthesis rewiring. Desaturation/elongation trends in the expected direction (not significant); FA-oxidation does not reproduce. This remains underpowered (n=11 patients total) and does not resolve cell-of-origin definitively, but the corrected direction of evidence is the opposite of what was previously reported.

3. **CAF/EMT cell-of-origin (REVISED — EMT reversal partially holds up, CAF descriptive only):** After patient-level correction, EMT score is still significantly LOWER in H-hi/A-lo malignant cells in the one testable cohort (Peng_et_al: 5 vs 6 patients, FDR=0.019) — opposite the expected direction, same qualitative conclusion as before but based on 1/1 testable cohorts rather than the previously claimed 2/2 (GSE154778 could not be tested with only 1 aggressive patient). CAF subtype proportions (recomputed as mean-of-per-patient proportions rather than pooled-cell proportions) remain broadly similar between aggressive and reference groups in Peng_et_al (5 vs 6 patients); GSE154778 CAF proportions are now explicitly descriptive only (1 vs 2 patients).

4. **Purity adjustment of bulk CAF/EMT (REVISED 2026-07-06 — circularity fixed; CAF/EMT signal is largely real, not a purity artifact):** The original purity proxy was rebuilt to be disjoint from the CAF/EMT signatures (see Section 6). With the corrected proxy, CAF and EMT effects survive purity adjustment essentially intact and remain highly significant in the two adequately powered cohorts (GSE71729 n=145: caf coef 0.53→0.46, p<0.0001; GSE62165 n=131: caf coef 0.44→0.43, p<0.0001; EMT similarly stable in both). Only in the smallest, already-underpowered cohort (GSE79668, n=49, unadjusted p=0.13–0.15) does the coefficient drop substantially. This **reverses** the earlier conclusion that the CAF signal was "largely explained by tumor purity variation" — that conclusion was itself substantially an artifact of the original proxy's circularity (60% gene overlap with the CAF signature). The corrected analysis instead supports a real, largely purity-independent CAF/EMT association in bulk PDAC.

5. **Protein-level validation (REAL DATA — partially concordant, independently replicated, unaffected by the fix):** ACADL replicated FDR<0.05 in both umich and BCM; FASN, ACACA, SQLE elevated and significant in umich (FDR<0.05), directionally concordant but underpowered in BCM.

6. **Survival (REVISED 2026-07-06 — real GSE21501 fit changes the trend from consistent to inconsistent):** GSE21501's entry was previously a hardcoded literature value (HR=1.40, p=0.21) despite the real SOFT file being on disk and unparsed. A real Cox fit (n=102, 66 events, same hypoxia-high/acinar-low definition as Phase 2) gives **HR=0.93, p=0.80** — essentially null, and in the opposite direction from the other two cohorts. The pooled estimate drops from HR=1.22 [0.89–1.68] to **HR=1.06 [0.77–1.47], p=0.72**. This is no longer "HR>1 in all cohorts, just underpowered" — it is 2/3 cohorts with HR>1 and one centered on the null, which is a materially weaker survival story than previously reported.

**Recommended wording (revised 2026-07-06 after purity-proxy fix and real GSE21501 survival fix):**

> Single-cell analysis indicates the bulk hypoxia-high/acinar-low signature is a composite of two largely independent malignant-cell axes (r ≈ 0 at the cell level in three cohorts) rather than a single coherent program. Patient-level (pseudobulk) re-analysis, correcting an earlier cell-level pseudoreplication error, finds a modest cell-intrinsic elevation of lipid synthesis genes in aggressive malignant cells in the one adequately powered cohort (n=5 vs 6 patients, FDR=0.019) — weak evidence in favor of, not against, tumor-intrinsic lipid rewiring — while EMT score is lower, not higher, in aggressive malignant cells in that same cohort (FDR=0.019). Contrary to the earlier report, the bulk CAF and EMT associations **survive adjustment for tumor purity** using a corrected, non-circular purity proxy, remaining highly significant in the two adequately powered bulk cohorts (GSE71729, GSE62165) with only 1–13% coefficient attenuation — the CAF/EMT signal in bulk PDAC is not primarily a tumor-purity artifact. A real Cox fit on the previously unparsed GSE21501 cohort gives HR=0.93 (p=0.80), pulling the pooled survival estimate to HR=1.06 [0.77–1.47], p=0.72 — the survival association across the three available cohorts should now be described as inconsistent in direction, not merely underpowered. All single-cell cell-of-origin conclusions are based on very small patient numbers (≤17 per cohort, ≤11 per comparison) and require replication in larger and, ideally, matched bulk+single-cell cohorts.

## 11. Limitations

1. **No matched bulk + single-cell patient data.** The three single-cell cohorts (GSE154778, GSE202051, Peng et al.) are independent from the bulk Phase 2 cohorts and share no patients. The within-single-cell-dataset H-hi/A-lo grouping is not equivalent to the between-dataset Phase 2 aggressive group definition. Lipid cell-of-origin conclusions require matched data.

2. **GSE202051 has only 1 unique patient ID.** The `pid` column in GSE202051 identifies a single patient, making patient-level aggressive/reference classification inapplicable. Cell-level H-hi/A-lo scoring is used instead.

3. **Dissociation stress artifacts.** GSE154778 (stress ratio 7.95) and Peng_et_al (stress ratio 9.84) show elevated dissociation-stress genes (FOS, JUN, HSP family). Hypoxia scores in single cells may be partially inflated by dissociation stress rather than true in vivo hypoxia.

4. **BCM proteomics replication uses an unbalanced reference group (n=23 vs n=46 in umich).** Of the 140-sample transcriptomics-defined groups (46 aggressive / 46 reference / 48 other), only 105 BCM samples overlap, yielding 41 aggressive and only 23 reference. This imbalance reduces power for BCM FDR thresholds. ELOVL6 and CPT1B are absent from both umich and BCM datasets.

5. **(Resolved 2026-07-06) GSE21501 survival now uses a real Cox fit** (HR=0.93, p=0.80, n=102, 66 events) instead of a hardcoded literature value. Residual caveats: only 102/132 samples have clinical annotation (30 lack it per the series description), and 2/11 acinar genes (AMY2A, CELA3A) are absent from this platform, so the acinar score rests on 9 genes rather than 11 for this cohort specifically.

6. **(Resolved 2026-07-06) Purity proxy circularity.** The proxy now uses 8 genes (PTPRC, THY1, SPARC, CD2, CD14, CD34, LAPTM5, SERPING1) verified disjoint from the CAF/EMT signatures, replacing the original 6-gene proxy that was a strict subset of the CAF signature. It remains an expression-mean heuristic rather than a validated ESTIMATE/CIBERSORTx purity estimate — that upgrade is still open (see Next Steps item 4) — but the collinearity that made the CAF-attenuation conclusion partly circular is fixed.

7. **GSE202051 is snRNA-seq (nuclei), not scRNA-seq.** Nuclear transcriptomes have lower cytoplasmic gene expression (especially mitochondrial). CAF and stromal gene scores may be underestimated.

8. **Small malignant cell counts in GSE202051 (n=2,039 total; 1 patient).** Statistical comparisons for this cohort are underpowered for patient-level inferences.

9. **(New) Cell-of-origin sample sizes are small even after the pseudoreplication fix.** Only Peng_et_al (17 patients) has enough aggressive/reference patients to test at all; GSE154778 (10 patients) and GSE202051 (1 patient) cannot support patient-level significance tests for lipid/EMT cell-of-origin. The corrected Section 4/5 conclusions rest on n=5 vs 6 patients in a single cohort and should be treated as suggestive pending more patients or matched bulk+single-cell data.

10. **(Resolved) Section 9 (Immune/Endothelial QC)** has been re-run at the patient level (see Section 9). With only Peng_et_al testable (5 vs 6 patients), the "no significant off-target lipid signal" conclusion is now correctly powered but rests on a single cohort — obtaining more patients per cohort would strengthen it.

## 12. Recommended Manuscript Wording

The recommended wording in Section 10 reflects what the real data support after the patient-level correction. Do not use the previous wording claiming malignant lipid rewiring "does not reproduce" or "reverses direction" — the corrected, adequately-powered comparison points the other way (weak positive support) for lipid synthesis specifically.

For the survival finding (revised 2026-07-06 — do not use the previous "consistent directional trend" wording, which relied on a hardcoded GSE21501 value):
> "Across the three cohorts with available survival data, the hypoxia-high/acinar-low state showed an inconsistent association with overall survival (HR=1.23 and 1.08 in two cohorts, HR=0.93 in a third), and the pooled estimate did not reach statistical significance (pooled HR = 1.06, 95% CI [0.77–1.47], p = 0.72). This is better characterized as no clear survival association in the available data than as a consistent trend obscured by low power."

## 13. Next Steps

1. **Obtain matched bulk + single-cell data** for the same patients (e.g., from PDAC cohorts with paired RNA-seq and scRNA-seq) to properly test lipid cell-of-origin with matched group assignments.

2. ~~**Replicate CPTAC protein results in the BCM proteomics source**~~ **COMPLETE.** BCM replication done: ACADL replicated FDR<0.05 in both umich and BCM; FASN/ACACA/SQLE directionally concordant but underpowered in BCM (n=23 reference). See Section 7.

3. ~~**Real GSE21501 survival analysis**~~ **COMPLETE.** `scripts/parse_gse21501_survival.py` now streams the real SOFT file and `run_expanded_survival_meta.py` fits a real Cox model (HR=0.93, p=0.80, n=102). See Section 8. Remaining open item: only 102/132 GSE21501 samples have clinical annotation and only 9/11 acinar genes are on this platform — a larger or better-annotated cohort would strengthen this estimate.

4. **Implement proper purity deconvolution:** Apply ESTIMATE (R package: `library(estimate); filterCommonGenes(); estimateScore()`) or TIMER2/CIBERSORTx for validated tumor purity estimates. The gene-overlap circularity is already fixed (Section 6), but a validated deconvolution method would still be more rigorous than the current 8-gene expression-mean heuristic.

5. **Spatial transcriptomics:** 10x Visium PDAC datasets (e.g., from the 10x Genomics public datasets) could directly test whether hypoxic and acinar-low regions spatially co-localize within the same tumor section, providing an orthogonal test of the composite-state vs. co-occurring-state question.

6. **Dissociation stress correction:** Apply computational correction for dissociation artifacts (e.g., van den Brink et al. 2017 gene list) before re-scoring hypoxia in GSE154778 and Peng et al.

7. ~~**Re-run Section 9 (Immune/Endothelial QC) at the patient level**~~ **COMPLETE.** See Section 9. Obtaining more single-cell cohorts (or more patients within existing cohorts) would still let GSE154778 and GSE202051 support patient-level lipid/EMT cell-of-origin tests, which they currently cannot (1–2 patients per arm).

---
*Report generated by Phase 3 pipeline. Single-cell data: REAL (GSE154778, GSE202051, Peng et al.). CPTAC proteomics: REAL (umich 145 + BCM 105 tumor samples; ACADL independently replicated FDR<0.05 in both). Survival cohort GSE21501: REAL Cox fit (HR=0.93, p=0.80, n=102) computed from the parsed SOFT file, replacing the previous hardcoded literature value (see Section 8/11). Purity adjustment: REAL Phase 2 bulk expression data with an 8-gene purity proxy verified disjoint from the CAF/EMT signatures (see Section 6/11). Sections 4, 5, and 9 re-analyzed 2026-07-06 using patient-level pseudobulk to correct a cell-level pseudoreplication error; Section 6 re-analyzed 2026-07-06 with the corrected purity proxy, reversing the earlier "CAF signal is purity-driven" conclusion; Section 8 re-analyzed 2026-07-06 with the real GSE21501 Cox fit, reversing the earlier "consistent survival trend" framing.*
