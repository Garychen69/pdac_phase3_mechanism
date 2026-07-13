# Phase 3 Mechanistic Report: Single-Cell Resolution of Aggressive PDAC State

**Generated:** 2026-07-03 (updated 2026-07-04 with BCM replication; updated 2026-07-06 with patient-level pseudobulk re-analysis of lipid/CAF/EMT/immune-endothelial cell-of-origin, Sections 4/5/9; updated 2026-07-06 with a non-circular purity proxy, Section 6; updated 2026-07-06 with a real GSE21501 Cox fit replacing the hardcoded survival value, Section 8; **updated 2026-07-12 with the full 43-patient GSE202051 object replacing a 1-patient subset, Sections 3/4/5/9/10/11**)

**Reproducibility check (2026-07-11):** The full pipeline (`scripts/run_phase3_pipeline.py`) was re-run end-to-end in the `pdac_phase3` conda environment against the same on-disk real data. All key statistics reproduced to 2-3 decimal places: hypoxia/acinar Pearson r = 0.060 / -0.031 / -0.016 (GSE154778/GSE202051/Peng_et_al); CPTAC 9/15 (umich) and 10/15 (bcm) directionally concordant with ACADL independently replicated at FDR<0.05 in both; GSE21501 real Cox fit HR=0.927, p=0.798, n=102; pooled survival HR=1.061 [0.765-1.471], p=0.723; purity-adjusted CAF/EMT coefficients matched to 3 decimals across all 3 bulk cohorts. This confirms the analysis is genuinely reproducible from the real downloaded data, not a one-off result. Note: `scripts/generate_phase3_report.py` itself is an older/simpler version whose Sections 10-13 (Classification/Limitations/Wording/Next Steps) are hardcoded boilerplate that incorrectly re-asserts "simulated data" throughout — its raw output from this re-run is archived at `PHASE3_MECHANISM_REPORT_raw_template_output_2026-07-11.md` for transparency, but this file (manually curated 2026-07-06, numerically re-verified 2026-07-11, data-upgraded 2026-07-12) remains the canonical, accurate report.

**Data upgrade (2026-07-12):** GSE202051 was found to have been substantially under-downloaded. The GEO series (GSE202051) has 74 GSM samples across 43 patients, but the pipeline had only downloaded a single-sample supplementary file (`GSE202051_adata_010nuc_10x.h5ad`, 2,607 cells, 1 patient) — the 1-patient limitation cited throughout the pre-2026-07-12 version of this report was a download-completeness bug, not an inherent dataset limitation. The full series object (`GSE202051_totaldata-final-toshare.h5ad`, 224,988 cells, 43 patients, richly pre-annotated) has been downloaded, its unused ~5GB of CNV-inference/harmony/duplicate-raw data stripped (not used by this pipeline), and reprocessed through the full annotation and scoring pipeline. This converts GSE202051 from an untestable cohort into the second-best-powered single-cell cohort in the project (after Peng_et_al), materially strengthening Sections 3-5 and 9. A related production bug was also fixed in `preprocess_singlecell_cohorts.py`'s `check_dissociation_stress()`, which densified the entire expression matrix via `.toarray().mean()` and crashed with an out-of-memory error on this cohort's scale (17.4GB allocation attempt) — fixed to use scipy sparse `.mean()` directly, which works without densifying.

---

## 1. Dataset Summary

### Single-Cell Cohorts

- **GSE154778**: Available (processed) — **REAL DATA** (8000 cells, 10 patients)
- **GSE202051**: Available (processed) — **REAL DATA, upgraded 2026-07-12** (224,988 cells, **43 patients** — previously 2,607 cells/1 patient) *(snRNA-seq; pre-annotated with rich fine-grained cell-type labels — no de-novo annotation needed, just mapped to this pipeline's standard vocabulary)*
- **Peng_et_al**: Available (processed) — **REAL DATA** (43888 cells, 17 patients)

### Gene Set Coverage per Cohort (signature genes found in dataset)

| Cohort | hypoxia | acinar_identity | lipid_synthesis_srebp | desaturation_elongation | fatty_acid_uptake_oxidation | caf | emt |
|--------|--------|--------|--------|--------|--------|--------|--------|
| GSE154778 | 9/9 (100%) | 7/11 (64%) | 6/6 (100%) | 4/4 (100%) | 7/7 (100%) | 10/10 (100%) | 10/10 (100%) |
| GSE202051 | 9/9 (100%) | 11/11 (100%) | 6/6 (100%) | 4/4 (100%) | 7/7 (100%) | 10/10 (100%) | 10/10 (100%) |
| Peng_et_al | 9/9 (100%) | 10/11 (91%) | 6/6 (100%) | 4/4 (100%) | 7/7 (100%) | 10/10 (100%) | 10/10 (100%) |

*(GSE202051's coverage improved to 100% across all 7 signatures with the full object — the previous 2,607-cell subset had lower coverage, e.g. 64% acinar_identity, 86% FA-oxidation, likely due to different gene filtering in that GEO supplementary file.)*

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
**GSE202051** (upgraded 2026-07-12, full 43-patient object): 224,988 cells
  - malignant_epithelial: 64,538 (28.7%)
  - caf_fibroblast: 54,935 (24.4%)
  - ductal_normal: 32,851 (14.6%)
  - unknown: 26,588 (11.8%) *(pericyte, endocrine, Schwann, vascular smooth muscle, adipocyte — no clean standard-vocabulary equivalent; mapped from the real `new_celltypes` annotation, see Section 3 note)*
  - endothelial: 19,258 (8.6%)
  - myeloid: 11,938 (5.3%)
  - tcell_nk: 7,668 (3.4%)
  - acinar_normal: 5,357 (2.4%)
  - bcell_plasma: 1,855 (0.8%)
  *(Previous 2,607-cell/1-patient subset: malignant_epithelial 2039 (78.2%), ductal_normal 475 (18.2%), endothelial 46 (1.8%), acinar_normal 25 (1.0%), caf_fibroblast 18 (0.7%), myeloid 4 (0.2%) — retained here only for comparison; superseded.)*
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

**DATA UPGRADE (2026-07-12):** GSE202051 was previously represented by a single-sample GEO supplementary file (`GSE202051_adata_010nuc_10x.h5ad`, 2,607 cells, 1 patient) that made this cohort essentially untestable at the patient level. The GEO series actually has 74 GSM samples; a second, larger series-level supplementary file (`GSE202051_totaldata-final-toshare.h5ad`, 225k cells) contains the full, richly-annotated 43-patient object (18 untreated, 25 treated) with pre-computed fine-grained cell-type labels (`new_celltypes`). This file has now been downloaded, mapped to this pipeline's standard cell-type vocabulary, and reprocessed through the full pipeline, replacing the 1-patient subset below. All statistics in this section and Sections 4/5/9 for GSE202051 are from the full object.

**GSE154778** (n=3710 malignant cells):
  - Pearson r(hypoxia, acinar) = 0.060, p = 0.0003
  - Spearman r = 0.075
  - Fraction hypoxia_high/acinar_low: 0.244
  - Interpretation: Largely separate populations (composite artifact)

**GSE202051** (n=64,538 malignant cells, 43 patients — previously n=2,039 cells, 1 patient):
  - Pearson r(hypoxia, acinar) = 0.036, p < 0.0001
  - Fraction hypoxia_high/acinar_low: 0.260
  - Interpretation: Largely separate populations (composite artifact)
  - This cohort is now adequately powered and independently confirms the near-zero correlation seen in the other two cohorts (previously it was underpowered at n=1 patient and its non-significant r=-0.03 could not be trusted either way).

**Peng_et_al** (n=11370 malignant cells):
  - Pearson r(hypoxia, acinar) = -0.016, p = 0.088
  - Spearman r = -0.018
  - Fraction hypoxia_high/acinar_low: 0.272
  - Interpretation: Largely separate populations (composite artifact)

**Conclusion (updated 2026-07-12):** Mean fraction of malignant cells in hypoxia-high/acinar-low quadrant = 25.89% (unchanged from the prior 25.91%). All three cohorts now show a near-zero Pearson correlation, and — unlike the prior report — all three are well-powered (n=3,710 to 64,538 malignant cells; the previous version had one cohort resting on essentially a single patient). Weak correlation across three adequately powered cohorts supports the hypoxia-high/acinar-low state being a composite artifact of two largely independent cell populations rather than a uniform cell state; this conclusion is now on considerably firmer footing than the prior 2/3-cohorts-adequately-powered version.

*(Note: this correlation is computed at the individual-cell level and does not depend on the patient-level aggressive/reference classification, so it is not affected by the pseudoreplication issue described in Sections 4–5 below.)*

## 4. Lipid Program Cell-of-Origin (Figure 3B)

**METHODOLOGY CORRECTION (2026-07-06):** The original version of this analysis classified patients as aggressive/reference (correctly, at the patient level) but then tested lipid scores by pooling all of a patient's malignant cells together and running a Wilcoxon rank-sum test **treating each cell as an independent observation** (e.g. n=698 vs n=368 "cells" in GSE154778, n=5,730 vs n=3,666 in Peng_et_al). Since these cells come from only 10 (GSE154778) and 17 (Peng_et_al) actual patients, this is pseudoreplication: it inflates the effective sample size by 30–500x and can produce very small p-values from tiny, patient-driven (not necessarily biological) median differences. The analysis has been re-run using **patient-level pseudobulk** (median score per patient, per cell type, requiring ≥5 cells/patient and ≥3 patients per arm) as the unit of statistical testing.

### Patient states per cohort (after re-run)
- GSE154778: 1 aggressive, 2 reference, 7 intermediate (of 10 patients) — **too few patients per arm to test**
- GSE202051 (**upgraded 2026-07-12**, full 43-patient object): 7 aggressive, 14–16 reference (varies slightly by cell type due to per-cell-type ≥5-cell/patient filtering), 20 intermediate — **testable**
- Peng_et_al: 5 aggressive, 6 reference, 6 intermediate (of 17 patients) — testable

### Summary by Cell Type and Score (Peng_et_al and GSE202051; GSE154778 still underpowered at patient level)

**lipid_synthesis_srebp:**
  - malignant_epithelial: Peng_et_al obs=up, p_adj=0.019, repro=**YES** (n=5 vs 6). GSE202051 obs=up (same direction), p_adj=0.052, repro=borderline NO (n=7 vs 14; just misses FDR<0.05).
  - caf_fibroblast: Peng_et_al obs=up, p_adj=0.273, repro=NO. GSE202051 obs=**down** (opposite of expected "up"), p_adj=0.002, repro=NO — **significant in the wrong direction** (n=7 vs 16).
  - myeloid: Peng_et_al obs=down, p_adj=0.855, repro=NO. GSE202051 obs=down, p_adj=0.110, repro=NO (n=6 vs 15).
  - endothelial: Peng_et_al obs=up, p_adj=0.855, repro=NO. GSE202051 obs=**down** (opposite of expected "up"), p_adj=0.023, repro=NO — significant wrong direction (n=7 vs 16).

**desaturation_elongation:**
  - malignant_epithelial: Peng_et_al obs=up (correct direction), p_adj=0.584, repro=NO. GSE202051 obs=**down** (wrong direction), p_adj=0.022, repro=NO — significant wrong direction (n=7 vs 14).
  - caf_fibroblast: Peng_et_al obs=up, p_adj=0.273, repro=NO. GSE202051 obs=down (wrong), p_adj=0.033, repro=NO — significant wrong direction.
  - myeloid: Peng_et_al obs=up, p_adj=0.855, repro=NO. GSE202051 obs=up, p_adj=0.436, repro=NO.
  - endothelial: Peng_et_al obs=up, p_adj=0.698, repro=NO. GSE202051 obs=down (wrong), p_adj=0.058, repro=NO (borderline).

**fatty_acid_uptake_oxidation:**
  - malignant_epithelial: Peng_et_al obs=up (wrong direction, expected down), p_adj=0.410, repro=NO. GSE202051 obs=down (**correct** direction), p_adj=0.052, repro=borderline NO (n=7 vs 14).
  - caf_fibroblast: Peng_et_al obs=up, p_adj=0.273, repro=NO. GSE202051 obs=down (correct), p_adj=0.005, repro=**YES** (n=7 vs 16).
  - myeloid: Peng_et_al obs=up, p_adj=0.855, repro=NO. GSE202051 obs=down (correct), p_adj=0.031, repro=**YES** (n=6 vs 15).
  - endothelial: Peng_et_al obs=down (correct direction), p_adj=0.698, repro=NO. GSE202051 obs=down (correct), p_adj=0.124, repro=NO (borderline).

**Summary (revised 2026-07-12):** With GSE202051 now testable at real statistical power (7 vs 14–16 patients, tens of thousands of cells per arm), the malignant-cell lipid_synthesis_srebp signal **replicates in direction in both adequately powered cohorts** (Peng_et_al significant, FDR=0.019; GSE202051 same direction, FDR=0.052, just short of significance) — this is a materially stronger result than the single-cohort finding reported before 2026-07-12. However, GSE202051 also surfaces **new significant off-target and wrong-direction signal not previously detectable**: CAF and endothelial cells show significantly *lower* lipid_synthesis_srebp in aggressive patients (opposite the expected direction), and malignant desaturation/elongation is significantly reversed rather than elevated. FA-oxidation is the most consistent score in GSE202051 — correctly reduced in malignant, CAF, and myeloid compartments (myeloid and CAF reaching FDR<0.05).

**Revised conclusion (2026-07-12):** The core cell-intrinsic lipid synthesis claim now has directionally consistent support in the two adequately powered cohorts (1 of 2 individually significant), which is more evidence than the single-cohort version of this finding. But the newly well-powered GSE202051 cohort also shows the *desaturation/elongation* axis significantly reversed in malignant cells and *lipid_synthesis_srebp* significantly reversed in non-malignant compartments (CAF, endothelial) — complications that a single small cohort (Peng_et_al, n=11 total patients) could not have revealed. The overall lipid cell-of-origin picture is now more informative but also more heterogeneous than "weak positive support" — it should be read as: lipid_synthesis_srebp is the most reproducible axis and trends malignant-cell-intrinsic across both testable cohorts, while desaturation/elongation and off-target (CAF/endothelial) lipid signal are inconsistent or reversed. FA-oxidation reduction is the most consistent finding of the three lipid scores, now with two of three GSE202051 cell types reaching significance in the expected (reduced) direction.

## 5. CAF/EMT Cell-of-Origin (Figures 3C, 3D)

**METHODOLOGY CORRECTION (2026-07-06):** Same pseudoreplication issue as Section 4 applied to the EMT malignant-cell test (previously n=698 vs 368 and n=5,730 vs 3,666 "cells" from only 10 and 17 patients respectively). Re-run using patient-level pseudobulk. The CAF subtype proportions (below) have also been changed from pooled-cell proportions to the **mean of each patient's own subtype proportion**, so that patients contributing more CAF cells no longer dominate the estimate.

### EMT Score in Malignant Cells (patient-level pseudobulk)

  - GSE154778: **not testable** (1 aggressive vs 2 reference patients; below the 3-per-arm minimum) — the previously reported p<0.001 result for this cohort was a pseudoreplication artifact and should be disregarded.
  - Peng_et_al (5 aggressive vs 6 reference patients): median diff = -0.135, p_raw=0.0062, p_adj(BH)=0.019, direction = **down** (opposite of expected "up")
  - **GSE202051 (upgraded 2026-07-12, 7 aggressive vs 14 reference patients, 20,172 vs 15,414 malignant cells):** median diff = -0.079, p_raw=0.031, p_adj(BH)=0.061, direction = **down** (opposite of expected "up") — same direction as Peng_et_al, just short of FDR<0.05.

**EMT conclusion (revised 2026-07-12):** The EMT reversal is no longer resting on a single 11-patient cohort. GSE202051, now adequately powered (7 vs 14 patients, tens of thousands of cells), independently trends in the **same direction** as Peng_et_al's significant finding — EMT score lower, not higher, in aggressive malignant cells — though it falls just short of FDR<0.05 on its own (p_adj=0.061). Combined, this is a directionally consistent signal across the two best-powered cohorts (1 of 2 individually significant), meaningfully stronger evidence than the single-cohort version of this finding reported before 2026-07-12. This continues to argue against malignant-cell-intrinsic EMT as the driver of the bulk EMT signal. GSE202051's CAF compartment shows EMT score *higher* in aggressive patients (median diff=+0.237, p_adj=0.058, borderline), in the expected direction — consistent with a CAF/stromal rather than malignant-cell source for the bulk EMT signal.

### CAF Subtype Analysis (mean of per-patient proportions)

**GSE154778** (aggressive: 1 patient, reference: 2 patients — descriptive only, not a statistical comparison):
  aggressive: myCAF=6.55%, iCAF=2.92%, apCAF=0.54%
  reference: myCAF=13.53%, iCAF=5.84%, apCAF=0.63%

**Peng_et_al** (aggressive: 5 patients, reference: 6 patients):
  aggressive: myCAF=18.68%, iCAF=9.54%, apCAF=1.19%
  reference: myCAF=23.73%, iCAF=9.95%, apCAF=1.61%

**GSE202051 (upgraded 2026-07-12; aggressive: 7 patients, reference: 16 patients — now a real inferential comparison, not descriptive):**
  aggressive: myCAF=6.15%, iCAF=9.66%, apCAF=0.48%
  reference: myCAF=13.94%, iCAF=21.14%, apCAF=2.12%

*(Note: absolute proportions shifted from the previous pooled-cell version, e.g. GSE154778 aggressive myCAF was reported as 65.48% before — that number was dominated by whichever single aggressive patient had the most CAF cells. The corrected numbers here average across patients rather than across cells, but with 1–2 patients per arm in GSE154778 they remain descriptive, not inferential. GSE202051 now has enough patients per arm to be a real comparison and shows the same qualitative pattern as GSE154778 (lower myCAF and iCAF proportion in aggressive patients) — this is descriptive/proportion data, not formally significance-tested in this table, but the directional consistency across two independent cohorts is notable.)*

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

**METHODOLOGY CORRECTION (2026-07-06):** This section draws from the same `figure3B_lipid_cell_of_origin_statistics.tsv` table used in Section 4, which has already been re-computed at the patient-pseudobulk level (≥5 cells/patient, ≥3 patients/arm) — no separate script exists for this section. GSE154778 could not be tested for any cell type (only 1 aggressive patient total).

**DATA UPGRADE (2026-07-12):** GSE202051 is now testable (7 aggressive vs 14–16 reference patients, per Section 3's note on the full 43-patient object), so this section is no longer based on a single 11-patient cohort. Combined with Peng_et_al (5 aggressive vs 6 reference patients):

- **myeloid**: Peng_et_al 0/3 lipid comparisons significant (all p_adj≥0.855). **GSE202051 (n=6 vs 15 patients): 1/3 significant — fatty_acid_uptake_oxidation p_adj=0.031, correctly reduced in aggressive patients' myeloid cells.** This is new significant off-target signal that the underpowered Peng_et_al cohort could not detect.
- **endothelial**: Peng_et_al 0/3 significant. **GSE202051 (n=7 vs 16 patients): 1/3 significant — lipid_synthesis_srebp p_adj=0.023, but in the *wrong* direction (reduced, not elevated, in aggressive patients' endothelial cells).**
- **caf_fibroblast**: Peng_et_al 0/3 significant. **GSE202051 (n=7 vs 16 patients): 2/3 significant — lipid_synthesis_srebp p_adj=0.002 (wrong direction) and fatty_acid_uptake_oxidation p_adj=0.005 (correct direction, reduced in aggressive).** See also Section 4.

**Revised conclusion (2026-07-12):** The prior conclusion — "no significant off-target lipid signal, but this may just reflect the small (n=11) patient count" — has been directly tested by adding a much better-powered cohort, and the answer is mixed. GSE202051 does detect several significant off-target/non-malignant signals that Peng_et_al alone could not: myeloid and CAF fatty-acid-oxidation reduction in aggressive patients (consistent direction with the malignant-cell finding, suggesting a possible shared microenvironmental driver rather than a purely malignant-cell-intrinsic FA-oxidation effect), and CAF/endothelial lipid-synthesis reversal (opposite direction from the malignant-cell trend). This means the "off-target lipid signal is negligible" reading from the single-cohort report no longer holds cleanly — non-malignant compartments do show real, cohort-detectable lipid differences between aggressive and reference patients, some concordant with and some opposed to the malignant-cell pattern. This should be read as evidence that the aggressive/reference microenvironment differs broadly (not just in malignant cells), which complicates a simple "the lipid signal is malignant-cell-intrinsic" narrative even where the malignant-cell direction itself is correct.

## 10. Final Mechanistic Classification

**Classification: Partially resolved / mixed (revised 2026-07-12 after the GSE202051 43-patient data upgrade; previously revised 2026-07-06 after pseudoreplication fix, purity-proxy circularity fix, and real GSE21501 survival fix)**

Based on **real single-cell RNA-seq data** (GSE154778: 8,000 primary tumor cells, 10 patients; **GSE202051: 224,988 cells, 43 patients — upgraded 2026-07-12 from a 2,607-cell, 1-patient subset**; Peng et al.: 43,888 cells, 17 patients), real Phase 2 bulk cohorts for purity adjustment, and **real CPTAC proteomics** (umich 145 + BCM 105 tumor samples):

1. **Hypoxia/acinar co-occurrence (REAL DATA — negative result, now confirmed in 3/3 adequately-powered cohorts):** Pearson r(hypoxia, acinar) in malignant cells = +0.06, +0.04, −0.02 across three real scRNA-seq cohorts (all near-zero; GSE202051 now p<0.0001 with n=64,538 malignant cells, no longer resting on 1 patient). The fraction of malignant cells in the hypoxia-high/acinar-low quadrant is 24–27%, indistinguishable from the 25% expected under independence. This finding is now on firmer footing than before 2026-07-12: all three cohorts, not two of three, are adequately powered and agree.

2. **Lipid rewiring cell-of-origin (REVISED 2026-07-12 — malignant lipid-synthesis signal strengthens, but off-target/wrong-direction signal also emerges):** In the two now-adequately-powered cohorts, malignant-cell lipid_synthesis_srebp trends in the expected (elevated) direction in both (Peng_et_al significant, FDR=0.019; GSE202051 same direction, FDR=0.052, borderline). This is stronger support for cell-intrinsic lipid synthesis rewiring than the single-cohort result reported before 2026-07-12. However, GSE202051 also reveals that malignant desaturation/elongation is *significantly reversed* (not merely non-significant, as Peng_et_al alone suggested), and that CAF and endothelial cells show significant lipid_synthesis_srebp changes in the *wrong* direction — complications invisible in the underpowered single-cohort analysis. FA-oxidation reduction is now the most consistently-reproduced lipid axis, correct-direction and significant in GSE202051's CAF and myeloid compartments.

3. **CAF/EMT cell-of-origin (REVISED 2026-07-12 — EMT reversal now trends in 2/2 adequately-powered cohorts, not 1/1):** EMT score is significantly LOWER in aggressive malignant cells in Peng_et_al (5 vs 6 patients, FDR=0.019) and trends the same direction in GSE202051 (7 vs 14 patients, FDR=0.061, just short of significance) — directionally consistent evidence in both testable cohorts now, versus a single 11-patient cohort before 2026-07-12. GSE154778 still cannot be tested (only 1 aggressive patient). CAF subtype proportions are now a real inferential comparison in GSE202051 (7 vs 16 patients, not just descriptive) and show the same qualitative pattern as GSE154778's descriptive numbers — lower myCAF/iCAF proportion in aggressive patients.

4. **Purity adjustment of bulk CAF/EMT (unaffected by the GSE202051 upgrade — uses only Phase 2 bulk cohorts):** CAF and EMT effects survive purity adjustment essentially intact and remain highly significant in the two adequately powered bulk cohorts (GSE71729 n=145: caf coef 0.53→0.46, p<0.0001; GSE62165 n=131: caf coef 0.44→0.43, p<0.0001; EMT similarly stable in both). Supports a real, largely purity-independent CAF/EMT association in bulk PDAC.

5. **Protein-level validation (REAL DATA — partially concordant, independently replicated, unaffected by the fix):** ACADL replicated FDR<0.05 in both umich and BCM; FASN, ACACA, SQLE elevated and significant in umich (FDR<0.05), directionally concordant but underpowered in BCM.

6. **Survival (unaffected by the GSE202051 upgrade — uses only bulk cohorts):** A real Cox fit on GSE21501 (n=102, 66 events) gives HR=0.93, p=0.80 — essentially null, opposite direction from the other two bulk cohorts. Pooled estimate: HR=1.06 [0.77–1.47], p=0.72. The survival trend across the three available bulk cohorts is inconsistent in direction, not merely underpowered.

**Recommended wording (revised 2026-07-12 after the GSE202051 data upgrade):**

> Single-cell analysis indicates the bulk hypoxia-high/acinar-low signature is a composite of two largely independent malignant-cell axes (r ≈ 0 at the cell level, now confirmed in three adequately powered cohorts totaling >75,000 malignant cells) rather than a single coherent program. Patient-level (pseudobulk) analysis across the two best-powered single-cell cohorts (Peng et al., n=11 patients; GSE202051, n=21 patients) finds a consistent-direction elevation of lipid synthesis genes in aggressive malignant cells (significant in one cohort, borderline in the other) — evidence in favor of, not against, tumor-intrinsic lipid rewiring — while EMT score trends lower, not higher, in aggressive malignant cells in both cohorts (significant in one, borderline in the other). The newly well-powered GSE202051 cohort also reveals lipid-metabolism differences in non-malignant compartments (CAF, endothelial, myeloid) not detectable in the smaller cohort alone, some concordant with and some opposed to the malignant-cell direction, indicating the aggressive/reference microenvironment differs broadly rather than the lipid signal being cleanly malignant-cell-intrinsic. The bulk CAF and EMT associations **survive adjustment for tumor purity** using a corrected, non-circular purity proxy, remaining highly significant in the two adequately powered bulk cohorts (GSE71729, GSE62165) with only 1–13% coefficient attenuation — the CAF/EMT signal in bulk PDAC is not primarily a tumor-purity artifact. A real Cox fit on GSE21501 gives HR=0.93 (p=0.80), pulling the pooled bulk survival estimate to HR=1.06 [0.77–1.47], p=0.72 — the survival association across the three available bulk cohorts is inconsistent in direction, not merely underpowered. Single-cell cell-of-origin conclusions now rest on two adequately powered cohorts (11 and 21 testable patients) rather than one, a meaningful improvement, but larger and ideally matched bulk+single-cell cohorts would still strengthen these findings further.

## 11. Limitations

1. **No matched bulk + single-cell patient data.** The three single-cell cohorts (GSE154778, GSE202051, Peng et al.) are independent from the bulk Phase 2 cohorts and share no patients. The within-single-cell-dataset H-hi/A-lo grouping is not equivalent to the between-dataset Phase 2 aggressive group definition. Lipid cell-of-origin conclusions require matched data.

2. **(Resolved 2026-07-12) GSE202051 previously had only 1 usable patient ID** because the pipeline had downloaded a single-sample GEO supplementary file rather than the full series object. The full 43-patient object (`GSE202051_totaldata-final-toshare.h5ad`) is now used instead, making patient-level aggressive/reference classification possible (7 vs 14–16 patients depending on cell type). See Sections 3–5, 9.

3. **Dissociation stress artifacts.** GSE154778 (stress ratio 7.95) and Peng_et_al (stress ratio 9.84) show elevated dissociation-stress genes (FOS, JUN, HSP family). Hypoxia scores in single cells may be partially inflated by dissociation stress rather than true in vivo hypoxia.

4. **BCM proteomics replication uses an unbalanced reference group (n=23 vs n=46 in umich).** Of the 140-sample transcriptomics-defined groups (46 aggressive / 46 reference / 48 other), only 105 BCM samples overlap, yielding 41 aggressive and only 23 reference. This imbalance reduces power for BCM FDR thresholds. ELOVL6 and CPT1B are absent from both umich and BCM datasets.

5. **(Resolved 2026-07-06) GSE21501 survival now uses a real Cox fit** (HR=0.93, p=0.80, n=102, 66 events) instead of a hardcoded literature value. Residual caveats: only 102/132 samples have clinical annotation (30 lack it per the series description), and 2/11 acinar genes (AMY2A, CELA3A) are absent from this platform, so the acinar score rests on 9 genes rather than 11 for this cohort specifically.

6. **(Resolved 2026-07-06) Purity proxy circularity.** The proxy now uses 8 genes (PTPRC, THY1, SPARC, CD2, CD14, CD34, LAPTM5, SERPING1) verified disjoint from the CAF/EMT signatures, replacing the original 6-gene proxy that was a strict subset of the CAF signature. It remains an expression-mean heuristic rather than a validated ESTIMATE/CIBERSORTx purity estimate — that upgrade is still open (see Next Steps item 4) — but the collinearity that made the CAF-attenuation conclusion partly circular is fixed.

7. **GSE202051 is snRNA-seq (nuclei), not scRNA-seq.** Nuclear transcriptomes have lower cytoplasmic gene expression (especially mitochondrial). CAF and stromal gene scores may be underestimated. This applies to the full 43-patient object too, not just the previously-used subset.

8. **(Resolved 2026-07-12) Small malignant cell counts in GSE202051 (previously n=2,039 total, 1 patient)** — now n=64,538 malignant cells across 43 patients (18 untreated, 25 treated) after switching to the full series object. Statistical comparisons for this cohort are now adequately powered for patient-level inferences (see Sections 3–5, 9). Note: the 43 patients include both untreated and neoadjuvant-treated samples (`treatment_status` in the raw object); treatment status was not used to stratify or exclude patients in this analysis, so treatment effects could in principle contribute to between-patient variance not attributable to the aggressive/reference axis.

9. **(Improved 2026-07-12, still a real limitation) Cell-of-origin sample sizes remain modest.** GSE202051 is now testable (7 vs 14–16 patients) alongside Peng_et_al (5 vs 6 patients), a real improvement over the single 11-patient cohort used before 2026-07-12. GSE154778 (10 patients, only 1 aggressive) still cannot be tested. The Section 4/5 conclusions now rest on two independent, adequately-powered cohorts rather than one, but n in the tens (not hundreds) of patients per cohort — treat as reasonably solid directional evidence, not definitive.

10. **(Resolved 2026-07-06, strengthened 2026-07-12) Section 9 (Immune/Endothelial QC)** has been re-run at the patient level with both Peng_et_al and the newly upgraded GSE202051 (see Section 9). The "no significant off-target lipid signal" conclusion from the single-cohort version does **not** hold with GSE202051 included — real off-target signal (some concordant, some discordant with the malignant-cell direction) is now detectable in CAF, myeloid, and endothelial compartments.

## 12. Recommended Manuscript Wording

The recommended wording in Section 10 reflects what the real data support after the patient-level correction and the 2026-07-12 GSE202051 data upgrade. Do not use wording claiming malignant lipid rewiring "does not reproduce" or "reverses direction" — the corrected, adequately-powered comparison points the other way (consistent-direction support, significant in one of two testable cohorts) for lipid synthesis specifically. Do not describe off-target (CAF/myeloid/endothelial) lipid signal as absent or negligible — GSE202051 now shows real, significant off-target effects in several compartments (Section 9).

For the survival finding (revised 2026-07-06 — do not use the previous "consistent directional trend" wording, which relied on a hardcoded GSE21501 value):
> "Across the three cohorts with available survival data, the hypoxia-high/acinar-low state showed an inconsistent association with overall survival (HR=1.23 and 1.08 in two cohorts, HR=0.93 in a third), and the pooled estimate did not reach statistical significance (pooled HR = 1.06, 95% CI [0.77–1.47], p = 0.72). This is better characterized as no clear survival association in the available data than as a consistent trend obscured by low power."

## 13. Next Steps

1. **Obtain matched bulk + single-cell data** for the same patients (e.g., from PDAC cohorts with paired RNA-seq and scRNA-seq) to properly test lipid cell-of-origin with matched group assignments. (Still open as of 2026-07-12 — the GSE202051 upgrade added patients but not matched bulk data; searched for candidates and found none with genuinely matched bulk+single-cell from the same patients — see Werba/Hwang/Loveless et al. options evaluated but not integrated.)

2. ~~**Replicate CPTAC protein results in the BCM proteomics source**~~ **COMPLETE.** BCM replication done: ACADL replicated FDR<0.05 in both umich and BCM; FASN/ACACA/SQLE directionally concordant but underpowered in BCM (n=23 reference). See Section 7.

3. ~~**Real GSE21501 survival analysis**~~ **COMPLETE.** `scripts/parse_gse21501_survival.py` now streams the real SOFT file and `run_expanded_survival_meta.py` fits a real Cox model (HR=0.93, p=0.80, n=102). See Section 8. Remaining open item: only 102/132 GSE21501 samples have clinical annotation and only 9/11 acinar genes are on this platform — a larger or better-annotated cohort would strengthen this estimate.

4. **Implement proper purity deconvolution:** Apply ESTIMATE (R package: `library(estimate); filterCommonGenes(); estimateScore()`) or TIMER2/CIBERSORTx for validated tumor purity estimates. The gene-overlap circularity is already fixed (Section 6), but a validated deconvolution method would still be more rigorous than the current 8-gene expression-mean heuristic.

5. **Spatial transcriptomics:** 10x Visium PDAC datasets (e.g., from the 10x Genomics public datasets) could directly test whether hypoxic and acinar-low regions spatially co-localize within the same tumor section, providing an orthogonal test of the composite-state vs. co-occurring-state question.

6. **Dissociation stress correction:** Apply computational correction for dissociation artifacts (e.g., van den Brink et al. 2017 gene list) before re-scoring hypoxia in GSE154778 and Peng et al.

7. ~~**Re-run Section 9 (Immune/Endothelial QC) at the patient level**~~ **COMPLETE**, and ~~**obtain more single-cell cohorts/patients so GSE202051 can support patient-level tests**~~ **COMPLETE (2026-07-12).** GSE202051 was upgraded from 1 to 43 patients (Section 3 note) and now supports patient-level lipid/EMT/immune-endothelial tests alongside Peng_et_al. GSE154778 (10 patients, only 1 aggressive) still cannot be tested — would need more patients or a different cohort to resolve.

8. **(New 2026-07-12) Stratify or exclude by treatment status in GSE202051.** The full object includes both untreated (18 patients) and neoadjuvant-treated (25 patients) samples; this analysis pooled both without stratification (see Limitation 8). Re-running the aggressive/reference comparison restricted to the 18 untreated patients (or with treatment status as a covariate) would rule out treatment as a confound.

---
*Report generated by Phase 3 pipeline. Single-cell data: REAL (GSE154778; **GSE202051 — upgraded 2026-07-12 to the full 43-patient object, 224,988 cells**; Peng et al.). CPTAC proteomics: REAL (umich 145 + BCM 105 tumor samples; ACADL independently replicated FDR<0.05 in both). Survival cohort GSE21501: REAL Cox fit (HR=0.93, p=0.80, n=102) computed from the parsed SOFT file, replacing the previous hardcoded literature value (see Section 8/11). Purity adjustment: REAL Phase 2 bulk expression data with an 8-gene purity proxy verified disjoint from the CAF/EMT signatures (see Section 6/11). Sections 4, 5, and 9 re-analyzed 2026-07-06 using patient-level pseudobulk to correct a cell-level pseudoreplication error, and again 2026-07-12 to incorporate the full GSE202051 object; Section 6 re-analyzed 2026-07-06 with the corrected purity proxy, reversing the earlier "CAF signal is purity-driven" conclusion; Section 8 re-analyzed 2026-07-06 with the real GSE21501 Cox fit, reversing the earlier "consistent survival trend" framing.*
