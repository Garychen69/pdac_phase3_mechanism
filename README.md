# PDAC Phase 3: Single-Cell Mechanism Analysis

Phase 3 of the PDAC hypoxia-high/acinar-low aggressive state study.

## Goals

1. Single-cell resolution: which cells drive lipid rewiring and CAF/EMT signals
2. Test if hypoxia-high AND acinar-low co-occur in the SAME malignant cells
3. Tumor purity adjustment of bulk CAF/EMT signals
4. Protein-level validation via CPTAC-PDA
5. Expanded survival meta-analysis

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run full pipeline
python scripts/run_phase3_pipeline.py

# Or run individual steps
python scripts/download_singlecell_data.py
python scripts/preprocess_singlecell_cohorts.py
python scripts/annotate_cell_types.py
python scripts/score_singlecell_signatures.py
python scripts/analyze_hypoxia_acinar_cooccurrence.py
python scripts/analyze_lipid_cell_of_origin.py
python scripts/analyze_caf_emt_cell_of_origin.py
python scripts/deconvolve_bulk_purity.py
python scripts/analyze_cptac_protein.py
python scripts/parse_gse21501_survival.py  # must run before run_expanded_survival_meta.py
python scripts/run_expanded_survival_meta.py
python scripts/generate_phase3_report.py

# Run tests
pytest tests/ -v
```

## Data Status

All datasets use **REAL data**, obtained outside the automated GEOparse download
(which cannot fetch these count matrices directly — see `data/raw/singlecell/DOWNLOAD_LOG.md`):

- **GSE154778**: real (30MB DGE CSV, ~8,000 primary tumor cells, 10 patients)
- **GSE202051**: real (h5ad, ~2,600 snRNA-seq nuclei, 1 patient, pre-annotated)
- **Peng et al. (GSE155698)**: real (919MB RAW.tar, ~44,000 cells, 17 patients)
- **CPTAC-PDA**: real, via the `cptac` Python package (`cptac.Pdac()`) — umich + BCM proteomics, WashU transcriptomics
- **GSE21501** (expanded survival cohort): real, SOFT file parsed by `scripts/parse_gse21501_survival.py`

A simulated fallback path still exists in the code (`simulate_anndata()` in
`preprocess_singlecell_cohorts.py`, `_simulated_fallback()` in
`run_expanded_survival_meta.py`) and only activates automatically if the real
files above are missing — check `PHASE3_MECHANISM_REPORT.md`'s "Dataset Summary"
section for the per-cohort REAL/SIMULATED status of any given run, not this file.

## Key Output Files

- `results/reports/PHASE3_MECHANISM_REPORT.md` — main findings report
- `results/tables/figure3A_hypoxia_acinar_cooccurrence_by_cohort.tsv`
- `results/tables/figure3B_lipid_cell_of_origin_statistics.tsv`
- `results/tables/figure3CD_caf_emt_cell_of_origin_statistics.tsv`
- `results/tables/figure3E_purity_adjusted_caf_emt_results_by_cohort.tsv`
- `results/tables/figure3F_cptac_lipid_protein_statistics.tsv`
- `results/tables/figure3G_expanded_survival_meta_analysis.tsv`

## Phase 2 Context

Phase 2 validated lipid rewiring, CAF/EMT elevation, and consistent HR > 1 across
GSE79668, GSE71729, GSE62165. Phase 3 zooms in on mechanisms using single-cell data.
