"""
test_signature_scoring_singlecell.py
Tests for single-cell signature scoring.
"""

import os
import sys
import numpy as np
import pandas as pd
import pytest
import warnings
warnings.filterwarnings("ignore")

np.random.seed(1234)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE_DIR, "scripts"))

import yaml
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix

# Import the actual production scoring function and signature list rather than
# reimplementing them here, so passing tests guarantee the shipped script
# (and the real gene_sets.yml marker genes) produce correct output.
from score_singlecell_signatures import SIGNATURES, score_cohort

with open(os.path.join(BASE_DIR, "config", "gene_sets.yml")) as f:
    GENE_SETS = yaml.safe_load(f)

REQUIRED_SCORE_COLS = [
    "cell_id", "patient_id", "cohort", "cell_type",
    "hypoxia_score", "acinar_identity_score",
    "lipid_synthesis_srebp_score", "desaturation_elongation_score",
    "fatty_acid_uptake_oxidation_score", "caf_score", "emt_score",
]


def make_test_adata(n_cells=150, include_signature_genes=True):
    """Create a test AnnData containing genes from the real gene_sets.yml
    signatures (not a hand-rolled marker list), so scoring exercises the
    actual production gene sets."""
    np.random.seed(1234)
    all_sig_genes = []
    if include_signature_genes:
        for sig_name in SIGNATURES:
            all_sig_genes.extend(GENE_SETS.get(sig_name, []))
    all_sig_genes = list(set(all_sig_genes))

    background = [f"BGENE_{i:04d}" for i in range(max(50, 100 - len(all_sig_genes)))]
    all_genes = list(set(all_sig_genes + background))
    n_genes = len(all_genes)

    X = np.abs(np.random.lognormal(0, 1, (n_cells, n_genes))).astype(np.float32)
    obs = pd.DataFrame({
        "patient_id": [f"TP{i % 5}" for i in range(n_cells)],
        "cell_type": np.random.choice(["malignant_epithelial", "caf_fibroblast", "myeloid"], n_cells),
    }, index=[f"T_cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=all_genes)
    adata = ad.AnnData(X=csr_matrix(X), obs=obs, var=var)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def run_score_cohort(tmp_path, n_cells=150):
    """Write a synthetic AnnData to disk and run it through the real
    score_cohort() pipeline, exactly as run_phase3_pipeline.py does."""
    sc.settings.verbosity = 0
    adata = make_test_adata(n_cells=n_cells)
    processed_file = str(tmp_path / "TEST_processed.h5ad")
    adata.write_h5ad(processed_file)

    cohort = {
        "name": "TEST",
        "processed_file": processed_file,
        "annotations_file": str(tmp_path / "TEST_annotations_does_not_exist.tsv"),
        "scores_file": str(tmp_path / "TEST_scores.tsv"),
    }
    return score_cohort(cohort, GENE_SETS)


def test_all_seven_signatures_scoreable(tmp_path):
    """Test that all 7 required signatures are scored by the production script."""
    out_df = run_score_cohort(tmp_path)
    score_cols = [c for c in out_df.columns if c.endswith("_score")]
    assert len(score_cols) == 7, f"Expected 7 signature scores, got {len(score_cols)}: {score_cols}"
    assert set(SIGNATURES) == {c.replace("_score", "") for c in score_cols}


def test_output_dataframe_has_required_columns(tmp_path):
    """Test that the output scores dataframe has all required columns."""
    out_df = run_score_cohort(tmp_path)
    for col in REQUIRED_SCORE_COLS:
        assert col in out_df.columns, f"Missing required column: {col}"


def test_scores_produce_finite_values(tmp_path):
    """Test that score_genes produces finite numeric scores for every signature."""
    out_df = run_score_cohort(tmp_path)
    for sig_name in SIGNATURES:
        col = f"{sig_name}_score"
        assert np.all(np.isfinite(out_df[col].values)), f"{col} contains non-finite values"


def test_scores_are_not_all_identical(tmp_path):
    """Test that scores vary across cells (not all same value) for a signature
    with multiple genes present."""
    out_df = run_score_cohort(tmp_path, n_cells=100)
    assert out_df["lipid_synthesis_srebp_score"].std() > 0, \
        "All lipid_synthesis_srebp scores are identical — scoring may have failed"


def test_scores_file_written_to_disk(tmp_path):
    """Test that score_cohort persists the scores TSV, matching the real pipeline."""
    out_df = run_score_cohort(tmp_path)
    scores_path = tmp_path / "TEST_scores.tsv"
    assert scores_path.exists()
    on_disk = pd.read_csv(scores_path, sep="\t")
    assert len(on_disk) == len(out_df)


def test_missing_signature_genes_assigned_zero(tmp_path):
    """Test that a signature with zero genes present falls back to a zero score
    (production behavior) rather than erroring."""
    sc.settings.verbosity = 0
    adata = make_test_adata(n_cells=50, include_signature_genes=False)
    processed_file = str(tmp_path / "TEST_nogenes_processed.h5ad")
    adata.write_h5ad(processed_file)
    cohort = {
        "name": "TEST_NOGENES",
        "processed_file": processed_file,
        "annotations_file": str(tmp_path / "does_not_exist.tsv"),
        "scores_file": str(tmp_path / "TEST_nogenes_scores.tsv"),
    }
    out_df = score_cohort(cohort, GENE_SETS)
    for sig_name in SIGNATURES:
        assert (out_df[f"{sig_name}_score"] == 0.0).all()


def test_score_direction_with_boosted_genes():
    """Test that cells with high marker expression score higher (sanity check
    on scanpy's score_genes, which the production code relies on)."""
    sc.settings.verbosity = 0
    np.random.seed(1234)
    n_cells = 100
    sig_genes = ["GENE_A", "GENE_B", "GENE_C"]
    background = [f"BG_{i}" for i in range(50)]
    all_genes = sig_genes + background
    n_genes = len(all_genes)
    gene_idx = {g: i for i, g in enumerate(all_genes)}

    X = np.abs(np.random.lognormal(0, 0.5, (n_cells, n_genes))).astype(np.float32)
    # Boost signature genes for first 50 cells
    for g in sig_genes:
        X[:50, gene_idx[g]] *= 10

    obs = pd.DataFrame(index=[f"c{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=all_genes)
    adata = ad.AnnData(X=csr_matrix(X), obs=obs, var=var)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.score_genes(adata, sig_genes, score_name="test_sig", random_state=1234)
    high_mean = adata.obs["test_sig"].values[:50].mean()
    low_mean = adata.obs["test_sig"].values[50:].mean()
    assert high_mean > low_mean, f"Expected higher scores for boosted cells: {high_mean:.3f} vs {low_mean:.3f}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
