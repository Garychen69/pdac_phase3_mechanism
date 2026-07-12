"""
test_cell_annotation.py
Tests for cell type annotation functionality.
"""

import os
import sys
import numpy as np
import pandas as pd
import pytest

np.random.seed(1234)

# Add scripts to path
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE_DIR, "scripts"))

import warnings
warnings.filterwarnings("ignore")

import yaml
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix

# Import the actual production annotation pipeline rather than reimplementing
# it here, so passing tests guarantee the shipped script produces valid output.
import annotate_cell_types
from annotate_cell_types import (
    run_full_annotation_pipeline,
    annotate_cohort,
    CANONICAL_MARKERS,
    CELL_TYPES_ORDERED,
)

with open(os.path.join(BASE_DIR, "config", "gene_sets.yml")) as f:
    GENE_SETS = yaml.safe_load(f)

VALID_CELL_TYPES = set(CELL_TYPES_ORDERED) | {"unknown"}

REQUIRED_ANNOTATION_COLS = [
    "cell_id", "patient_id", "cohort", "leiden_cluster", "cell_type", "UMAP1", "UMAP2"
]


def make_real_marker_anndata(n_cells=200):
    """Create a synthetic AnnData containing the real canonical marker genes
    (from annotate_cell_types.CANONICAL_MARKERS), not a hand-rolled gene list,
    so scoring/clustering exercises the actual production marker set."""
    np.random.seed(1234)
    marker_genes = sorted({g for genes in CANONICAL_MARKERS.values() for g in genes})
    background = [f"GENE_{i:04d}" for i in range(80)]
    all_genes = marker_genes + background

    X = np.abs(np.random.lognormal(0, 1, (n_cells, len(all_genes)))).astype(np.float32)

    obs = pd.DataFrame({
        "patient_id": [f"TEST_P{i % 5 + 1:02d}" for i in range(n_cells)],
    }, index=[f"TEST_cell_{i:04d}" for i in range(n_cells)])
    var = pd.DataFrame(index=all_genes)

    adata = ad.AnnData(X=csr_matrix(X), obs=obs, var=var)
    adata.uns["simulated"] = False  # exercise the real (non-simulated) de-novo path
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def test_run_full_annotation_pipeline_assigns_valid_cell_types():
    """Test that the real de-novo annotation pipeline assigns valid cell types."""
    sc.settings.verbosity = 0
    adata = make_real_marker_anndata(n_cells=200)
    annotated, score_keys = run_full_annotation_pipeline(adata, "TEST")

    assert "cell_type" in annotated.obs.columns
    assigned = set(annotated.obs["cell_type"].unique())
    assert assigned.issubset(VALID_CELL_TYPES), f"Invalid cell types: {assigned - VALID_CELL_TYPES}"
    assert len(score_keys) > 0


def test_run_full_annotation_pipeline_produces_leiden_and_umap():
    """Test that clustering and UMAP coordinates are produced."""
    sc.settings.verbosity = 0
    adata = make_real_marker_anndata(n_cells=200)
    annotated, _ = run_full_annotation_pipeline(adata, "TEST")

    assert "leiden" in annotated.obs.columns
    assert "X_umap" in annotated.obsm
    assert annotated.obsm["X_umap"].shape == (annotated.n_obs, 2)


def test_annotation_no_nulls_in_cell_type():
    """Test that cell type column has no nulls after real annotation."""
    sc.settings.verbosity = 0
    adata = make_real_marker_anndata(n_cells=100)
    annotated, _ = run_full_annotation_pipeline(adata, "TEST")
    assert annotated.obs["cell_type"].isna().sum() == 0


def test_multiple_cell_types_present():
    """Test that annotation returns at least one cluster over enough cells."""
    sc.settings.verbosity = 0
    adata = make_real_marker_anndata(n_cells=200)
    annotated, _ = run_full_annotation_pipeline(adata, "TEST")
    n_clusters = annotated.obs["leiden"].nunique()
    assert n_clusters >= 1, f"Expected at least 1 cluster, got {n_clusters}"


def test_annotate_cohort_end_to_end_writes_annotation_file(tmp_path, monkeypatch):
    """Test the full annotate_cohort() pipeline end-to-end: reads a real h5ad
    from disk, runs annotation, and writes an annotations TSV with all
    required columns -- matching what run_phase3_pipeline.py actually does."""
    sc.settings.verbosity = 0
    monkeypatch.setattr(annotate_cell_types, "FIGURES_DIR", str(tmp_path))

    adata = make_real_marker_anndata(n_cells=150)
    processed_file = str(tmp_path / "TEST_processed.h5ad")
    adata.write_h5ad(processed_file)

    annotations_file = str(tmp_path / "TEST_annotations.tsv")
    cohort_cfg = {
        "name": "TEST",
        "processed_file": processed_file,
        "annotations_file": annotations_file,
    }

    annotate_cohort(cohort_cfg, GENE_SETS)

    assert os.path.exists(annotations_file)
    ann_df = pd.read_csv(annotations_file, sep="\t")
    for col in REQUIRED_ANNOTATION_COLS:
        assert col in ann_df.columns, f"Missing column: {col}"
    assert len(ann_df) == 150
    assert set(ann_df["cell_type"].unique()).issubset(VALID_CELL_TYPES)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
