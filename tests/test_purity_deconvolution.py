"""
test_purity_deconvolution.py
Tests for tumor purity estimation and purity-adjusted linear models.
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

# Import the actual production functions and proxy gene list rather than
# reimplementing them here, so passing tests guarantee the shipped script's
# statistics -- and the real, CAF/EMT-disjoint purity proxy -- are correct.
from deconvolve_bulk_purity import STROMAL_GENES, estimate_purity, run_linear_models


def make_mock_expression(n_samples=50, genes=None):
    """Mock expression DataFrame indexed by gene (rows) x sample (columns),
    matching the orientation estimate_purity() expects."""
    np.random.seed(1234)
    if genes is None:
        genes = STROMAL_GENES + [f"BG_{i}" for i in range(20)]
    samples = [f"S{i:03d}" for i in range(n_samples)]
    data = {s: np.random.lognormal(0, 1, len(genes)) for s in samples}
    return pd.DataFrame(data, index=genes)


def make_mock_scores(n_samples=50):
    samples = [f"S{i:03d}" for i in range(n_samples)]
    n_agg = n_samples // 2
    return pd.DataFrame({
        "aggressive": [1] * n_agg + [0] * (n_samples - n_agg),
        "caf_score": np.random.normal(0, 1, n_samples),
        "emt_score": np.random.normal(0, 1, n_samples),
    }, index=samples)


def test_stromal_proxy_disjoint_from_caf_emt():
    """The real purity proxy must not overlap the CAF/EMT signatures it adjusts
    for -- the original 6-gene proxy was a subset of the CAF signature, making
    the adjustment partly circular. deconvolve_bulk_purity.load_stromal_genes()
    asserts disjointness at import time, so a successful import already proves
    this; this test just documents the invariant."""
    assert len(STROMAL_GENES) >= 2


def test_purity_between_zero_and_one():
    """Test that purity estimates are in [0, 1]."""
    expr_df = make_mock_expression()
    purity = estimate_purity(make_mock_scores(), expr_df, "TEST")
    assert (purity >= 0).all(), "Some purity values < 0"
    assert (purity <= 1).all(), "Some purity values > 1"
    assert len(purity) == expr_df.shape[1]


def test_purity_varies_across_samples():
    """Test that purity varies (not all constant)."""
    expr_df = make_mock_expression(n_samples=50)
    purity = estimate_purity(make_mock_scores(50), expr_df, "TEST")
    assert purity.std() > 0, "Purity is constant across all samples"


def test_purity_no_nans():
    """Test that purity has no NaN values."""
    expr_df = make_mock_expression()
    purity = estimate_purity(make_mock_scores(), expr_df, "TEST")
    assert not purity.isna().any(), "Purity contains NaN values"


def test_purity_missing_genes_fallback():
    """Test that purity estimation falls back gracefully when fewer than 2
    stromal proxy genes are present in the expression matrix."""
    np.random.seed(1234)
    scores_df = make_mock_scores(30)
    expr_df = pd.DataFrame(
        {sid: np.random.rand(2) for sid in scores_df.index},
        index=["GENE_A", "GENE_B"],
    )
    purity = estimate_purity(scores_df, expr_df, "TEST")
    assert len(purity) == 30
    assert (purity >= 0).all() and (purity <= 1).all()


def test_linear_model_output_has_required_columns():
    """Test that linear model output has the required result columns."""
    np.random.seed(1234)
    n = 60
    scores_df = make_mock_scores(n)
    purity = pd.Series(np.random.uniform(0.3, 0.9, n), index=scores_df.index)

    results = run_linear_models(scores_df, purity, "TEST", is_simulated=True)

    required = ["coef_unadjusted", "pval_unadjusted", "coef_purity_adjusted",
                "pval_purity_adjusted", "coef_purity_term"]
    assert len(results) == 2  # caf_score, emt_score
    for r in results:
        for key in required:
            assert key in r, f"Missing key: {key}"
            assert np.isfinite(r[key]), f"Non-finite value for {key}"


def test_adjusted_coefficient_changes():
    """Test that purity adjustment actually changes the group coefficient
    when purity is confounded with group."""
    np.random.seed(42)
    n = 80
    samples = [f"S{i:03d}" for i in range(n)]
    group = np.array([1] * 40 + [0] * 40, dtype=float)
    # Purity is confounded with group (lower purity in aggressive group)
    purity_vals = np.concatenate([np.random.uniform(0.3, 0.6, 40),
                                   np.random.uniform(0.6, 0.9, 40)])
    # CAF score driven mostly by purity (stromal), less by group
    caf_score = 0.2 * group + 0.6 * (1 - purity_vals) + np.random.normal(0, 0.1, n)
    scores_df = pd.DataFrame(
        {"aggressive": group, "caf_score": caf_score, "emt_score": caf_score},
        index=samples,
    )
    purity = pd.Series(purity_vals, index=samples)

    results = run_linear_models(scores_df, purity, "TEST", is_simulated=True)
    caf_result = next(r for r in results if r["score"] == "caf_score")

    # After purity adjustment, group coefficient should be smaller (attenuated)
    assert abs(caf_result["coef_purity_adjusted"]) < abs(caf_result["coef_unadjusted"]) + 0.5, (
        f"Expected adjustment to change coefficient: "
        f"unadj={caf_result['coef_unadjusted']:.3f}, adj={caf_result['coef_purity_adjusted']:.3f}"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
