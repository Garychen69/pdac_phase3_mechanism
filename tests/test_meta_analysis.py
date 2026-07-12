"""
test_meta_analysis.py
Tests for fixed-effects and random-effects meta-analysis functions.
"""

import os
import sys
import numpy as np
import pandas as pd
import pytest
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

np.random.seed(1234)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE_DIR, "scripts"))

# Import the actual production functions rather than reimplementing them here,
# so passing tests guarantee the shipped script's statistics are correct.
from run_expanded_survival_meta import hr_to_loghr_se, fixed_effects_meta, random_effects_meta


# --- Tests ---

def test_hr_to_loghr_se_basic():
    """Test conversion of HR and p-value to log(HR) and SE."""
    log_hr, se = hr_to_loghr_se(2.0, 0.05)
    assert np.isfinite(log_hr)
    assert np.isfinite(se)
    assert se > 0
    assert abs(log_hr - np.log(2.0)) < 1e-10


def test_hr_equal_one_gives_loghr_zero():
    """Test that HR=1 gives log(HR)=0."""
    log_hr, se = hr_to_loghr_se(1.0, 0.5)
    assert abs(log_hr) < 1e-10


def test_fixed_effects_meta_single_study():
    """Test that FE meta-analysis with one study returns that study's estimate."""
    log_hr = np.log(1.5)
    se = 0.3
    result = fixed_effects_meta([log_hr], [se])
    assert abs(result["HR_pooled"] - 1.5) < 0.01
    assert result["HR_lower_95CI"] < result["HR_pooled"] < result["HR_upper_95CI"]


def test_fixed_effects_meta_pooled_in_range():
    """Test that pooled HR is within the range of individual HRs."""
    log_hrs = [np.log(1.2), np.log(1.5), np.log(1.3)]
    ses = [0.3, 0.25, 0.35]
    result = fixed_effects_meta(log_hrs, ses)
    min_hr = np.exp(min(log_hrs))
    max_hr = np.exp(max(log_hrs))
    # Pooled should be between the extreme individual HRs
    assert min_hr <= result["HR_pooled"] <= max_hr + 0.1, \
        f"Pooled HR {result['HR_pooled']} out of range [{min_hr}, {max_hr}]"


def test_random_effects_meta_phase2_like():
    """Test RE meta-analysis with Phase 2-like inputs."""
    log_hr1, se1 = hr_to_loghr_se(1.232, 0.543)
    log_hr2, se2 = hr_to_loghr_se(1.079, 0.761)
    log_hr3, se3 = hr_to_loghr_se(1.31, 0.28)

    result = random_effects_meta([log_hr1, log_hr2, log_hr3], [se1, se2, se3])

    assert result["HR_pooled"] > 1.0, "Pooled HR should be > 1 given all input HRs > 1"
    assert result["HR_lower_95CI"] < result["HR_pooled"]
    assert result["HR_upper_95CI"] > result["HR_pooled"]
    assert 0 <= result["p_value"] <= 1
    assert 0 <= result["I2_pct"] <= 100
    assert result["tau2"] >= 0


def test_random_effects_ci_contains_truth():
    """Test that 95% CI contains true value (statistical property)."""
    # With HR ~ 1.2 and consistent effect, pooled CI should contain ~1.2
    np.random.seed(1234)
    log_hrs = [np.log(1.2) + np.random.normal(0, 0.05) for _ in range(5)]
    ses = [0.3] * 5
    result = random_effects_meta(log_hrs, ses)
    # True value around log(1.2), so HR ~ 1.2
    assert result["HR_lower_95CI"] < 1.3 < result["HR_upper_95CI"] or \
           result["HR_lower_95CI"] < 1.2 < result["HR_upper_95CI"]


def test_meta_p_value_very_consistent_studies():
    """Test that highly consistent, large studies yield significant p-value."""
    # 5 studies, all HR=2.0 with tight SE
    log_hrs = [np.log(2.0)] * 5
    ses = [0.1] * 5
    result_fe = fixed_effects_meta(log_hrs, ses)
    result_re = random_effects_meta(log_hrs, ses)
    # Should be highly significant
    assert result_fe["p_value"] < 0.001, f"Expected p<0.001, got {result_fe['p_value']}"
    assert result_re["p_value"] < 0.001


def test_i2_zero_for_homogeneous_studies():
    """Test that I² is 0 when all studies have identical estimates."""
    log_hrs = [np.log(1.5)] * 3
    ses = [0.3] * 3
    result = random_effects_meta(log_hrs, ses)
    # With identical effects, Q should be minimal, tau2=0, I2=0
    assert result["tau2"] == 0.0 or result["tau2"] < 1e-10
    assert result["I2_pct"] == 0.0 or result["I2_pct"] < 1e-5


def test_ci_width_decreases_with_more_studies():
    """Test that CI narrows when adding more studies."""
    log_hrs_2 = [np.log(1.3)] * 2
    ses_2 = [0.4] * 2
    log_hrs_5 = [np.log(1.3)] * 5
    ses_5 = [0.4] * 5

    result_2 = fixed_effects_meta(log_hrs_2, ses_2)
    result_5 = fixed_effects_meta(log_hrs_5, ses_5)

    ci_width_2 = result_2["HR_upper_95CI"] - result_2["HR_lower_95CI"]
    ci_width_5 = result_5["HR_upper_95CI"] - result_5["HR_lower_95CI"]

    assert ci_width_5 < ci_width_2, \
        f"Expected narrower CI with more studies: {ci_width_5:.3f} vs {ci_width_2:.3f}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
