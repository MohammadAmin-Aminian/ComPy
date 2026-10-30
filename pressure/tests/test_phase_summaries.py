"""Phase summaries must respect wrapping without changing bin-level QC."""

import numpy as np
import pytest

from compy_pressure.calibration_experiment import evaluate
from compy_pressure.heldout_validation import frozen_gain_check
from compy_pressure.microseism_diagnostic import diagnostic_metrics
from compy_pressure.phase import circular_median_degrees

PHASE = np.array([179.0, -179.0, 178.0, -178.0])


def test_circular_median_is_wrap_and_order_invariant():
    assert abs(circular_median_degrees(PHASE)) == pytest.approx(180)
    assert circular_median_degrees([-3, 0, 4]) == pytest.approx(0)
    assert circular_median_degrees(PHASE[::-1] + 720) == pytest.approx(
        circular_median_degrees(PHASE)
    )
    assert circular_median_degrees([]) is None
    assert circular_median_degrees([np.nan]) is None
    assert circular_median_degrees([0, 180]) is None
    with pytest.raises(ValueError, match="one-dimensional"):
        circular_median_degrees([[1, 2]])


def test_circular_median_minimizes_shortest_distance_against_brute_force():
    random = np.random.default_rng(813)
    for _ in range(40):
        sample = random.uniform(-180, 180, 31)
        actual = circular_median_degrees(sample)

        def distances(v):
            return np.abs((sample - v + 180) % 360 - 180).sum()

        assert distances(actual) <= min(map(distances, sample)) + 1e-8


def test_undefined_phase_is_not_a_fabricated_family():
    from compy_pressure.microseism_compare import phase_family

    assert phase_family(None) == ("indeterminate", None)
    assert phase_family(np.nan) == ("indeterminate", None)


def test_trial_summary_straddling_180_is_not_zero():
    spectra = {
        "coherence_squared": np.full(4, 0.99),
        "cross_pa": np.exp(1j * np.radians(PHASE)),
        "psd_ratio_pa": np.ones(4),
        "expected_pa": np.ones(4),
    }
    settings = {"coherence_squared_threshold": 0.95, "phase_limit_deg": 20}
    result = evaluate(spectra, settings, False, np.ones(4), "absolute")
    assert abs(result["median_phase_deg"]) == pytest.approx(180)
    assert result["phase_valid_bins"] == 0
    assert result["status"] == "polarity_review"  # unchanged per-bin decision


def test_heldout_summary_straddling_180_does_not_change_frozen_gain():
    spectra = {
        "coherence_squared": np.full(4, 0.99),
        "cross_pa": np.exp(1j * np.radians(PHASE)),
        "psd_ratio_pa": np.ones(4),
        "expected_pa": np.ones(4),
    }
    settings = {"coherence_squared_threshold": 0.95, "phase_limit_deg": 20}
    result = frozen_gain_check(spectra, settings, 0.67, np.ones(4))
    assert abs(result["median_phase_deg"]) == pytest.approx(180)
    assert result["status"] == "phase_review"
    assert result["frozen_gain"] == 0.67 and result["gain_refitted"] is False


def test_model_free_diagnostic_uses_wrapped_phase_summary():
    arrays = {
        "frequency_hz": np.linspace(0.11, 0.19, 4),
        "transfer_dp": np.exp(1j * np.radians(PHASE)),
        "coherence_squared": np.full(4, 0.99),
        "valid": np.ones(4, dtype=bool),
    }
    result = diagnostic_metrics(arrays, [], 0.67)["summary"]
    assert abs(result["median_nominal_transfer_phase_deg"]) == pytest.approx(180)
    assert result["fixed_pressure_multiplier"] == 0.67
