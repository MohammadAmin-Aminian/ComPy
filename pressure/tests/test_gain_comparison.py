"""Display-only transforms never create an accepted broadband calibration."""

import numpy as np
import pytest

from compy_pressure.quality_review import pressure_reference_warnings
from scripts.plot_rr40_gain_comparison import band_centre_shifts, pressure_gain_db


def test_amplitude_gain_uses_twenty_log10():
    assert pressure_gain_db(1) == 0
    assert pressure_gain_db(2) == pytest.approx(6.020599913279624)
    assert pressure_gain_db(0.5) == pytest.approx(-6.020599913279624)


@pytest.mark.parametrize("gain", [0, -1, np.nan, np.inf])
def test_invalid_gain_refused(gain):
    with pytest.raises(ValueError, match="finite and positive"):
        pressure_gain_db(gain)


def test_local_p_does_not_extrapolate_or_mutate_inputs():
    f = np.array([0.01, 2, 4, 7.9, 8, 15.9, 16, 30])
    before = f.copy()
    shift = band_centre_shifts(
        f, {"4-8": {"gain_median": 0.5}, "8-16": {"gain_median": 0.25}}
    )
    assert np.isnan(shift[[0, 1, 6, 7]]).all()
    np.testing.assert_allclose(shift[[2, 3]], pressure_gain_db(0.5))
    np.testing.assert_allclose(shift[[4, 5]], pressure_gain_db(0.25))
    np.testing.assert_array_equal(f, before)


def test_overlapping_local_p_display_bands_refused():
    with pytest.raises(ValueError, match="Overlapping"):
        band_centre_shifts(
            np.array([6.0]),
            {"4-8": {"gain_median": 1}, "5-9": {"gain_median": 1}},
        )


def test_historical_reference_statistics_flagged_not_rewritten():
    old = {
        "median_in_range_fraction": 16 / 19,
        "references": {"webb": {}, "brown": {}},
        "reference_coverage": {"webb_instructions": 4},
    }
    assert pressure_reference_warnings(old) == [
        "historical_unverified_handout_included",
        "historical_brown_reference_support",
    ]
    assert old["median_in_range_fraction"] == 16 / 19
    current = {
        "references": {"brown": {"data_file": "brown_hydroacoustic_fullband.csv"}},
        "reference_coverage": {"webb_instructions": 0},
    }
    assert pressure_reference_warnings(current) == []
