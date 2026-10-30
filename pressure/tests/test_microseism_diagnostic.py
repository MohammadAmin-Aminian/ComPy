"""Tests for the model-free saved-spectrum microseism diagnostic."""

import numpy as np
import pytest

from compy_pressure import microseism_diagnostic as diagnostic


def synthetic_arrays():
    f = np.array([0.05, 0.1, 0.15, 0.2, 0.25])
    pressure = np.ones((4, 5))
    transfer = np.array(
        [
            [1, 2, 3, 4, 5],
            [1, 2, 3, 4, 5],
            [1, 4, 6, 8, 5],
            [1, 4, 6, 8, 5],
        ],
        dtype=complex,
    )
    displacement = abs(transfer) ** 2
    cross = transfer.copy()
    mean_cross = cross.mean(axis=0)
    mean_pressure = pressure.mean(axis=0)
    mean_displacement = displacement.mean(axis=0)
    coherence = abs(mean_cross) ** 2 / (mean_pressure * mean_displacement)
    return {
        "frequency_hz": f,
        "pressure_psd": mean_pressure,
        "displacement_psd": mean_displacement,
        "cross_dp": mean_cross,
        "transfer_dp": mean_cross / mean_pressure,
        "coherence_squared": coherence,
        "valid": np.ones(5, dtype=bool),
        "pressure_window_psd": pressure,
        "displacement_window_psd": displacement,
        "cross_window_dp": cross,
        "retained_window_starts": np.array(
            [
                "2020-01-01T00:00:00Z",
                "2020-01-01T01:00:00Z",
                "2020-01-02T00:00:00Z",
                "2020-01-02T01:00:00Z",
            ]
        ),
    }


def test_daily_complex_averaging_and_published_reference_screen():
    arrays = synthetic_arrays()
    daily = diagnostic.daily_products(arrays)
    assert len(daily) == 2
    np.testing.assert_allclose(daily[0]["transfer"][1:4], [2, 3, 4])
    np.testing.assert_allclose(daily[1]["transfer"][1:4], [4, 6, 8])
    assert all(item["passes_saikia_reference"] for item in daily)


def test_fixed_multiplier_is_only_reciprocal_transfer_comparison():
    arrays = synthetic_arrays()
    daily = diagnostic.daily_products(arrays)
    metrics = diagnostic.diagnostic_metrics(arrays, daily, pressure_multiplier=0.5)
    np.testing.assert_allclose(metrics["fixed_transfer"], arrays["transfer_dp"] / 0.5)
    summary = metrics["summary"]
    assert summary["median_nominal_transfer_amplitude_m_pa"] == pytest.approx(4.5)
    assert summary["median_fixed_reference_transfer_amplitude_m_pa"] == pytest.approx(9)
    with pytest.raises(ValueError, match="positive"):
        diagnostic.diagnostic_metrics(arrays, daily, pressure_multiplier=0)


def test_saved_aggregate_reproduction_and_mutation_refusal(tmp_path):
    arrays = synthetic_arrays()
    path = tmp_path / "spectra.npz"
    np.savez_compressed(path, **arrays)
    loaded = diagnostic.load_and_check_spectra(path, expected_segments=4)
    np.testing.assert_array_equal(loaded["transfer_dp"], arrays["transfer_dp"])

    arrays["transfer_dp"] = arrays["transfer_dp"].copy()
    arrays["transfer_dp"][2] += 1
    np.savez_compressed(path, **arrays)
    with pytest.raises(ValueError, match="does not reproduce"):
        diagnostic.load_and_check_spectra(path, expected_segments=4)
