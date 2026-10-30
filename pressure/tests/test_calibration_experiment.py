"""Analytic checks for opt-in response/fit trials, separate from production defaults."""

from types import SimpleNamespace

import numpy as np
import pytest

from compy_pressure.calibration_experiment import (
    apply_response_hypothesis,
    evaluate,
    highpass_diagnostic,
    log_gain,
    pole_factor,
    real_pressure_poles,
)


def test_rc_phase_reflection_matches_causal_transfer_without_amplitude_change():
    f = np.linspace(0.03, 0.07, 9)
    omega, pole = 2j * np.pi * f, 0.1
    causal = omega / (omega + pole)
    nominal = omega / (omega - pole)
    np.testing.assert_allclose(nominal / causal, pole_factor(f, [pole]))
    np.testing.assert_allclose(abs(pole_factor(f, [pole])), 1, atol=1e-15)
    # Nominal deconvolution of a true causal signal yields causal/nominal.
    np.testing.assert_allclose((causal / nominal) * pole_factor(f, [pole]), 1)


def test_pole_trial_rejects_wrong_units_and_does_not_edit_response():
    stage = SimpleNamespace(
        poles=[0.1 + 0j, -2 + 0j], pz_transfer_function_type="LAPLACE (RADIANS/SECOND)"
    )
    response = SimpleNamespace(response_stages=[stage])
    assert real_pressure_poles(response) == [0.1]
    assert stage.poles == [0.1 + 0j, -2 + 0j]
    stage.pz_transfer_function_type = "LAPLACE (HERTZ)"
    with pytest.raises(ValueError, match="rad/s"):
        real_pressure_poles(response)


def test_log_gain_is_weighted_geometric_mean_and_invariant_to_units():
    r, t, c = (
        np.array([2.0, 4.0, 8.0]),
        np.array([1.0, 4.0, 16.0]),
        np.array([0.95, 0.98, 1.0]),
    )
    w = np.minimum(c, 0.999) / (1 - np.minimum(c, 0.999))
    expected = np.prod((t / r) ** (w / w.sum()))
    assert log_gain(r, t, c, 9) == pytest.approx(expected)
    assert log_gain(r * 1e6, t * 1e6, c, 9) == pytest.approx(expected)


def fixture_spectra(phase=0):
    f = np.linspace(0.03, 0.07, 9)
    return {
        "frequency_hz": f,
        "coherence_squared": np.full(9, 0.99),
        "cross_pa": np.exp(1j * np.radians(phase)) * np.ones(9),
        "psd_ratio_pa": np.full(9, 2.0),
        "expected_pa": np.ones(9),
        "segment_count": 9,
    }


SETTINGS = {
    "coherence_squared_threshold": 0.95,
    "phase_limit_deg": 20,
    "gain_grid_bounds": [0.1, 3.0],
}


@pytest.mark.parametrize("method", ["absolute", "weighted_log"])
def test_reflection_rechecks_phase_and_retains_overlap_coherence_guards(method):
    spectra = fixture_spectra()
    factor = pole_factor(spectra["frequency_hz"], [0.1])
    spectra["cross_pa"] = 1 / factor
    assert (
        evaluate(spectra, SETTINGS, False, np.ones(9), method)["status"]
        == "polarity_review"
    )
    fixed = evaluate(spectra, SETTINGS, False, factor, method)
    assert fixed["status"] == "candidate"
    assert fixed["pressure_multiplier"] == pytest.approx(0.5)
    assert (
        evaluate(spectra, SETTINGS, True, factor, method)["status"]
        == "arrival_overlap_review"
    )
    spectra["coherence_squared"][0:5] = 0.5
    assert evaluate(spectra, SETTINGS, False, factor, method)["status"] == "rejected_qc"


def test_joint_highpass_recovers_analytic_response_but_is_diagnostic_only():
    spectra = fixture_spectra()
    f = spectra["frequency_hz"]
    h = 0.28 * (1j * f / 0.08) / (1 + 1j * f / 0.08)
    spectra["psd_ratio_pa"] = abs(h)
    spectra["cross_pa"] = h
    result = highpass_diagnostic(spectra, np.ones(9), 0.95)
    assert result["diagnostic_only"] is True
    assert result["transfer_gain_G"] == pytest.approx(0.28, rel=1e-6)
    assert result["corner_hz"] == pytest.approx(0.08, rel=1e-6)
    assert result["phase_rms_deg"] < 1e-5
    assert (
        evaluate(spectra, SETTINGS, False, np.ones(9), "weighted_log")[
            "pressure_multiplier"
        ]
        is None
    )


def test_out_of_bounds_log_fit_does_not_relax_search_bounds():
    spectra = fixture_spectra()
    spectra["psd_ratio_pa"][:] = 0.1
    result = evaluate(spectra, SETTINGS, False, np.ones(9), "weighted_log")
    assert result["status"] == "unresolved"
    assert result["pressure_multiplier"] is None


def test_response_hypothesis_transform_is_exact_identity_and_keeps_coherence():
    spectra = fixture_spectra(phase=10)
    before = {k: np.copy(v) for k, v in spectra.items()}
    same = apply_response_hypothesis(spectra, 1.0, 1.0)
    for key in ("cross_pa", "psd_ratio_pa", "coherence_squared", "expected_pa"):
        assert np.array_equal(same[key], spectra[key])
    rp = 2 * np.exp(1j * np.radians(30))
    rz = 0.5 * np.exp(1j * np.radians(10))
    moved = apply_response_hypothesis(spectra, rp, rz)
    np.testing.assert_allclose(moved["phase_deg"], 30)
    np.testing.assert_allclose(moved["psd_ratio_pa"], 8.0)
    assert np.array_equal(moved["coherence_squared"], spectra["coherence_squared"])
    for key, value in before.items():
        assert np.array_equal(spectra[key], value)
    for bad in (0.0, np.nan):
        with pytest.raises(ValueError, match="finite and nonzero"):
            apply_response_hypothesis(spectra, bad, 1.0)


def test_response_hypothesis_path_reproduces_the_pole_reflection_trial_path():
    spectra = fixture_spectra()
    factor = pole_factor(spectra["frequency_hz"], [0.1])
    spectra["cross_pa"] = 1 / factor
    trial = evaluate(spectra, SETTINGS, False, factor, "absolute")
    moved = apply_response_hypothesis(spectra, factor, 1.0)
    mine = evaluate(moved, SETTINGS, False, np.ones(9), "absolute")
    assert mine["status"] == trial["status"] == "candidate"
    assert mine["pressure_multiplier"] == pytest.approx(
        trial["pressure_multiplier"], rel=0, abs=1e-12
    )
    # An amplitude-changing vertical hypothesis rescales the fitted gain.
    halved = apply_response_hypothesis(spectra, factor, 2.0)
    assert evaluate(halved, SETTINGS, False, np.ones(9), "absolute")[
        "pressure_multiplier"
    ] == pytest.approx(2 * trial["pressure_multiplier"])
