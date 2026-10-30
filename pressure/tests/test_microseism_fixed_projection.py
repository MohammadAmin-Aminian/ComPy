import numpy as np
import pytest

from compy_pressure.microseism_fixed_projection import (
    ideal_rayleigh_distortion,
    prior_projections,
    projection_products,
    validate_weights,
)


def saved_fit(date, event="event", status="rejected_qc", weights=(0, 0, 1)):
    return {
        "id": event,
        "quiet_start_utc": date,
        "status": status,
        "rotation_coefficients": list(weights),
        "tilt": {"method": "tiskitpy.CleanRotator", "response_output": "VEL"},
    }


def test_primary_is_latest_completed_prior_not_best_event_qc():
    events = [
        saved_fit("2013-04-16", "good", "candidate"),
        saved_fit("2013-04-18", "bad", "rejected_qc"),
        saved_fit("2013-07-15", "future", "candidate"),
    ]
    selected = prior_projections(events, "2013-06-10")
    assert selected[-1]["event_id"] == "bad"
    assert selected[-1]["role"] == "primary_latest_prior"
    assert selected[0]["role"] == "earlier_sensitivity"
    assert len(selected) == 2
    assert prior_projections(list(reversed(events)), "2013-06-10") == selected


def test_training_overlap_duplicate_and_missing_projection_are_refused():
    with pytest.raises(ValueError, match="nonoverlapping"):
        prior_projections([saved_fit("2013-06-10T00:00:00")], "2013-06-10T01:00:00")
    with pytest.raises(ValueError, match="duplicate"):
        prior_projections([saved_fit("2013-04-18")] * 2, "2013-06-10")
    assert prior_projections([saved_fit("2013-06-10T00:00:00")], "2013-06-10T02:00:00")


@pytest.mark.parametrize(
    "weights", [[1, 0], [0, 0, -1], [0.5, 0, np.sqrt(0.75)], [np.nan, 0, 1], [0, 0, 2]]
)
def test_invalid_projection_cannot_be_silently_normalized_or_reversed(weights):
    with pytest.raises(ValueError):
        validate_weights(weights)


def random_complex(rng, shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def test_native_count_training_weights_cannot_be_used_on_physical_test_traces():
    event = saved_fit("2013-04-18")
    event["tilt"]["response_output"] = "COUNTS"
    with pytest.raises(ValueError, match="physical velocity"):
        prior_projections([event], "2013-06-10")


def test_identity_projection_preserves_input_and_complex_products():
    rng = np.random.default_rng(92)
    x = random_complex(rng, (23, 4, 10))
    before = x.copy()
    result = projection_products(x, [0, 0, 1])
    np.testing.assert_array_equal(x, before)
    np.testing.assert_allclose(result["complex_transfer_ratio"], 1)
    np.testing.assert_allclose(result["vertical_power_ratio"], 1)
    np.testing.assert_allclose(result["projected_transfer"], result["nominal_transfer"])


def test_fixed_weights_match_direct_time_series_linear_combination():
    rng = np.random.default_rng(93)
    raw = rng.normal(size=(12, 4, 64))
    weights = np.array([-0.002, -0.007, 1])
    weights /= np.linalg.norm(weights)
    x = np.fft.rfft(raw, axis=2)[:, :, 1:]
    r = projection_products(x, weights)
    projected = np.einsum("i,wit->wt", weights, raw[:, 1:])
    z = np.fft.rfft(projected, axis=1)[:, 1:]
    expected = (z * x[:, 0].conj()).mean(0) / (abs(x[:, 0]) ** 2).mean(0)
    np.testing.assert_allclose(r["projected_transfer"], expected, rtol=1e-12)


def test_velocity_and_displacement_application_commute():
    rng = np.random.default_rng(94)
    x = random_complex(rng, (23, 4, 10))
    w = [-0.002, -0.007, np.sqrt(1 - 0.002**2 - 0.007**2)]
    velocity = x.copy()
    velocity[:, 1:] *= 1j * np.linspace(0.1, 1, 10)
    disp = projection_products(x, w)
    vel = projection_products(velocity, w)
    np.testing.assert_allclose(
        disp["complex_transfer_ratio"], vel["complex_transfer_ratio"], rtol=1e-12
    )


def test_synthetic_leakage_recovery_does_not_require_pressure_as_regressor():
    rng = np.random.default_rng(95)
    clean = random_complex(rng, (23, 4, 10))
    clean[:, 3] = 2 * clean[:, 0]
    contaminated = clean.copy()
    a, b = 0.003, 0.006
    contaminated[:, 3] += a * clean[:, 1] + b * clean[:, 2]
    w = np.array([-a, -b, 1])
    w /= np.linalg.norm(w)
    result = projection_products(contaminated, w)
    np.testing.assert_allclose(result["projected_transfer"], 2 * w[2], rtol=1e-12)


def test_ideal_quadrature_injection_is_not_a_bound_on_general_signals():
    w = np.array([-0.003, -0.007, 1])
    w /= np.linalg.norm(w)
    ideal = ideal_rayleigh_distortion(w, [5])[0]
    assert ideal["maximum_amplitude_change_pct"] < 0.1
    # An in-phase physical horizontal signal can instead suffer first-order loss.
    physical = np.ones((6, 4, 10), complex)
    physical[:, 1] = 5
    physical[:, 2] = 5
    result = projection_products(physical, w)
    assert np.min(abs(result["complex_transfer_ratio"]) - 1) < -0.04


def test_too_few_windows_are_refused():
    with pytest.raises(ValueError, match=">=6"):
        projection_products(np.ones((5, 4, 10)), [0, 0, 1])
