"""Behavioral regression tests, importing the real public modules."""

import numpy as np
import pytest

import compy
import ffplot
import inv_compy
import Pressure_calibration as calibration
from compy_numerics import log_likelihood


def test_gaussian_likelihood_uses_squared_normalized_residual():
    d, m, s = np.array([1.0, 2.0]), np.zeros(2), np.array([1.0, 2.0])
    expected = -0.5 * np.sum(((d - m) / s) ** 2)
    assert log_likelihood(d, m, s) == expected
    assert np.isclose(inv_compy.liklihood(d, m, s=s), np.exp(expected))


@pytest.mark.parametrize("s", [0, -1, np.nan, [1, 0], [1, 2, 3]])
def test_invalid_uncertainty_rejected(s):
    with pytest.raises(ValueError):
        inv_compy.liklihood([1, 2], [0, 0], s=s)


def test_log_likelihood_remains_usable_when_probability_underflows():
    assert np.isfinite(log_likelihood([1e3], [0], 1))
    assert inv_compy.liklihood([1e3], [0], s=1) == 0


def test_roughness_energy():
    assert inv_compy.Roughness(np.ones(20), 1) == 0
    assert inv_compy.Roughness(np.linspace(0, 1, 20) ** 2, 2) > 0
    assert inv_compy.Roughness([1], 2) == 0
    with pytest.raises(ValueError):
        inv_compy.Roughness([1, 2], 1.5)


@pytest.mark.parametrize("depth", [1.0, 4000.0, 10000.0])
def test_wavenumber_satisfies_dispersion_relation(depth):
    omega = np.r_[0, np.geomspace(1e-8, 10, 40)]
    k = compy.wavenumber(omega, depth)
    np.testing.assert_allclose(9.79329 * k * np.tanh(k * depth), omega**2, rtol=1e-11)
    assert compy.wavenumber([0.1], depth)[0] > 0  # first bin need not be DC
    assert np.isclose(compy.wavenumber(0.1, depth), compy.wavenumber([0.1], depth)[0])
    np.testing.assert_equal(k, ffplot.wavenumber(omega, depth))
    np.testing.assert_equal(k, inv_compy.gravd(omega, depth))


@pytest.mark.parametrize("omega,depth", [([1], 0), ([-1], 4), ([np.nan], 4)])
def test_invalid_wavenumber_arguments(omega, depth):
    with pytest.raises(ValueError):
        compy.wavenumber(omega, depth)


def test_dtanh_zero_and_negative_extremes():
    x = np.array([-1e3, -1, 0, 1, 1e3])
    np.testing.assert_equal(inv_compy.dtanh(x), np.tanh(x))


def test_gravity_ratio_avoids_overflow():
    ratio, acceleration, _ = compy.gravitational_attraction(
        np.ones(3), 4000, [0, 0.01, 1]
    )
    assert np.all(np.isfinite(ratio))
    np.testing.assert_equal(ratio, acceleration)
    assert np.isclose(ratio[0], 2 * np.pi * 6.6743e-11 / 9.8)
    assert np.isclose(ratio[-1], np.pi * 6.6743e-11 / 9.8)


def test_amplitude_coherence_uncertainty():
    expected = 2 * np.sqrt(1 - 0.8**2) / (0.8 * np.sqrt(20))
    assert np.isclose(compy.Comliance_uncertainty(2, 0.8, 10), expected)
    assert np.isclose(ffplot.compliance_uncertainty(2, 0.8, 10), expected)
    assert compy.Comliance_uncertainty(2, 1, 10) == 0
    assert np.isinf(compy.Comliance_uncertainty(2, 0, 10))
    with pytest.raises(ValueError):
        compy.Comliance_uncertainty(2, 1.01, 10)


def test_calibration_gain_matches_known_scale_without_quantization():
    d = np.array([1.0, 2.0, 3.0])
    assert np.isclose(calibration.grid_search(d, d * 0.6637), 0.6637)
    assert calibration.grid_search(d, d * 7) == 7
    with pytest.raises(ValueError):
        calibration.grid_search(np.zeros(3), d)


def test_mackenzie_reference_values():
    assert calibration.calculate_speed_of_sound_in_water(0, 35, 0) == 1448.96
    # Direct substitution into Mackenzie's published nine-term polynomial.
    t, s, d = 10, 37, 1000
    expected = (
        1448.96
        + 4.591 * t
        - 0.05304 * t * t
        + 0.0002374 * t**3
        + 1.34 * (s - 35)
        + 0.0163 * d
        + 1.675e-7 * d * d
        - 0.01025 * t * (s - 35)
        - 7.139e-13 * t * d**3
    )
    assert np.isclose(calibration.calculate_speed_of_sound_in_water(t, s, d), expected)
