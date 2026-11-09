"""Tests of compy.response_check (pressure-gauge response from phase and tide)."""

import numpy as np
import pytest

from compy import response_check as rc

TWO_PI = 2 * np.pi
#: SIO-style DPG: one zero at DC, leak pole at 2 mHz (rad/s).
SIO = [(np.array([0j]), np.array([-TWO_PI * 0.002 + 0j]), True)]
#: INSU-style DPG: two zeros at DC, poles at 2.9 and 7.7 mHz.
INSU = [(np.array([0j]), np.array([-TWO_PI * 0.0029 + 0j]), True),
        (np.array([0j, -11656 + 0j]), np.array([-TWO_PI * 0.0077 + 0j, -57.471 + 0j]), True)]  # fmt: skip
F = np.linspace(0.0073, 0.017, 9)


def test_unchanged_response_has_no_residual():
    assert np.allclose(rc.predicted_residual(SIO, F, 1.0), 0)
    assert np.allclose(rc.compliance_bias(SIO, F, 1.0), 0)
    assert rc.predicted_ratio(SIO, 1.0) == pytest.approx(1.0)


def test_lower_true_corner_gives_positive_residual_and_larger_tide():
    res = rc.predicted_residual(SIO, F, 0.5)
    assert np.all(res > 0) and res[0] > res[-1]  # decays with frequency
    # far below the corner a first-order high-pass scales as f/fc: ratio 1/k
    assert rc.predicted_ratio(SIO, 0.5) == pytest.approx(2.0, rel=1e-3)
    # compliance ~ |H_meta/H_true| < 1 at low f, relative to the calibration band
    bias = rc.compliance_bias(SIO, F, 0.5)
    assert bias[0] < bias[-1] < 0.5


@pytest.mark.parametrize("stages,k", [(SIO, 0.45), (SIO, 1.8), (INSU, 0.6)])
def test_fit_recovers_the_factor(stages, k):
    res = rc.predicted_residual(stages, F, k)
    fit = rc.fit_leak_factor(
        stages, F, res + np.random.default_rng(1).normal(0, 0.1, len(F))
    )
    assert fit["k"] == pytest.approx(k, rel=0.05)
    assert fit["k_low"] <= k <= fit["k_high"]
    assert not fit["at_grid_edge"]


def test_mode_all_scales_every_low_pole():
    lowest = rc.predicted_ratio(INSU, 0.5, mode="lowest")
    both = rc.predicted_ratio(INSU, 0.5, mode="all")
    assert lowest == pytest.approx(2.0, rel=1e-2)
    assert both == pytest.approx(4.0, rel=1e-2)
    with pytest.raises(ValueError):
        rc.analog_response(SIO, F, mode="some")


def _hourly_tide(days=120, m2=50.0, k1=20.0, noise=2.0, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(0, days * 24) * 3600.0 + 1.6e9
    h = (t - t[0]) / 3600
    x = (m2 * np.cos(np.radians(rc.TIDAL_SPEED_DEG_H["M2"] * h) + 0.3)
         + k1 * np.cos(np.radians(rc.TIDAL_SPEED_DEG_H["K1"] * h) - 1.0)
         + 0.02 * h + rng.normal(0, noise, len(h)))  # fmt: skip
    return t, x


def test_tidal_harmonics_recovers_amplitudes_with_drift():
    t, x = _hourly_tide()
    out = rc.tidal_harmonics(t, x)
    assert out["amplitude"]["M2"] == pytest.approx(50.0, rel=0.02)
    assert out["amplitude"]["K1"] == pytest.approx(20.0, rel=0.05)
    assert out["noise"] < 1.0


def test_tidal_harmonics_survives_gaps_and_nan():
    # Regression: a NaN running median (next to a gap) used to make np.median
    # NaN and exclude every hour (pressure atlas tide check, s03 first run).
    t, x = _hourly_tide()
    keep = np.ones(len(t), bool)
    keep[500:530] = False  # a gap
    keep[2000:2003] = False  # a short isolated stretch follows
    x = x.copy()
    x[100] = np.nan
    out = rc.tidal_harmonics(t[keep], x[keep])
    assert out["amplitude"]["M2"] == pytest.approx(50.0, rel=0.03)
    assert out["n_used"] > 0.9 * keep.sum()


def test_paz_stages_reads_poles_and_zeros_stages():
    class PolesZerosResponseStage:
        zeros, poles = [0j], [-0.0126 + 0j]
        pz_transfer_function_type = "LAPLACE (RADIANS/SECOND)"

    class CoefficientsTypeResponseStage:
        pass

    class Response:
        response_stages = [PolesZerosResponseStage(), CoefficientsTypeResponseStage()]

    class Channel:
        response = Response()

    stages = rc.paz_stages(Channel())
    assert len(stages) == 1 and stages[0][2] is True
    assert stages[0][1][0] == pytest.approx(-0.0126)


@pytest.mark.parametrize("phase,sign", [([5.0, 3.0, 1.0], -1), ([-176.0, 178.0, 179.0], 1),
                                        ([60.0, 50.0, 40.0], None), ([np.nan, 1.0, 2.0], None)])  # fmt: skip
def test_pressure_sign_from_phase(phase, sign):
    assert rc.pressure_sign_from_phase(phase) == sign


def test_correct_compliance_inverts_the_bias():
    eta = np.full(len(F), 2e-11)
    k = 0.4
    corrected = rc.correct_compliance(eta, SIO, F, k)
    bias = rc.compliance_bias(SIO, F, k)
    assert np.allclose(corrected * (1 + bias / 100), eta)
    assert np.all(corrected > eta)  # k < 1 biases compliance low


def test_scale_leak_poles_keeps_unit_gain_at_normalisation():
    class PolesZerosResponseStage:
        def __init__(self):
            self.zeros = [0j]
            self.poles = [-TWO_PI * 0.002 + 0j, -TWO_PI * 5.0 + 0j]
            self.pz_transfer_function_type = "LAPLACE (RADIANS/SECOND)"
            self.normalization_frequency = 1.0
            self.normalization_factor = 1.0

    class Response:
        response_stages = [PolesZerosResponseStage()]

    class Channel:
        response = Response()

    ch = Channel()
    assert rc.scale_leak_poles(ch, 0.5) == 1
    s = ch.response.response_stages[0]
    assert s.poles[0] == pytest.approx(-TWO_PI * 0.001)
    assert s.poles[1] == pytest.approx(-TWO_PI * 5.0)  # untouched
    sn = 1j * TWO_PI * 1.0
    h = (
        s.normalization_factor
        * np.prod([sn - z for z in s.zeros])
        / np.prod([sn - p for p in s.poles])
    )
    assert abs(h) == pytest.approx(1.0)
