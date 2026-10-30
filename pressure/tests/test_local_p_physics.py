import numpy as np
import pytest

from compy_pressure import local_p_physics as lp
from compy_pressure.local_p_reference import simulate_plane_wave
from compy_pressure.rayleigh_calibration import water_ratio


def test_slowness_conversions():
    # 1 s/deg is 1/111.195 s/km; a TauP ray parameter in s/rad divides by R.
    assert lp.slowness_from_s_per_degree(1.0) == pytest.approx(1 / 111194.93, rel=1e-6)
    p_rad = 10.0 * 180 / np.pi  # 10 s/deg expressed per radian
    assert lp.slowness_from_ray_parameter(p_rad) == pytest.approx(
        lp.slowness_from_s_per_degree(10.0), rel=1e-12
    )
    with pytest.raises(ValueError):
        lp.slowness_from_ray_parameter(-1)


def test_vertical_cosine_and_p_bound():
    assert lp.vertical_cosine(0.0) == 1.0
    assert lp.vertical_cosine(1 / 3000.0, 1500.0) == pytest.approx(np.sqrt(0.75))
    with pytest.raises(ValueError, match="evanescent"):
        lp.vertical_cosine(1 / 1400.0, 1500.0)
    # Any P wave with v_P >= 5 km/s stays within 17.5 degrees of vertical in water.
    assert lp.eta_min_for_p() == pytest.approx(np.sqrt(1 - 0.09))
    assert np.degrees(np.arccos(lp.eta_min_for_p())) < 17.5


def test_incompressible_limit_and_series():
    f = np.array([0.0, 1e-4, 0.01])
    ratio = lp.steady_state_p_over_a(f, 4000.0, density=1000.0)
    assert ratio[0] == pytest.approx(1000.0 * 4000.0)  # p = rho h a (An 2017 Eq. 6)
    theta = lp.water_phase(f, 4000.0)
    assert ratio / 4e6 == pytest.approx(1 + theta**2 / 3 + 2 * theta**4 / 15, rel=1e-5)


def test_resonance_guard_and_poles():
    depth = 152.0
    f0 = lp.resonance_frequencies(depth)[0]
    assert f0 == pytest.approx(1500.0 / (4 * depth))
    assert lp.resonance_frequencies(depth, q=1e-4, count=2)[1] == pytest.approx(
        3 * 1500.0 / (4 * depth * lp.vertical_cosine(1e-4))
    )
    lp.steady_state_p_over_a([0.89 * f0], depth)
    with pytest.raises(ValueError, match="resonance"):
        lp.steady_state_p_over_a([0.91 * f0], depth)
    with pytest.raises(ValueError):
        lp.steady_state_p_over_a([1.0], -5.0)


def test_consistency_with_production_water_ratio():
    """Shared physics with the Rayleigh code (a consistency check, not independent)."""
    f = np.linspace(0.01, 1.5, 20)
    q = 1 / 8000.0
    ours = lp.steady_state_p_over_a(f, 152.0, q, 1028.0, 1500.0)
    theirs = water_ratio(f, np.full_like(f, 8000.0), 152.0, 1028.0, 1500.0)
    assert ours == pytest.approx(theirs, rel=1e-10)


def test_regime_classification_uses_bounds():
    eta = (lp.eta_min_for_p(), 1.0)
    c = lp.SOUND_SPEED_RANGE
    regime, tau_min, tau_max = lp.classify_regime(3.0, (2900, 3100), eta, c)
    assert regime == "pre_multiple" and tau_min > 3.0
    assert lp.classify_regime(3.0, (140, 165), eta, c)[0] == "steady_state"
    assert lp.classify_regime(1.0, (140, 165), eta, c)[0] == "mixed"
    assert tau_max >= tau_min


def test_model_variation_guard():
    # At 152 m, F changes by ~3.5% over 0.5+-0.2 Hz but ~50% over 2.0+-0.2 Hz.
    variation = lp.model_variation(np.array([0.5, 2.0]), 0.2, 152.0)
    assert variation[0] < 0.05 < variation[1]
    assert np.isinf(lp.model_variation([2.3], 0.4, 152.0)[0])


def test_polarity_and_circular_mean():
    assert list(lp.polarity_class([5, -170, 179, 90, -30])) == [1, -1, -1, 0, 1]
    assert lp.circular_mean_deg([170, -170]) == pytest.approx(180, abs=1e-9) or (
        lp.circular_mean_deg([170, -170]) == pytest.approx(-180, abs=1e-9)
    )


def test_coherence_zero_level():
    assert lp.coherence_zero_level(4) == pytest.approx(1 - 0.05 ** (1 / 3))
    with pytest.raises(ValueError):
        lp.coherence_zero_level(1)


def test_multitaper_psd_scaling_and_transfer():
    rng = np.random.default_rng(1)
    fs, n = 50.0, 4000
    v = rng.standard_normal(n)
    p = -2.5 * v  # known real transfer with negative polarity
    f, pk, vk = lp.multitaper_spectra(p, v, fs, nw=3.0)
    out = lp.cross_spectral_products(pk, vk)
    band = (f > 1) & (f < 20)
    # White noise of unit variance has one-sided PSD 2/fs.
    assert np.median(out["svv"][band]) == pytest.approx(2 / fs, rel=0.1)
    assert np.allclose(out["h1"][band], -2.5)
    assert np.allclose(out["coherence"][band], 1.0)
    with pytest.raises(ValueError):
        lp.multitaper_spectra(p, v[:-1], fs)


def test_aic_onset_finds_step_in_variance():
    rng = np.random.default_rng(2)
    x = np.concatenate([0.1 * rng.standard_normal(500), rng.standard_normal(300)])
    assert abs(lp.aic_onset(x) - 500) <= 3


@pytest.mark.parametrize("q", [0.0, 0.15e-3])
@pytest.mark.parametrize("solid", [(2500.0, 2000.0), (4500.0, 2600.0)])
def test_independent_simulation_reproduces_steady_state(q, solid):
    """Full transient: P/A equals Deng Eq. 5 and does not depend on the seafloor."""
    depth = 150.0
    t, p, v, _ = simulate_plane_wave(
        depth=depth,
        q=q,
        solid_speed=solid[0],
        solid_density=solid[1],
        dz=1.0,
        duration=30.0,
        peak_frequency=2.0,
    )
    f = np.fft.rfftfreq(len(t), t[1] - t[0])
    ratio = np.fft.rfft(p) / (1j * 2 * np.pi * f * np.fft.rfft(v) + 1e-300)
    band = (f > 0.2) & (f < 0.85 * lp.resonance_frequencies(depth, q)[0])
    model = lp.steady_state_p_over_a(f[band], depth, q, 1028.0, 1500.0)
    assert ratio[band].real / model == pytest.approx(1.0, abs=0.025)
    assert np.max(np.abs(ratio[band].imag / ratio[band].real)) < 0.01


@pytest.mark.parametrize("q", [0.0, 0.18e-3])
def test_independent_simulation_first_transmission_only_before_multiple(q):
    depth = 3000.0
    t, p, v, onset = simulate_plane_wave(
        depth=depth,
        q=q,
        solid_speed=3000.0,
        solid_density=2200.0,
        dz=4.0,
        duration=12.0,
        peak_frequency=6.0,
        source_below_seafloor=800.0,
    )
    tau = lp.round_trip_time(depth, q)
    model = lp.first_transmission_p_over_v(q)
    gated = (t > onset - 0.3) & (t < onset + 0.9 * tau)
    fit = np.dot(p[gated], v[gated]) / np.dot(v[gated], v[gated])
    assert fit / model == pytest.approx(1.0, abs=0.01)
    # The same relation over the full window is badly wrong (SOC-012 correction).
    full = np.dot(p, v) / np.dot(v, v)
    assert abs(full / model - 1) > 0.5
    # Multitaper transfer inside the gate is flat and in phase.
    fs = 1 / (t[1] - t[0])
    f, pk, vk = lp.multitaper_spectra(p[gated], v[gated], fs, nw=2.5)
    out = lp.cross_spectral_products(pk, vk)
    band = (f > 3) & (f < 10)
    assert np.abs(out["h1"][band]) / model == pytest.approx(1.0, abs=0.03)
    assert np.all(lp.polarity_class(np.angle(out["h1"][band], deg=True)) == 1)
