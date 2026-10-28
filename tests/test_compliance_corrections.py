"""Regression tests for the 2.2 compliance corrections (synthetic, offline)."""

import numpy as np
import pytest
from obspy import Stream, Trace, UTCDateTime

import compy
from compy_numerics import (
    FREE_AIR_GRADIENT,
    GRAVITATIONAL_CONSTANT,
    GRAVITY,
    HANN_HALF_OVERLAP_DOF,
    compliance_uncertainty,
    welch_averages,
)
from compy_processing import (
    GJI2025,
    SelectionCriteria,
    _covered_samples,
    coherent_band,
    estimate,
    gravity_corrected_compliance,
    stacked_compliance,
)
from compy_streams import split_stream

DEPTH = 4000.0


def make_stream(channels, n, fs=1.0):
    return Stream(
        [
            Trace(
                np.zeros(n),
                header={
                    "network": "XX",
                    "station": "TEST",
                    "channel": channel,
                    "sampling_rate": fs,
                    "starttime": UTCDateTime(0),
                },
            )
            for channel in channels
        ]
    )


def coupled_stream(n, coupling, noise, seed=0):
    """Z = coupling * P + white noise, P white, both in SI units."""
    rng = np.random.default_rng(seed)
    st = make_stream(("BHZ", "BDH"), n)
    p = rng.normal(0, 100, n)
    st.select(channel="BDH")[0].data = p
    st.select(channel="BHZ")[0].data = coupling * p + rng.normal(0, noise, n)
    return st


def expected(f, coupling, pressure_amp):
    return gravity_corrected_compliance(f, DEPTH, abs(coupling), pressure_amp)


def test_gravity_ratio_uses_shared_constants():
    ratio, _, _ = compy.gravitational_attraction(np.ones(1), DEPTH, [0.0])
    assert ratio[0] == pytest.approx(
        2 * np.pi * GRAVITATIONAL_CONSTANT / GRAVITY, rel=1e-12, abs=0
    )


def test_gravity_correction_formula_and_sign():
    f = np.array([0.008, 0.012])
    k = compy.wavenumber(2 * np.pi * f, DEPTH)
    ratio, _, _ = compy.gravitational_attraction(np.ones(2), DEPTH, f)
    omega2 = (2 * np.pi * f) ** 2
    eta = gravity_corrected_compliance(f, DEPTH, 1e-8, np.ones(2))
    np.testing.assert_allclose(
        eta, k * (omega2 * 1e-8 + ratio) / (omega2 + FREE_AIR_GRADIENT), rtol=1e-14
    )
    # Wave attraction opposes loading, so the correction increases |compliance|.
    assert np.all(eta > k * 1e-8 * omega2 / (omega2 + FREE_AIR_GRADIENT))


def test_legacy_estimator_does_not_scale_gravity_by_coherence():
    from compy_processing import spectra

    coupling = 6.6e-9  # omega² * coupling ~ gravity term near 0.01 Hz
    st = coupled_stream(3600, coupling, 40 * coupling, seed=1)
    curves, _, windows, fc, f, *_ = estimate(
        st, DEPTH, 1.0, 1, 0.005, 0.02, True, nseg=512
    )
    _, zz, pp, _, coh = spectra(windows[0], 512)
    keep = (f >= 0.001) & (f <= 0.1)
    k = compy.wavenumber(2 * np.pi * f, DEPTH)
    grav, _, _ = compy.gravitational_attraction(np.sqrt(pp), DEPTH, f)
    omega2 = (2 * np.pi * f) ** 2
    ratio = np.sqrt(zz / pp)
    fixed = k * (omega2 * coh * ratio + grav) / (omega2 + FREE_AIR_GRADIENT)
    old = k * coh * (omega2 * ratio + grav) / (omega2 + FREE_AIR_GRADIENT)
    np.testing.assert_allclose(curves[0], fixed[keep], rtol=1e-12)
    low = (fc >= 0.005) & (fc <= 0.012)
    assert np.all(np.abs(curves[0][low] / old[keep][low] - 1) > 1e-3)


def test_stacked_estimator_is_unbiased_at_low_coherence():
    # Noise lowers coherence to ~0.7 where the gravity term is ~half of the
    # signal, so scaling it by coherence (the 2.1 behaviour) would bias ~15 %.
    coupling = 6.6e-9
    st = coupled_stream(86400 * 4, coupling, 100 * coupling, seed=2)
    open_gates = SelectionCriteria(0.0, None, None, "open")
    result = stacked_compliance(
        split_stream(st, 7200), DEPTH, 1.0, (0.008, 0.02), open_gates, nseg=1024
    )
    assert 0.6 < np.median(result.coherence) < 0.85
    sel = (result.all_freqs >= 0.008) & (result.all_freqs <= 0.02)
    truth = expected(result.freqs, coupling, np.sqrt(result.stack["S_pp"][sel]))
    np.testing.assert_allclose(result.compliance, truth, rtol=0.06)
    omega2 = (2 * np.pi * result.freqs) ** 2
    k = compy.wavenumber(2 * np.pi * result.freqs, DEPTH)
    gravity_part = truth - k * omega2 * coupling / (omega2 + FREE_AIR_GRADIENT)
    old = truth - (1 - result.coherence) * gravity_part
    assert np.median(np.abs(old / truth - 1)) > 0.06


def test_stacked_uncertainty_follows_bendat_piersol():
    coupling = 1e-8
    st = coupled_stream(86400, coupling, 50 * coupling, seed=3)
    windows = split_stream(st, 2100)
    result = stacked_compliance(
        windows,
        DEPTH,
        1.0,
        (0.005, 0.02),
        SelectionCriteria(0.0, None, None, "open"),
        nseg=512,
        gravity_correction=False,
    )
    assert result.n_windows == len(windows)
    assert result.n_averages == pytest.approx(len(windows) * welch_averages(2100, 512))
    np.testing.assert_allclose(
        result.uncertainty,
        compliance_uncertainty(result.compliance, result.coherence, result.n_averages),
        rtol=1e-12,
    )
    # The scatter of per-window curves about the stack agrees with the
    # analytical error to within sampling variability.
    sem = result.window_std / np.sqrt(result.n_windows)
    assert 0.5 < np.median(sem / result.uncertainty) < 2


def test_stacked_phase_reports_negative_coupling():
    st = coupled_stream(86400, -1e-8, 1e-7, seed=4)
    result = stacked_compliance(
        split_stream(st, 3600),
        DEPTH,
        1.0,
        (0.005, 0.02),
        SelectionCriteria(0.0, None, None, "open"),
        nseg=512,
    )
    assert np.all(np.abs(np.abs(result.phase_deg) - 180) < 10)


def test_gain_scales_pressure_not_vertical():
    coupling = 1e-8
    st = coupled_stream(86400, coupling, 0.0, seed=5)
    gates = SelectionCriteria(0.0, None, None, "open")
    args = (split_stream(st, 3600), DEPTH)
    one = stacked_compliance(
        *args, 1.0, (0.005, 0.02), gates, 512, gravity_correction=False
    )
    half = stacked_compliance(
        *args, 0.5, (0.005, 0.02), gates, 512, gravity_correction=False
    )
    np.testing.assert_allclose(half.compliance, 2 * one.compliance, rtol=1e-12)


def test_selection_criteria_report_reasons():
    f = np.linspace(0.001, 0.05, 50)
    band = (f >= 0.006) & (f <= 0.02)
    coh = np.full_like(f, 0.95)
    pp = np.full_like(f, 1e3)  # 30 dB
    zz = 10 ** (-160 / 10) / (2 * np.pi * f) ** 4  # -160 dB acceleration
    assert GJI2025.check(f, zz, pp, coh, band) is None
    assert GJI2025.check(f, zz, pp, coh * 0.5, band) == "coherence"
    assert GJI2025.check(f, zz, pp * 1e3, coh, band) == "pressure_psd"
    assert GJI2025.check(f, zz * 1e3, pp, coh, band) == "acceleration_psd"


def test_failed_gates_raise_with_counts():
    st = coupled_stream(36000, 0.0, 1e-6, seed=6)
    with pytest.raises(ValueError, match="coherence"):
        stacked_compliance(split_stream(st, 3600), DEPTH, 1.0, (0.005, 0.02))


def test_coherent_band_picks_widest_run_in_search_range():
    f = np.round(np.arange(0, 0.05, 0.001), 6)
    coh = np.where((f >= 0.006) & (f <= 0.018), 0.9, 0.3)
    coh[f == 0.03] = 0.95
    assert coherent_band(f, coh, 0.8, (0.002, 0.04)) == pytest.approx((0.006, 0.018))
    assert coherent_band(f, coh * 0, 0.8, (0.002, 0.04)) is None


def test_welch_averages_and_overlap_coverage():
    assert welch_averages(1024, 1024) == pytest.approx(HANN_HALF_OVERLAP_DOF)
    assert welch_averages(2048, 1024) == pytest.approx(3 * HANN_HALF_OVERLAP_DOF)
    st = make_stream(("BHZ", "BDH"), 600)
    windows = compy.cut_stream_with_overlap(st, 100, 50)
    assert _covered_samples(windows, 100) == 600
    assert _covered_samples(windows[:1], 100) == 100


def test_legacy_theory_uses_all_windows_not_the_first():
    coupling = 1e-8
    st = coupled_stream(4 * 3600, coupling, 30 * coupling, seed=7)
    *_, fc, f, scatter, theory, _, _ = estimate(
        st, DEPTH, 1.0, 1, 0.005, 0.02, True, 512
    )
    curves = estimate(st, DEPTH, 1.0, 1, 0.005, 0.02, True, 512)[0]
    assert theory.shape == fc.shape
    # Overlapping 1-h windows over 4 h: the error reflects 4 h of data, which
    # is smaller than a single window's error.
    single = compliance_uncertainty(
        curves[0], np.full_like(fc, 0.5), welch_averages(3600, 512)
    )
    assert np.median(theory) < np.median(single)


def test_rotate_without_per_window_cleaning(monkeypatch):
    calls = []

    class Spy:
        def __init__(self, *args, **kwargs):
            calls.append(1)

    class Rotator:
        azimuth, angle = 10.0, 1.0

        def __init__(self, *args, **kwargs):
            pass

        def apply(self, stream, **kwargs):
            return stream.copy()

    monkeypatch.setattr(compy.tiskit, "DataCleaner", Spy)
    monkeypatch.setattr(compy.tiskit, "CleanRotator", Rotator)
    rng = np.random.default_rng(8)
    st = make_stream(("BHZ", "BH1", "BH2", "BDH"), 7200)
    for tr in st:
        tr.data = rng.normal(size=7200) * (0.01 if tr.stats.channel[-1] in "12" else 1)
    compy.Rotate(st, 1, plot=False, clean_horizontals=False)
    assert not calls


def test_remove_coherent_noise_reduces_tilt_leakage():
    rng = np.random.default_rng(9)
    n = 2 * 86400
    st = make_stream(("BHZ", "BH1", "BH2", "BDH"), n)
    signal_z = rng.normal(0, 1, n)
    h1 = np.cumsum(rng.normal(0, 1, n))
    h1 -= np.convolve(h1, np.ones(601) / 601, mode="same")
    st.select(channel="BH1")[0].data = h1
    st.select(channel="BH2")[0].data = rng.normal(0, 1, n)
    st.select(channel="BDH")[0].data = rng.normal(0, 1, n)
    st.select(channel="BHZ")[0].data = signal_z + 0.5 * h1
    cleaned = compy.remove_coherent_noise(st, 24, spectral_window_s=1024)
    z = cleaned.select(component="Z")[0].data
    residual = np.var(z - signal_z[: len(z)]) / np.var(0.5 * h1)
    assert residual < 0.05
    np.testing.assert_array_equal(
        cleaned.select(channel="BDH")[0].data,
        st.select(channel="BDH")[0].data[: len(z)],
    )


def test_remove_coherent_noise_accepts_spans_crossing_window_edges():
    import tiskitpy

    rng = np.random.default_rng(10)
    n = 2 * 86400
    st = make_stream(("BHZ", "BH1", "BH2", "BDH"), n)
    for tr in st:
        tr.data = rng.normal(size=n)
    t0 = st[0].stats.starttime
    spans = tiskitpy.TimeSpans(
        spans=[(t0 + 86400 - 1800, t0 + 86400 + 1800), (t0 + 3600, t0 + 7200)]
    )
    cleaned = compy.remove_coherent_noise(
        st, 24, avoid_spans=spans, spectral_window_s=1024, cover_all=True
    )
    assert len(cleaned.select(component="Z")[0]) == n
