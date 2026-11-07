"""ComPy 3 audit: anomaly detectors, time-frequency helpers, references and figures."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from scipy import signal  # noqa: E402

from compy import anomaly as an  # noqa: E402
from compy import noise_models, results, timefreq  # noqa: E402
from compy.plotting import diagnostics as dg  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
YV = ROOT / "docs/results/yv_compliance_2026-10-28/results"
HAVE_ARCHIVE = (YV / "RR52" / "RR52_compliance_compy.npz").exists()
T0 = 1.35e9


# --------------------------------------------------------------------------
# Robust statistics and detectors
# --------------------------------------------------------------------------
def test_robust_z_ignores_outliers():
    x = np.r_[np.random.default_rng(0).normal(0, 1, 1000), [50.0, -40.0]]
    z = an.robust_z(x)
    assert abs(an.mad_sigma(x) - 1) < 0.1
    assert z[-2] > 40 and z[-1] < -35


def test_changepoint_found_at_the_step_and_not_in_red_noise():
    rng = np.random.default_rng(1)
    t = T0 + np.arange(3000) * 2100.0
    ar = signal.lfilter([1], [1, -0.6], rng.normal(0, 0.03, 3000))
    x = ar + np.where(np.arange(3000) >= 2000, np.log10(1.07), 0.0)
    cps = an.changepoints(t, x, min_size=100)
    assert len(cps) == 1
    assert abs(cps[0].time - t[2000]) < 30 * 2100
    assert abs(an._pct_from_log10(cps[0].shift) - 7.0) < 1.5
    assert an.changepoints(t, ar, min_size=100) == []


def test_stack_shift_weights_windows_by_pressure_power():
    t = T0 + np.arange(400) * 2100.0
    r = np.where(np.arange(400) >= 300, np.log10(1.1), 0.0)
    p = np.zeros(400)
    cps = an.changepoints(
        t, r + 1e-4 * np.random.default_rng(2).standard_normal(400), min_size=50
    )
    assert abs(an.stack_shift(t, r, p, cps) - 2.5) < 0.1  # 25 % of windows, +10 %
    p[300:] = 10.0  # ten times the power after the step
    # (300 x 1 x 1 + 100 x 10 x 1.1) / (300 + 1000) = 1.0769
    assert abs(an.stack_shift(t, r, p, cps) - 7.69) < 0.1


def test_periodicity_detects_m2_and_not_white_noise():
    rng = np.random.default_rng(3)
    t = T0 + np.sort(rng.uniform(0, 300 * 86400, 4000))
    noise = rng.normal(0, 0.05, t.size)
    m2 = 0.02 * np.cos(2 * np.pi * t / 3600 / an.TIDAL_PERIODS_H["M2"])
    hit = an.periodicity(t, noise + m2)["M2"]
    assert hit["p_value"] <= 0.01 and hit["p_tail"] < 1e-6
    assert abs(hit["amplitude"] - 0.02) < 0.005
    assert an.periodicity(t, noise)["M2"]["p_tail"] > 1e-3


def test_rayleigh_periodicity_detects_tidal_clustering_of_events():
    rng = np.random.default_rng(4)
    grid = T0 + np.arange(0, 300 * 86400, 2100.0)
    keep = rng.random(grid.size) < 0.5 * (
        1 + 0.6 * np.cos(2 * np.pi * grid / 3600 / 12.4206012)
    )
    out = an.rayleigh_periodicity(grid[keep])
    assert out["M2"]["p_tail"] < 1e-10 and abs(200 * out["M2"]["R"] - 60) < 10
    flat = an.rayleigh_periodicity(grid[rng.random(grid.size) < 0.5])
    assert flat["M2"]["p_tail"] > 1e-3


def test_level_dependence_slope():
    rng = np.random.default_rng(5)
    level = rng.uniform(20, 50, 3000)
    value = 0.002 * level + rng.normal(0, 0.01, 3000)
    dep = an.level_dependence(level, value)
    assert abs(dep["slope"] - 0.002) < 3e-4 and dep["spearman"] > 0.8


def test_kmeans_separates_two_regimes():
    rng = np.random.default_rng(6)
    x = np.r_[rng.normal(0, 1, (200, 2)), rng.normal(8, 1, (100, 2))]
    labels, _, _ = an.kmeans(x, k=2)
    assert {np.sum(labels == j) for j in (0, 1)} == {200, 100}


def test_common_mode_correlation():
    rng = np.random.default_rng(7)
    keys = [f"2013-{m:02d}" for m in range(1, 13)]
    shared = rng.normal(0, 1, 12)
    table = {s: dict(zip(keys, shared + 0.3 * rng.normal(0, 1, 12))) for s in "ABCD"}
    assert an.common_mode(table)["network_mean_corr"] > 0.7


@pytest.mark.parametrize("model", ["seismometer", "dpg"])
def test_phase_models_identify_a_response_mismatch(model):
    f = np.linspace(0.005, 0.015, 9)
    if model == "seismometer":
        ratio = an._hp2(f, 1 / 200.0, 0.707) / an._hp2(f, 1 / 240.0, 0.707)
        true_bias = np.abs(ratio)
        cal = np.abs(an._hp2(np.geomspace(0.02, 0.06, 9), 1 / 200.0, 0.707)
                     / an._hp2(np.geomspace(0.02, 0.06, 9), 1 / 240.0, 0.707))  # fmt: skip
    else:
        ratio = an._hp1(f, 0.0005) / an._hp1(
            f, 0.002
        )  # P residual; Z/P phase is minus its phase
        true_bias = 1 / np.abs(ratio)
        cal = 1 / np.abs(an._hp1(np.geomspace(0.02, 0.06, 9), 0.0005)
                         / an._hp1(np.geomspace(0.02, 0.06, 9), 0.002))  # fmt: skip
        ratio = np.conj(ratio)
    phase = 180 + np.degrees(np.angle(ratio))
    m = an.phase_models(f, phase)
    assert m["polarity"] == 180.0
    key = "seismometer_corner" if model == "seismometer" else "dpg_corner"
    assert m[key]["rms_deg"] < 0.3 < m["timing"]["rms_deg"]
    bias = 100 * (true_bias / cal.mean() - 1)
    assert np.all(m[key]["bias_low_pct"] <= bias + 0.2)
    assert np.all(m[key]["bias_high_pct"] >= bias - 0.2)


def test_score_and_ranking():
    a = an.Anomaly("A", "temporal", "x", "big", 9.0, 12.0)
    b = an.Anomaly("B", "noise", "y", "small", 3.0, 1.0)
    c = an.Anomaly("C", "noise", "z", "weak", 1.0, 50.0)
    ranked = an.rank([b, a, c])
    assert [x.station for x in ranked] == ["A", "B"]
    assert "| 1 | A |" in an.to_markdown(ranked)


# --------------------------------------------------------------------------
# Time-frequency helpers
# --------------------------------------------------------------------------
def _record(n=40000, fs=0.5, azimuth=70.0, seed=8):
    rng = np.random.default_rng(seed)
    b, a = signal.butter(2, [0.004, 0.05], btype="band", fs=fs)
    p = signal.lfilter(b, a, rng.standard_normal(n))
    tilt = signal.lfilter(b, a, rng.standard_normal(n))
    other = 0.5 * signal.lfilter(b, a, rng.standard_normal(n))
    az = np.radians(azimuth)
    h1 = np.cos(az) * tilt - np.sin(az) * other
    h2 = np.sin(az) * tilt + np.cos(az) * other
    z = p + 0.4 * tilt + 0.1 * rng.standard_normal(n)
    return {"Z": z, "P": p, "1": h1, "2": h2}, fs


def test_window_spectra_shapes_and_coherence():
    ch, fs = _record()
    spec = timefreq.window_spectra(ch, fs, 7200, 512, t0=T0)
    n_win = len(spec["times"])
    assert spec["ZZ"].shape == (n_win, 257) and spec["ZP"].dtype.kind == "c"
    coh = timefreq.band_average(
        spec["f"], timefreq.coherence(spec, "P", "Z"), (0.01, 0.03)
    )
    assert np.all(coh > 0.8)
    assert np.allclose(
        timefreq.coherence(spec, "Z", "P"), timefreq.coherence(spec, "P", "Z")
    )


def test_window_spectra_marks_windows_with_gaps():
    ch, fs = _record(n=20000)
    ch["Z"][5000] = np.nan
    spec = timefreq.window_spectra(ch, fs, 7200, 512)
    bad = np.isnan(spec["ZZ"][:, 10])
    assert bad.sum() == 1 and bad[1]


def test_tilt_direction_recovers_azimuth_and_admittance():
    ch, fs = _record(n=80000, azimuth=250.0)
    spec = timefreq.window_spectra(ch, fs, 14400, 512)
    tilt = timefreq.tilt_direction(spec, (0.005, 0.035))
    err = (tilt["azimuth"] - 250.0 + 180) % 360 - 180
    assert np.median(np.abs(err)) < 5
    assert abs(np.median(tilt["admittance"]) - 0.4) < 0.1


# --------------------------------------------------------------------------
# Reference models and archived products
# --------------------------------------------------------------------------
def test_peterson_models_spot_values():
    assert noise_models.nlnm(1 / 100.0) == pytest.approx(-216.47 + 15.70 * 2, abs=1e-9)
    assert noise_models.nhnm(1.0) == pytest.approx(-116.85, abs=1e-9)
    assert noise_models.nlnm(1 / 45.0) == pytest.approx(-187.5, abs=1e-9)
    assert np.isnan(noise_models.nlnm(20.0))
    f = np.geomspace(1e-4, 9.0, 400)
    assert np.all(noise_models.nhnm(f) > noise_models.nlnm(f))


def test_webb_bounds_load():
    webb = noise_models.webb_pressure_bounds()
    assert webb is not None
    f, lo, hi = webb
    assert np.all(np.diff(f) > 0) and np.all(hi > lo)


@pytest.mark.skipif(not HAVE_ARCHIVE, reason="archived results not present")
def test_load_station_and_scan_rr52():
    s = results.load_station(YV / "RR52", label="compy")
    assert s.depth == 2918.0 and s.in_band.sum() == 9
    assert s.window_levels.shape[1] == 4 and 0 < s.acceptance < 1
    keys, counts, months = s.monthly_matrix(30)
    assert months.shape == (len(keys), 9)
    found, info = an.scan_station(s)
    kinds = {a.kind for a in found}
    assert "phase_residual" in kinds and "changepoint" in kinds
    step = [c for c in info["levels"]["changepoints"]["transfer function"]]
    assert any(abs(c.time - 1.3795e9) < 20 * 86400 for c in step)  # September 2013


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------
def test_time_frequency_figures_render():
    ch, fs = _record(n=30000)
    spec = timefreq.window_spectra(ch, fs, 7200, 512, t0=T0)
    tilt = timefreq.tilt_direction(spec, (0.005, 0.035))
    figs = [
        dg.plot_spectrogram(spec["times"], spec["f"], {"Z": spec["ZZ"], "P": spec["PP"]}, band=(0.006, 0.02)),
        dg.plot_coherogram(spec["times"], spec["f"], {"Z–P": timefreq.coherence(spec, "P", "Z")}, threshold=0.9),
        dg.plot_tilt(spec["times"], tilt["azimuth"], tilt["admittance"], strength=tilt["coherence"]),
    ]  # fmt: skip
    for fig in figs:
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


@pytest.mark.skipif(not HAVE_ARCHIVE, reason="archived results not present")
def test_archive_figures_render():
    yv = results.discover(YV, label="compy")
    s = yv["RR52"]
    found, info = an.scan_station(s)
    models = {
        c: an.phase_models(v.freqs[v.in_band], v.phase_deg[v.in_band])
        for c, v in yv.items()
    }
    monthly = {c: an.scan_monthly(v)[1] for c, v in yv.items()}
    cm = an.common_mode(
        {c: dict(zip(m["keys"], m["band_mean_pct"])) for c, m in monthly.items() if m}
    )
    figs = [
        dg.plot_compliance_confidence(s),
        dg.plot_compliance_stability(s, steps=[1.3795e9]),
        dg.plot_window_levels(s, info["levels"]),
        dg.plot_noise_models(yv, highlight=("RR52",)),
        dg.plot_noise_models(yv, kind="pressure"),
        dg.plot_multistation(yv),
        dg.plot_phase_diagnostic(models),
        dg.plot_periodicity(
            {c: an.rayleigh_periodicity(v.window_starts) for c, v in yv.items()}
        ),
        dg.plot_acceptance_phase(yv),
        dg.plot_common_mode(cm),
        dg.plot_anomaly_ranking(an.rank(found)),
    ]
    for fig in figs:
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
