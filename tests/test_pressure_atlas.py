"""Tests of compy_pressure.atlas: the native-rate estimator, the screens, the models and the whole build."""

import json

import numpy as np
import pandas as pd
import pytest
from obspy import Stream, Trace, UTCDateTime
from obspy.core.inventory.response import Response
from scipy import signal

from compy_pressure.atlas import classify, config, cube, models, native, screens

CFG = config.load()


def _flat_response():
    import warnings

    warnings.filterwarnings("ignore", message="ObsPy can not map unit")
    return Response.from_paz(
        [],
        [],
        1.0,
        input_units="PA",
        output_units="COUNTS",
        normalization_frequency=1.0,
    )


def _synthetic_channel(tmp_path, fs, sigma=10.0, days=2, tide=0.0, name="t"):
    n = int(days * 86400 * fs)
    rng = np.random.default_rng(0)
    x = rng.normal(0, sigma, n) + tide * np.sin(2 * np.pi * np.arange(n) / fs / 44714.0)
    t0 = UTCDateTime(2020, 1, 1)
    tr = Trace(
        x,
        header=dict(
            network="XX",
            station="T",
            location="00",
            channel="HDH",
            sampling_rate=fs,
            starttime=t0,
        ),
    )

    def reader(p):
        d = int(str(p).split(".")[-1]) - 1
        return Stream(
            [
                tr.slice(
                    t0 + d * 86400, t0 + (d + 1) * 86400 - 1 / fs, nearest_sample=False
                ).copy()
            ]
        )

    files = [f"XX.T.00.HDH.D.2020.{d + 1:03d}" for d in range(days)]
    out = tmp_path / f"{name}.npz"
    native.process_channel(files, _flat_response(), "XX.T.00.HDH", fs, out, read=reader)
    return native.load(out), 10 * np.log10(2 * sigma**2 / fs)


def test_grid_shares_the_atlas_nodes():
    g = native.GRID
    assert g[0] == pytest.approx(0.02098, rel=1e-3) and g[-1] == pytest.approx(
        1060 * 2**0.25
    )
    assert np.any(np.isclose(g, 1060.0)) and np.allclose(np.diff(np.log2(g)), 1 / 8)


def test_welch_matches_scipy():
    rng = np.random.default_rng(1)
    for fs in (1.0, 62.5):
        x = np.cumsum(rng.normal(size=int(10800 * fs)))
        f1, p1 = native.welch_psd(x, fs)
        n = int(2700 * fs)
        f2, p2 = signal.welch(
            x, fs, "hann", n, int(round(0.75 * n)), detrend="linear", scaling="density"
        )
        assert np.array_equal(f1, f2) and np.max(np.abs(p1 - p2) / p2) < 1e-9


def test_band_weights_are_power_means():
    f = np.linspace(0, 10, 1001)
    w, resolved = native.band_weights(f, [1.0, 0.5])
    assert resolved.all() and np.allclose(np.asarray(w.sum(axis=1)).ravel(), 1.0)
    assert (w @ np.full(len(f), 3.0)) == pytest.approx([3.0, 3.0])


@pytest.mark.parametrize("fs", [1.0, 50.0])
def test_white_noise_level_is_recovered(tmp_path, fs):
    z, expect = _synthetic_channel(tmp_path, fs)
    med = np.nanmedian(z["levels_db"], axis=0)
    ok = np.isfinite(med)
    assert ok.sum() > 50 and np.all(np.abs(med[ok] - expect) < 1.5)
    assert z["periods"][ok].max() > 1200
    assert z["periods"][ok].min() == pytest.approx(1 / (0.45 * fs), rel=0.2)


def test_tide_does_not_leak_into_long_periods(tmp_path):
    z, expect = _synthetic_channel(tmp_path, 10.0, tide=1e5)
    med = np.nanmedian(z["levels_db"], axis=0)
    long = (z["periods"] > 300) & np.isfinite(med)
    assert long.any() and np.all(np.abs(med[long] - expect) < 1.5)


def test_histogram_and_shift_keep_counts():
    edges = np.arange(-10, 11, 1.0)
    lv = np.array([[0.2, 5.5], [0.7, np.nan], [-5.0, -9.5]])
    h = cube.histogram(lv, np.array([0, 0, 5]), edges)
    assert h.sum() == 5 and h[0, 0, 10] == 2 and h[5, 1, 0] == 1
    s = models.shift_hist(h, 2.5)
    assert s.sum() == pytest.approx(5) and s[0, 0, 12] == pytest.approx(1.0)


def test_edge_days_cut_a_descent():
    days = np.repeat(np.arange(30), 8)
    level = np.zeros(len(days))
    level[days < 2] = 40.0
    assert cube.edge_days(days, level, 15.0, 20) == {0, 1}


def test_percentiles_of_a_known_distribution():
    edges = np.arange(0, 101, 1.0)
    P = np.ones((1, 100)) / 100
    pc, _ = models.percentiles(P, np.ones(1), edges, [5, 50, 95])
    assert pc[50][0] == pytest.approx(50, abs=0.01) and pc[5][0] == pytest.approx(
        5, abs=0.01
    )


def test_piecewise_fit_recovers_two_lines():
    T = np.logspace(0, 3, 40)
    k = T[20]  # knots sit on sample periods
    y = np.where(
        T < k, 10 + 20 * np.log10(T), 10 + 20 * np.log10(k) - 10 * np.log10(T / k)
    )
    f = models.fit_piecewise(T, y, 0.1, 6)
    assert f["segments"] == 2 and f["rms_db"] < 0.1
    assert f["table"][0]["B_db_per_decade"] == pytest.approx(20, abs=0.5)


def test_band_taper():
    u = models.band_taper(0.01, 1.0, periods=np.array([0.5, 1.0, 10.0, 100.0, 200.0]))
    assert (
        u[0] == 0 and u[2] == 1 and u[3] == pytest.approx(0.0, abs=1e-12) and u[4] == 0
    )


def test_classification_rules():
    assert classify.sensor_type("Cox-Webb differential pressure gauge", "BDH") == "DPG"
    assert classify.sensor_type("Paroscientific 8B7000", "HDH") == "APG"
    assert classify.sensor_type("", "HDH") == "Hydrophone"
    assert classify.basin(-30, 60) == "Indian" and classify.basin(70, 0) == "Arctic"
    assert (
        classify.local_season(1, 10) == "winter"
        and classify.local_season(7, -10) == "winter"
    )
    assert list(classify.season_months("winter", -5)) == [5, 6, 7]


def test_day_screen_drops_a_day_peers_do_not_see():
    days = np.arange(10)
    daily = {
        k: pd.DataFrame(
            {"sec_microseism_3_10s": 20.0 + 0.1 * np.arange(10)}, index=days
        )
        for k in "abcd"
    }
    daily["a"].loc[4, "sec_microseism_3_10s"] = 60.0
    t = pd.DataFrame(
        dict(
            key=list("abcd"),
            net="N",
            depth_bin="3-4 km",
            depth_m=3500.0,
            f_lo=0.001,
            f_hi=1.0,
        )
    )
    cfg = config.load(bands={"sec_microseism_3_10s": [3.0, 10.0]})
    cfg["bands"] = {"sec_microseism_3_10s": [3.0, 10.0]}
    drops, _ = screens.day_screen(daily, t, cfg)
    assert drops["a"] == {4} and drops["b"] == set()


def test_calibration_screen_excludes_a_network():
    t = pd.DataFrame(dict(key=[f"k{i}" for i in range(6)], net=["A"] * 3 + ["B"] * 3, type=["DPG"] * 3 + ["Hydrophone"] * 3,
                          depth_bin="3-4 km"))  # fmt: skip
    off = pd.Series([0.5, -1.0, 2.0, 25.0, 30.0, 18.0], index=t.key)
    cal = screens.calibration_screen(t, off, CFG)
    assert not cal.cal_excluded[:3].any() and cal.cal_excluded[3:].all()
    assert set(cal.cal_reason[3:]) == {"network fails P-Z level"}


def test_geography_distance_to_a_line():
    from compy_pressure.atlas import geography

    near = geography.NearestLine(
        geography._densify(np.array([[0.0, 0.0], [10.0, 0.0]]), 1.0)
    )
    d, _ = near.distance_km(1.0, 5.0)
    assert d[0] == pytest.approx(111.2, abs=0.6)


def test_build_end_to_end(tmp_path):
    """Twelve synthetic gauges in two networks: one with a +40 dB day, one 30 dB off its peers."""
    from compy_pressure.atlas import pipeline

    rows = []
    nd = tmp_path / "native"
    nd.mkdir()
    for i in range(12):
        sigma = 10.0 * (1 + 0.05 * i)
        if i == 11:
            sigma *= 10 ** (30 / 20)
        z, _ = _synthetic_channel(tmp_path, 2.0, sigma=sigma, days=8, name=f"g{i}")
        if i == 3:  # one day 40 dB louder on gauge 3
            day = np.floor(z["starts"] / 86400).astype(int)
            z["levels_db"][day == day[0] + 4] += 40
        p = nd / f"N{i // 6}.S{i}.00.HDH.npz"
        np.savez_compressed(
            p, **{k: v for k, v in z.items() if k != "meta"}, meta=json.dumps(z["meta"])
        )
        rows.append(
            dict(
                seed_id=f"N{i // 6}.S{i}.00.HDH",
                net=f"N{i // 6}",
                sta=f"S{i}",
                cha="HDH",
                lat=-10.0 + i,
                lon=150.0,
                depth_m=3500.0,
                fs=2.0,
                type="DPG",
                brand="SIO Cox-Webb DPG",
                product=str(p),
                gain_db=0.0,
                coast_km=100.0 + 50 * i,
                shelf_break_km=60.0 + 40 * i,
                gauge=f"g{i}" if i % 2 else None,
            )
        )  # fmt: skip  (half of the gauge ids are missing)
    pd.DataFrame(rows).to_csv(tmp_path / "channels.csv", index=False)
    cfg = config.load(
        population={"min_windows": 50, "level_gate_bands_s": [[2.5, 4.0]]},
        windows={"edge_days": 2},
    )
    A = pipeline.build(
        tmp_path / "channels.csv", tmp_path / "atlas", cfg, workers=2, make_figures=True
    )
    pop = pd.read_csv(tmp_path / "atlas" / "tables" / "population.csv")
    assert pop.set_index("sta").drop_reason["S11"].startswith("2.5-4 s level")
    assert pop.in_atlas.sum() == 11
    days = pd.read_csv(tmp_path / "atlas" / "tables" / "screen_days.csv")
    assert (days.key == "N0.S3.00.HDH").sum() == 1
    tab = pd.read_csv(tmp_path / "atlas" / "models" / "all" / "all.csv")
    ok = tab.gauges >= 8
    assert ok.sum() > 20 and tab.p50_db[ok].between(15, 25).all()
    assert (tmp_path / "atlas" / "figures" / "global_density.png").exists()
    assert len(A.index_table()) > 5


def test_edge_days_spare_a_loud_day_inside_a_short_record():
    days = np.repeat(np.arange(30), 8)
    level = np.zeros(len(days))
    level[days == 15] = 40.0
    assert cube.edge_days(days, level, 15.0, 20) == set()
    level[days == 0] = 40.0  # a real edge is still cut, and only that edge
    assert cube.edge_days(days, level, 15.0, 20) == {0}


def test_band_taper_without_lower_limit():
    u = models.band_taper(0.0, 1.0, periods=np.array([10.0, 1000.0]))
    assert np.all(u == 1)


def test_usable_limits_keep_the_longest_periods_when_no_limit_is_given():
    from compy_pressure.atlas import pipeline

    t = pd.DataFrame(dict(fs=[50.0, 50.0, 50.0]))
    for col in (None, [np.nan, 0.0, 0.0]):
        if col is not None:
            t["f_lo_hz"] = col
        f_lo, f_hi = pipeline.usable_limits(t, CFG)
        for lo, hi in zip(f_lo, f_hi):
            assert models.band_taper(lo, hi)[np.argmax(models.GRID)] == pytest.approx(
                1.0
            )


def test_piecewise_fit_keeps_the_best_fit_when_segments_run_out():
    rng = np.random.default_rng(1)
    T = np.logspace(0, 3, 30)
    y = rng.normal(0, 8.0, 30)
    few = models.fit_piecewise(T, y, 1.0, 15)
    many = models.fit_piecewise(T, y, 1.0, 20)
    assert many["segments"] > 1 and many["rms_db"] <= few["rms_db"] + 1e-9


def test_window_takes_the_response_of_its_epoch():
    first, second = np.arange(0.0, 11.0), np.arange(5.0, 16.0)
    rows = native.pick_epoch_windows([first, second], [(0.0, 8.0), (8.0, np.inf)])
    assert [s for s, _, _ in rows] == list(np.arange(0.0, 16.0))
    assert all(k == (0 if s < 8 else 1) for s, k, _ in rows)
    # a window before the first epoch is kept from the first epoch that computed it
    rows = native.pick_epoch_windows([np.array([-3.0, 1.0])], [(0.0, 8.0)])
    assert rows == [(-3.0, 0, 0), (1.0, 0, 1)]


def _two_trace_reader(fs, traces):
    t0 = UTCDateTime(2020, 1, 1)

    def reader(p):
        d = int(str(p).split(".")[-1]) - 1
        return Stream([tr.slice(t0 + d * 86400, t0 + (d + 1) * 86400 - 1 / tr.stats.sampling_rate, nearest_sample=False).copy()
                       for tr in traces])  # fmt: skip

    return reader


def _trace(fs, sigma, days, loc="00", seed=0):
    x = np.random.default_rng(seed).normal(0, sigma, int(days * 86400 * fs))
    return Trace(x, header=dict(network="XX", station="T", location=loc, channel="HDH", sampling_rate=fs,
                                starttime=UTCDateTime(2020, 1, 1)))  # fmt: skip


def test_traces_at_another_sampling_rate_are_left_out(tmp_path):
    fs = 2.0
    files = [f"XX.T.00.HDH.D.2020.{d + 1:03d}" for d in range(2)]
    reader = _two_trace_reader(fs, [_trace(2.5, 10.0, 2)])
    meta = native.process_channel(
        files, _flat_response(), "XX.T.00.HDH", fs, tmp_path / "o.npz", read=reader
    )
    assert (
        meta["traces_other_rate"] > 0
        and len(native.load(tmp_path / "o.npz")["starts"]) == 0
    )


def test_two_location_codes_are_not_mixed(tmp_path):
    fs = 2.0
    files = [f"XX.T.00.HDH.D.2020.{d + 1:03d}" for d in range(2)]
    reader = _two_trace_reader(
        fs, [_trace(fs, 10.0, 2, "00", 0), _trace(fs, 1000.0, 2, "10", 1)]
    )
    native.process_channel(
        files, _flat_response(), "XX.T.00.HDH", fs, tmp_path / "a.npz", read=reader
    )
    med = np.nanmedian(native.load(tmp_path / "a.npz")["levels_db"])
    assert med == pytest.approx(10 * np.log10(2 * 10.0**2 / fs), abs=1.5)
    with pytest.raises(ValueError):  # neither code is the channel's own
        native.process_channel(
            files, _flat_response(), "XX.T.20.HDH", fs, tmp_path / "b.npz", read=reader
        )


def test_polygon_cut_lines_are_dropped():
    from compy_pressure.atlas import geography as geo

    ring = np.array(
        [
            [170.0, -10.0],
            [180.0, -10.0],
            [180.0, 10.0],
            [170.0, 10.0],
            [-180.0, 0.0],
            [20.0, 90.0],
            [20.0, -90.0],
        ]
    )
    kept = geo._without_cut_lines(ring)
    assert kept.tolist() == [[170.0, -10.0], [170.0, 10.0]]


def test_shelf_break_ignores_polygon_cut_lines(monkeypatch):
    shp = pytest.importorskip(
        "cartopy.io.shapereader"
    )  # cartopy is an optional dependency of the distances
    Polygon = pytest.importorskip("shapely.geometry").Polygon

    from compy_pressure.atlas import geography as geo

    class Reader:
        def __init__(self, path):
            pass

        def geometries(self):
            return [
                Polygon([(170.0, -10.0), (180.0, -10.0), (180.0, 10.0), (170.0, 10.0)])
            ]

    monkeypatch.setattr(shp, "natural_earth", lambda **kw: "x")
    monkeypatch.setattr(shp, "Reader", Reader)
    pts = geo.line_points("shelf_break", step_km=50.0)
    assert len(pts) and np.all(np.abs(pts[:, 0]) < 179.999)
