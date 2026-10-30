from types import SimpleNamespace

import numpy as np
import pytest
from obspy import Stream, Trace, UTCDateTime

from compy_pressure import local_p_check as lc
from compy_pressure import local_p_physics as lp
from compy_pressure.local_p_reference import simulate_plane_wave


def flat_responses(pressure_gain, vertical_gain):
    """Counts per Pa and counts per (m/s), flat and real (signs allowed)."""

    def evaluate(f):
        rp = np.full(f.shape, pressure_gain, dtype=complex)
        rv = np.full(f.shape, vertical_gain, dtype=complex)
        return rp, rv, np.ones(f.shape), np.ones(f.shape)

    return evaluate


def synthetic_event(depth, dz, peak, duration, source_below, pad_seconds, seed):
    """FDTD seafloor p and v on a 100 Hz grid with a quiet pre-event pad."""
    t, p, v, onset = simulate_plane_wave(
        depth=depth,
        dz=dz,
        duration=duration,
        peak_frequency=peak,
        source_below_seafloor=source_below,
        solid_speed=3000.0,
        solid_density=2200.0,
    )
    rate = 100.0
    grid = np.arange(0, t[-1], 1 / rate)
    p, v = np.interp(grid, t, p), np.interp(grid, t, v)
    pad = int(pad_seconds * rate)
    rng = np.random.default_rng(seed)
    scale_p, scale_v = 1e-4 * np.max(np.abs(p)), 1e-4 * np.max(np.abs(v))
    p = np.concatenate([np.zeros(pad), p]) + scale_p * rng.standard_normal(pad + p.size)
    v = np.concatenate([np.zeros(pad), v]) + scale_v * rng.standard_normal(pad + v.size)
    return rate, p, v, int(round(onset * rate)) + pad


@pytest.mark.parametrize("true_gain,sign", [(0.7, 1.0), (1.4, -1.0)])
def test_forward_model_recovers_gain_and_polarity_pre_multiple(true_gain, sign):
    """Pressure recorded with a wrong nominal sensitivity: p_nominal = p_true/g."""
    depth = 3000.0
    rate, p, v, i_on = synthetic_event(depth, 4.0, 6.0, 14.0, 800.0, 8.0, 3)
    k_p, k_v = 500.0, 6e8
    # Nominal response says k_p counts/Pa, the instrument delivers k_p/g * sign.
    p_counts = sign * (k_p / true_gain) * p
    v_counts = k_v * v
    n = int(round(0.9 * lp.round_trip_time(depth) * rate))
    signal = slice(i_on - 20, i_on - 20 + n)
    noise = slice(i_on - 120 - n, i_on - 120)
    bands, arrays, _ = lc.analyse_event(
        p_counts,
        v_counts,
        signal,
        noise,
        rate,
        flat_responses(k_p, k_v),
        "pre_multiple",
        depth,
        1.0,
        lc.DEFAULTS,
    )
    good = [b for b in bands if "gain_noise_corrected" in b and b["band_hz"][0] >= 2]
    assert good
    for band in good:
        assert band["gain_noise_corrected"] == pytest.approx(true_gain, rel=0.05)
        assert band["gain_h2_lower"] <= band["gain_noise_corrected"] * 1.01
        assert band["polarity"] == int(sign)
        assert band["gain_model_low"] < true_gain < band["gain_model_high"]


def test_forward_model_recovers_gain_steady_state_shallow():
    depth = 150.0
    rate, p, v, i_on = synthetic_event(depth, 1.0, 1.5, 40.0, 600.0, 30.0, 4)
    k_p, k_v, true_gain = 100.0, 1.25e8, 0.8
    p_counts = (k_p / true_gain) * p
    v_counts = k_v * v
    n = int(20 * rate)
    signal = slice(i_on - 20, i_on - 20 + n)
    noise = slice(i_on - 120 - n, i_on - 120)
    bands, arrays, diag = lc.analyse_event(
        p_counts,
        v_counts,
        signal,
        noise,
        rate,
        flat_responses(k_p, k_v),
        "steady_state",
        depth,
        1.0,
        lc.DEFAULTS,
    )
    assert diag["model_passband_hz"][3] < lp.resonance_frequencies(depth)[0]
    good = [b for b in bands if "gain_noise_corrected" in b]
    assert good
    for band in good:
        assert band["gain_noise_corrected"] == pytest.approx(true_gain, rel=0.05)
        assert band["polarity"] == 1


def test_model_transfer_convention_and_limits():
    o = lc.DEFAULTS
    f = np.array([0.01, 0.5])
    steady = lc.model_transfer(f, "steady_state", 4000.0, 1.0, o)
    # p/v = i 2 pi f rho h F (numpy d/dt = +i 2 pi f), F -> 1 at low frequency.
    theta = lp.water_phase(0.01, 4000.0)
    expected = 1j * 2 * np.pi * 0.01 * 1028 * 4000 * np.tan(theta) / theta
    assert steady[0] == pytest.approx(expected, rel=1e-4)
    assert np.isnan(steady[1])  # beyond the first-pole margin at 4 km
    flat = lc.model_transfer(f, "pre_multiple", 4000.0, 0.96, o)
    assert np.allclose(flat, 1028 * 1500 / 0.96)


def test_choose_window_regimes():
    o = dict(lc.DEFAULTS)
    deep = SimpleNamespace(p_travel_s=20.0, s_travel_s=40.0)
    regime, length, tau_min, _ = lc.choose_window(deep, 0.0, 4000.0, o)
    assert regime == "pre_multiple"
    assert length <= o["pre_multiple_safety"] * tau_min + o["lead_s"] + 1e-9
    shallow = SimpleNamespace(p_travel_s=10.0, s_travel_s=25.0)
    regime, length, _, tau_max = lc.choose_window(shallow, 0.0, 152.0, o)
    assert regime == "steady_state" and length >= 10 * tau_max
    close = SimpleNamespace(p_travel_s=2.0, s_travel_s=3.0)
    with pytest.raises(lc.EventRejected) as info:
        lc.choose_window(close, 0.0, 4000.0, o)
    assert info.value.reason == "s_too_close"
    # 500 m: too long for the pre-multiple gate, too short for steady state.
    mid = SimpleNamespace(p_travel_s=5.0, s_travel_s=9.5)
    with pytest.raises(lc.EventRejected) as info:
        lc.choose_window(mid, 0.0, 500.0, o)
    assert info.value.reason == "window_regime_mixed"


def test_duration_quantum_holds_integer_samples():
    q = lc.duration_quantum([125.0, 100.0, 62.5, 50.0, 40.0])
    for rate in (125.0, 100.0, 62.5, 50.0, 40.0):
        assert q * rate == pytest.approx(round(q * rate))
    assert lc.duration_quantum([100.0]) == pytest.approx(0.01)


def test_envelope_lag_sign_convention():
    rate = 50.0
    rng = np.random.default_rng(5)
    burst = np.zeros(6000)
    burst[2000:2400] = rng.standard_normal(400) * np.hanning(400)
    v = burst + 0.01 * rng.standard_normal(6000)
    p = np.roll(burst, -int(12.3 * rate)) + 0.01 * rng.standard_normal(6000)
    lag, corr = lc.envelope_lag(p, v, rate, (0.5, 10.0), 30.0)
    assert lag == pytest.approx(-12.3, abs=0.2)  # pressure earlier: negative lag
    assert corr > 0.5


def test_measure_onset_and_clipping():
    rate = 100.0
    rng = np.random.default_rng(6)
    x = np.concatenate([0.01 * rng.standard_normal(3000), rng.standard_normal(1000)])
    onset, snr = lc.measure_onset(x, rate, 1000, 3500, (2.0, 10.0), 3.0)
    assert abs(onset - 3000) < 15 and snr > 10
    clipped = np.sin(np.linspace(0, 20, 500)) * 2
    assert lc.flat_topped(np.clip(clipped, -1, 1))
    assert not lc.flat_topped(clipped)


def test_station_summary_outcomes():
    o = dict(lc.DEFAULTS)
    band = {
        "band_hz": [1.0, 2.0],
        "gain_noise_corrected": 1.0,
        "gain_h1_upper": 1.05,
        "gain_h2_lower": 0.95,
        "gain_model_low": 0.9,
        "gain_model_high": 1.1,
        "phase_deg": 3.0,
        "polarity": 1,
    }
    events = [
        {
            "event_id": f"e{i}",
            "status": "accepted",
            "regime": "steady_state",
            "bands": [band],
        }
        for i in range(3)
    ]
    assert lc.station_summary(events, o, 3)["outcome"] == "conditional_gain_candidate"
    gated = [dict(e, regime="pre_multiple") for e in events]
    result = lc.station_summary(gated, o, 3)
    assert result["outcome"] == "diagnostic_only" and result["score"] is None
    assert lc.station_summary([], o, 5)["outcome"] == "unresolved"
    assert lc.station_summary([], o, 0)["outcome"] == "no_suitable_events"


def test_catalog_csv_parsing(tmp_path):
    path = tmp_path / "cat.csv"
    path.write_text(
        "time,latitude,longitude,depth,mag,magType,net,id,type,horizontalError,depthError\n"
        "2018-06-01T18:13:07.100Z,55.1,-156.2,30.5,3.5,ml,ak,ak1,earthquake,2.1,\n"
    )
    (event,) = lc.load_local_catalog(path)
    assert event["time"] == UTCDateTime("2018-06-01T18:13:07.100Z")
    assert event["horizontal_error_km"] == 2.1 and event["depth_error_km"] is None


def test_day_file_index_uses_headers_and_prefers_sds(tmp_path):
    start = UTCDateTime("2018-07-20T00:00:00")
    trace = Trace(np.arange(100, dtype=np.int32), header={"sampling_rate": 1.0})
    trace.stats.network, trace.stats.station = "XO", "LT02"
    trace.stats.location, trace.stats.channel = "", "HDH"
    trace.stats.starttime = start
    sds = tmp_path / "HDH.D"
    sds.mkdir()
    daily = tmp_path / "waveforms"
    daily.mkdir()
    Stream([trace]).write(str(sds / "XO.LT02.00.HDH.D.2018.201.mseed"), format="MSEED")
    Stream([trace]).write(str(daily / "XO.LT02.--.HDH.20180720.mseed"), format="MSEED")
    index = lc.DayFileIndex([tmp_path], "XO.LT02..HDH")
    records = index.records(start, start + 50)
    assert len(records) == 1 and records[0].path.parent.name == "HDH.D"
    assert records[0].seed_id == "XO.LT02..HDH"
    assert "2018.201" in index.duplicate_layout_days
    assert index.location_mismatch is False  # header location is blank


def test_output_protection(tmp_path):
    with pytest.raises(ValueError, match="LOCAL_P"):
        lc.assert_fresh_output(tmp_path / "out", [])
    escape = tmp_path / "proj/results/LOCAL_P/../../NETWORKS/run"
    with pytest.raises(ValueError, match="LOCAL_P"):
        lc.assert_fresh_output(escape, [tmp_path / "proj/results"])
    allowed = tmp_path / "proj/results/LOCAL_P/run1"
    assert (
        lc.assert_fresh_output(allowed, [tmp_path / "proj/results"])
        == allowed.resolve()
    )
    allowed.mkdir(parents=True)
    (allowed / "x").write_text("1")
    with pytest.raises(FileExistsError):
        lc.assert_fresh_output(allowed, [])
