"""Tilt correction accepts float32 traces (stage A and MiniSEED FLOAT32 output)."""

import numpy as np
from obspy import Stream, Trace, UTCDateTime

import compy.preprocessing as pre


def _stream(dtype):
    rng = np.random.default_rng(0)
    n = int(2 * 3600 * 1.25)
    h1, h2 = rng.standard_normal(n), rng.standard_normal(n)
    z = 0.2 * rng.standard_normal(n) + 0.05 * h1 - 0.03 * h2  # tilt leakage
    traces = []
    for comp, data in (("Z", z), ("1", h1), ("2", h2), ("H", rng.standard_normal(n))):
        header = dict(network="XX", station="TEST", channel=f"BH{comp}" if comp != "H" else "BDH",
                      sampling_rate=1.25, starttime=UTCDateTime(2020, 1, 1))  # fmt: skip
        tr = Trace(data.astype(dtype), header=header)
        traces.append(tr)
    return Stream(traces)


def test_correct_tilt_float32_matches_float64():
    st32 = _stream(np.float32)
    st64 = st32.copy()
    for tr in st64:
        tr.data = tr.data.astype(np.float64)
    out32, az32, ang32, vr32 = pre.correct_tilt(st32, 1, plot=False)
    out64, az64, ang64, vr64 = pre.correct_tilt(st64, 1, plot=False)
    np.testing.assert_array_equal(ang32, ang64)
    np.testing.assert_array_equal(vr32, vr64)
    for a, b in zip(out32, out64):
        np.testing.assert_array_equal(a.data, b.data)
    assert st32[0].data.dtype == np.float32  # input not modified


def test_correct_tilt_skip_failed_keeps_window_uncorrected(monkeypatch):
    st = _stream(np.float64)
    real = pre.tiskit.CleanRotator
    calls = {"n": 0}

    def flaky(window, **kwargs):
        calls["n"] += 1
        if calls["n"] == 2:
            raise ValueError("angle > rot_limit")
        return real(window, **kwargs)

    monkeypatch.setattr(pre.tiskit, "CleanRotator", flaky)
    try:
        pre.correct_tilt(st, 1, plot=False, clean_horizontals=False)
        raise AssertionError("should raise without skip_failed")
    except RuntimeError:
        pass
    calls["n"] = 0
    with np.testing.suppress_warnings() as sup:
        sup.filter(UserWarning)
        out, az, ang, vr = pre.correct_tilt(
            st, 1, plot=False, clean_horizontals=False, skip_failed=True
        )
    assert np.isfinite(ang[0]) and np.isnan(ang[1]) and np.isnan(vr[1])
    # The failed hour passes through unchanged.
    z_in = st.select(component="Z")[0]
    z_out = out.select(component="Z")[0]
    i0 = int(3600 * 1.25)
    np.testing.assert_array_equal(z_out.data[i0 : 2 * i0], z_in.data[i0 : 2 * i0])


def test_correct_tilt_spectral_window_option():
    st = _stream(np.float64)
    out, _, ang, _ = pre.correct_tilt(st, 1, plot=False, spectral_window_s=600)
    assert np.all(np.isfinite(ang)) and len(out) == 4
