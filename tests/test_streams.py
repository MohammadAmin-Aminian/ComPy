import numpy as np
import pytest
from obspy import Stream, Trace, UTCDateTime

import compy
from compy_processing import spectra


def make_stream(n=100, fs=10):
    return Stream(
        [
            Trace(
                np.arange(n, dtype=float),
                header={
                    "channel": channel,
                    "sampling_rate": fs,
                    "starttime": UTCDateTime(0),
                    "station": "TEST",
                    "network": "XX",
                },
            )
            for channel in ("BHZ", "BDH")
        ]
    )


def test_windows_cover_samples_once_and_copy_input():
    original = make_stream()
    parts = compy.split_stream(original, 2)
    assert len(parts) == 5
    assert all(len(part[0]) == 20 for part in parts)
    np.testing.assert_equal(
        np.concatenate([part[0].data for part in parts]), original[0].data
    )
    parts[0][0].data[0] = -99
    assert original[0].data[0] == 0


@pytest.mark.parametrize(
    "length,overlap", [(0, 0), (1, 1), (1, 2), (1, -1), (1, np.nan)]
)
def test_bad_window_parameters(length, overlap):
    with pytest.raises(ValueError):
        compy.cut_stream_with_overlap(make_stream(), length, overlap)


def test_overlap_and_common_coverage():
    st = make_stream()
    st[1].trim(st[1].stats.starttime + 1, st[1].stats.endtime - 1)
    parts = compy.cut_stream_with_overlap(st, 2, 1)
    assert parts[0][0].stats.starttime == st[1].stats.starttime
    assert parts[-1][0].stats.endtime <= st[1].stats.endtime
    assert parts[1][0].stats.starttime - parts[0][0].stats.starttime == 1
    assert all(len(part[0]) == len(part[1]) == 20 for part in parts)


def test_save_preserves_fractional_second_samples(tmp_path):
    from obspy import read

    st = make_stream(n=105)
    paths = compy.split_and_save_stream(st, 2 / 60, tmp_path)
    restored = np.concatenate(
        [read(str(path)).select(channel="BHZ")[0].data for path in paths]
    )
    np.testing.assert_equal(restored, st[0].data)
    with pytest.raises(FileExistsError):
        compy.split_and_save_stream(st, 2 / 60, tmp_path)


def test_sliding_window_tapers_integer_data_in_float():
    windows, n = compy.sliding_window(np.ones(8, dtype=int), 4, 2)
    assert n == 3
    np.testing.assert_allclose(windows[0], [0, 0.75, 0.75, 0])
    for ws, ss in [(0, 1), (4, 0), (1.5, 1)]:
        with pytest.raises(ValueError):
            compy.sliding_window(np.ones(8), ws, ss)


def test_spectral_ratio_and_coherence_known_linear_transfer():
    st = make_stream(8192, 2)
    p = np.random.default_rng(4).normal(size=8192)
    st.select(channel="BDH")[0].data = p
    st.select(channel="BHZ")[0].data = p * 3
    f, zz, pp, zp, coh = spectra(st, 1024)
    np.testing.assert_allclose(zz, pp * 9, rtol=1e-12)
    np.testing.assert_allclose(zp / pp, 3, atol=1e-12)
    np.testing.assert_allclose(coh, 1, atol=1e-12)
    assert f[0] == 0


def test_empty_quality_selection_has_clear_error():
    st = make_stream(7200, 2)
    rng = np.random.default_rng(12)
    for tr in st:
        tr.data = rng.normal(size=len(tr))
    with pytest.raises(ValueError, match="no windows passed"):
        compy.Calculate_Compliance_beta(
            st, time_window=1, depth=4000, nseg=1024, plot=False
        )


def test_compliance_known_transfer_and_gravity_units():
    st = make_stream(7200, 2)
    p = np.random.default_rng(3).normal(0, 100, 7200)
    st.select(channel="BDH")[0].data = p
    st.select(channel="BHZ")[0].data = p * 1e-7
    curves, coherences, windows, fc, f, scatter = compy.Calculate_Compliance_beta(
        st, time_window=1, depth=4000, gain_factor=1, nseg=1024, plot=False
    )
    omega2 = (2 * np.pi * fc) ** 2
    gravity, _, _ = compy.gravitational_attraction(np.ones_like(fc), 4000, fc)
    expected = (
        compy.wavenumber(2 * np.pi * fc, 4000)
        * (omega2 * 1e-7 + gravity)
        / (omega2 + 3.07e-6)
    )
    np.testing.assert_allclose(curves[0], expected, rtol=1e-12)
    assert len(windows) == 1
    assert np.all(scatter == 0)


def test_rotation_reports_failure_instead_of_returning_unclean_data(monkeypatch):
    def fail(*args, **kwargs):
        raise ValueError("synthetic failure")

    monkeypatch.setattr(compy.tiskit, "CleanRotator", fail)
    st = make_stream(7200, 2)
    with pytest.raises(RuntimeError, match="rotation failed in window 0"):
        compy.Rotate(st, time_window=1, plot=False)
    assert len(st[0]) == 7200


def test_rotation_with_real_tiskitpy_api():
    rng = np.random.default_rng(3)
    horizontal = rng.normal(size=7200)
    traces = []
    for channel, data in [
        ("BH1", horizontal),
        ("BH2", rng.normal(size=7200)),
        ("BHZ", 0.01 * horizontal + 0.001 * rng.normal(size=7200)),
        ("BDH", rng.normal(size=7200)),
    ]:
        traces.append(
            Trace(
                data,
                header={
                    "channel": channel,
                    "sampling_rate": 2.0,
                    "starttime": UTCDateTime(0),
                },
            )
        )
    stream = Stream(traces)
    original = stream.copy()
    rotated, azimuth, angle, variance = compy.Rotate(stream, 1, plot=False)
    assert stream == original
    assert len(rotated.select(component="Z")[0]) == 7200
    assert variance.shape == (1,)
    assert 0 < variance[0] < 1
