"""Unequal-rate APG checks against analytic waveforms, with source preservation."""

import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from obspy import Stream, Trace, UTCDateTime
from obspy.core.inventory import Channel

from compy_pressure import rayleigh_calibration as ray
from compy_pressure.archive import index_waveform_files


@pytest.fixture
def apg_pair(tmp_path):
    rows, channels, paths = [], [], []
    start = UTCDateTime("2018-06-01")
    for code, rate, offset, gain in [("HDH", 125, 0.002, 4), ("HHZ", 100, 0.004, 8)]:
        native_start = start - 100 + offset
        t = np.arange(2200 * rate) / rate + float(native_start - start)
        data = np.zeros_like(t)
        for f in [0.03, 0.045, 0.065]:
            if code == "HDH":
                data += (
                    2
                    * np.pi
                    * f
                    * ray.water_ratio(f, 4000, 1000, 1028, 1500)
                    / 0.67
                    * np.cos(2 * np.pi * f * t)
                )
            else:
                data += np.sin(2 * np.pi * f * t)
        if code == "HHZ":
            # Would alias into .05 Hz if decimated without an anti-alias filter.
            data += 100 * np.sin(2 * np.pi * 1.95 * t)
        trace = Trace(
            gain * data,
            header=dict(
                network="XX",
                station="TEST",
                location="",
                channel=code,
                sampling_rate=rate,
                starttime=native_start,
            ),
        )
        path = tmp_path / f"XX.TEST..{code}.D.2018.151"
        trace.write(str(path), format="MSEED")
        paths.append(path)
        rows.append(pd.Series(dict(Seed_id=trace.id)))
        channels.append(
            Channel(
                code=code,
                location_code="",
                latitude=0,
                longitude=0,
                elevation=-1000,
                depth=0,
                sample_rate=rate,
                dip=-90 if code == "HHZ" else 0,
                start_date=start - 200,
                end_date=start + 3000,
            )
        )
    return start, rows, channels, paths


@pytest.mark.parametrize("target_rate", [1, 2])
def test_resampled_apg_retains_gain_phase_and_rejects_alias(
    apg_pair, monkeypatch, target_rate
):
    start, rows, channels, paths = apg_pair
    before = [hashlib.sha256(p.read_bytes()).hexdigest() for p in paths]
    records = index_waveform_files(paths)
    with pytest.raises(ValueError, match="Unequal sample rates"):
        ray.load_aligned(records, rows, channels, start, 1800)
    arrays, actual, rate = ray.load_aligned(
        records, rows, channels, start, 1800, target_rate
    )
    assert (
        actual == start
        and rate == target_rate
        and arrays.shape == (2, 1800 * target_rate)
    )
    monkeypatch.setattr(
        ray,
        "checked_response",
        lambda row, channel, f, output: np.full(
            len(f), 4 if channel.code == "HDH" else 8, complex
        ),
    )
    f, coefficients, _ = ray.corrected_fft(arrays, rows, channels, rate, 0.02, 0.08)
    for frequency in [0.03, 0.045, 0.065]:
        i = np.argmin(abs(f - frequency))
        measured = coefficients[0, i] / (2j * np.pi * f[i] * coefficients[1, i])
        expected = ray.water_ratio(f[i], 4000, 1000, 1028, 1500) / 0.67
        assert abs(measured) == pytest.approx(expected, rel=2e-4)
        assert abs(np.angle(measured, deg=True)) < 0.02
    alias = np.argmin(abs(f - 0.05))
    signal = np.argmin(abs(f - 0.045))
    assert abs(coefficients[1, alias] / coefficients[1, signal]) < 0.001
    assert [c.sample_rate for c in channels] == [125, 100]
    assert before == [hashlib.sha256(p.read_bytes()).hexdigest() for p in paths]


def test_resampling_refuses_real_gaps(apg_pair):
    start, rows, channels, paths = apg_pair
    from obspy import read

    trace = read(str(paths[0]))[0]
    left = trace.copy().trim(endtime=start + 500)
    right = trace.copy().trim(starttime=start + 501)
    Stream([left, right]).write(str(paths[0]), format="MSEED")
    with pytest.raises(ValueError, match="gapped"):
        ray.load_aligned(index_waveform_files(paths), rows, channels, start, 1800, 2)


def test_resampling_refuses_epoch_crossing_and_upsampling(apg_pair):
    start, rows, channels, paths = apg_pair
    records = index_waveform_files(paths)
    with pytest.raises(ValueError, match="not upsample"):
        ray.load_aligned(records, rows, channels, start, 1800, 126)
    channels[0].start_date = start - 1
    with pytest.raises(ValueError, match="response epoch"):
        ray.load_aligned(records, rows, channels, start, 1800, 2)


def test_native_rate_roundoff_requires_explicit_tolerance(apg_pair):
    from obspy import read

    start, rows, channels, paths = apg_pair
    trace = read(str(paths[0]))[0]
    trace.stats.sampling_rate = 125.00000762939453
    trace.write(str(paths[0]), format="MSEED")
    records = index_waveform_files(paths)
    with pytest.raises(ValueError, match="ppm tolerance"):
        ray.load_aligned(records, rows, channels, start, 1800, 2)
    arrays, actual, rate = ray.load_aligned(
        records, rows, channels, start, 1800, 2, 0.1
    )
    assert np.isfinite(arrays).all() and actual == start and rate == 2


def test_conflicting_overlaps_are_not_interpolated(apg_pair):
    from obspy import read

    start, rows, channels, paths = apg_pair
    original = read(str(paths[0]))[0]
    conflict = original.copy().trim(starttime=start + 300, endtime=start + 400)
    conflict.data += 1
    Stream([original, conflict]).write(str(paths[0]), format="MSEED")
    with pytest.raises(ValueError, match="Conflicting native overlap"):
        ray.load_aligned(index_waveform_files(paths), rows, channels, start, 1800, 2)


@pytest.mark.parametrize("offset", [-0.01, 0.01])
def test_native_fractional_boundary_is_not_interpolated(monkeypatch, offset):
    left = UTCDateTime("2018-06-01")
    headers = dict(
        network="XX",
        station="TEST",
        location="",
        channel="HDH",
        sampling_rate=50,
        starttime=left,
    )
    first = Trace(np.arange(100, dtype=float), header=headers)
    second = Trace(
        np.arange(100, 200, dtype=float),
        header=headers | {"starttime": left + 2 + offset},
    )
    monkeypatch.setattr(ray, "read", lambda *a, **k: Stream([first, second]))
    records = [
        SimpleNamespace(
            path=Path("synthetic.mseed"), seed_id=first.id, start=left, end=left + 4
        )
    ]
    with pytest.raises(ValueError, match="gapped|Fractional"):
        ray.native_resampling_samples(
            records, pd.Series({"Seed_id": first.id}), left, 3, 50, 0
        )
