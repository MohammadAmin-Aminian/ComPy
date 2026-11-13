"""Crashes found while running every OBS network of the archive through stages A and B."""

import importlib
import sys
from pathlib import Path

import numpy as np
import pytest
from obspy import Stream, Trace, UTCDateTime


@pytest.fixture
def stage_a(tmp_path, monkeypatch):
    monkeypatch.setenv("COMPY_RUN_DIR", str(tmp_path))
    monkeypatch.setenv("COMPY_SDS", str(tmp_path / "sds"))
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "pipeline"))
    for name in ("config", "a_preprocess"):
        sys.modules.pop(name, None)
    module = importlib.import_module("a_preprocess")
    yield module
    for name in ("config", "a_preprocess"):
        sys.modules.pop(name, None)


def trace(data, rate=50.0, start=UTCDateTime(2020, 1, 1, 0, 0, 0.01)):
    return Trace(
        np.asarray(data, float),
        {"sampling_rate": rate, "starttime": start, "channel": "BHZ"},
    )


def test_piece_ending_before_next_grid_point_is_skipped(stage_a):
    assert stage_a.decimate_segment(trace(np.arange(5))) is None


def test_spans_less_than_a_sample_apart_do_not_lose_the_chunk(stage_a):
    # XE.CC10, 8 June 2018: the bad spans of two channels were 2 ms apart, and
    # slicing the gap between them raised ValueError for the whole 6-day chunk.
    tr = trace(np.arange(50 * 400, dtype=float))
    t0 = tr.stats.starttime
    spans = [[t0 + 100, t0 + 200, {"BHZ:dead"}], [t0 + 200.002, t0 + 300, {"BDH:clip"}]]
    pieces = stage_a.cut_spans(tr, spans)
    assert [p.stats.npts for p in pieces] == [50 * 100, 50 * 100]
    assert pieces[0].stats.endtime < t0 + 100 and pieces[1].stats.starttime >= t0 + 300


def test_invalid_day_files_are_skipped(stage_a, tmp_path):
    folder = stage_a.SDS / "2013" / stage_a.NETWORK / "CC11" / f"{stage_a.Z}.D"
    folder.mkdir(parents=True)
    for day in ("365", "366"):
        (folder / f"x.D.2013.{day}").touch()
    (folder / "x.D.2013.001.part").touch()
    monkeypatch_years = stage_a.YEARS
    stage_a.YEARS = (2013,)
    try:
        assert stage_a.station_days("CC11") == [UTCDateTime(year=2013, julday=365)]
    finally:
        stage_a.YEARS = monkeypatch_years


def test_clip_threshold_follows_the_digitizer(stage_a):
    width = 100
    quiet = np.random.default_rng(0).normal(0, 1000, 10 * width)
    assert stage_a.clip_threshold(quiet, width) == stage_a.CLIP_COUNTS
    assert stage_a.clip_threshold(quiet + 2**28, width) == stage_a.CLIP_COUNTS_32
    # A 24-bit channel at full scale stays a 24-bit channel, and so is flagged as clipped.
    assert (
        stage_a.clip_threshold(np.full(10 * width, 2.0**23), width)
        == stage_a.CLIP_COUNTS
    )


def test_absolute_gauge_is_not_cut_as_clipped(stage_a):
    noisy = np.random.default_rng(0).normal(2**28, 1000, 200000)
    assert stage_a.bad_spans(trace(noisy)) == []


def test_written_cleaned_stream_has_no_masked_samples():
    a = trace(np.ones(4000), 1.25)
    b = trace(np.ones(4000), 1.25, a.stats.endtime + 800)
    merged = Stream([a, b]).merge(method=0)
    assert np.ma.is_masked(merged[0].data)
    split = merged.split()
    assert len(split) == 2 and not any(np.ma.is_masked(tr.data) for tr in split)
