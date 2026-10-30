"""Read adjacent daily files while keeping genuine gaps out of calculations."""

import numpy as np
import pytest
from obspy import Trace, UTCDateTime

from compy_pressure.archive import discover_waveforms, read_window


def test_downloaded_daily_names_and_header_identity(tmp_path):
    from compy_pressure.network import scan_files

    trace = Trace(np.arange(60, dtype=np.int32))
    trace.stats.update(
        dict(
            network="YL",
            station="A02W",
            location="",
            channel="LDH",
            sampling_rate=1,
            starttime=UTCDateTime("2010-03-11"),
        )
    )
    path = tmp_path / "YL.A02W.--.LDH.20100311.mseed"
    trace.write(str(path), format="MSEED")
    seed = "YL.A02W..LDH"
    assert scan_files(tmp_path, networks=["YL"], stations=["A02W"]) == {seed: [path]}
    records = discover_waveforms(tmp_path, seed_ids={seed})
    assert records[0].seed_id == seed
    np.testing.assert_array_equal(
        read_window(records, trace.stats.starttime, trace.stats.starttime + 60)[0].data,
        trace.data,
    )
    trace.stats.location = "01"
    trace.write(str(path), format="MSEED")
    with pytest.raises(ValueError, match="Filename/header SEED ID mismatch"):
        discover_waveforms(tmp_path)


def test_sds_scan_skips_unselected_network_folders(tmp_path, monkeypatch):
    from compy_pressure import network

    for net in ("YV", "7D"):
        folder = tmp_path / "2012" / net / "STA" / "BDH.D"
        folder.mkdir(parents=True)
        (folder / f"{net}.STA..BDH.D.2012.001").write_bytes(b"")
    walked = []
    original = network._walk_files

    def spy(root, networks=None):
        walked.append(root.name)
        yield from original(root, networks)

    monkeypatch.setattr(network, "_walk_files", spy)
    files = network.scan_files(tmp_path, networks=["7D"])
    assert files == {
        "7D.STA..BDH": [tmp_path / "2012/7D/STA/BDH.D/7D.STA..BDH.D.2012.001"]
    }
    assert "YV" not in walked
    assert set(network.scan_files(tmp_path)) == {"YV.STA..BDH", "7D.STA..BDH"}


def test_filename_location_00_resolves_to_empty_header_location(tmp_path):
    from compy_pressure.archive import index_waveform_files
    from compy_pressure.network import apply_location_aliases, scan_files

    trace = Trace(np.arange(60, dtype=np.int32))
    trace.stats.update(dict(network="7D", station="J28A", channel="BDH"))
    folder = tmp_path / "2012/7D/J28A/BDH.D"
    folder.mkdir(parents=True)
    path = folder / "7D.J28A.00.BDH.D.2012.001"
    trace.write(str(path), format="MSEED")
    files = apply_location_aliases(scan_files(tmp_path), {"7D.J28A..BDH"})
    assert files == {"7D.J28A..BDH": [path]}
    assert index_waveform_files(files["7D.J28A..BDH"])[0].seed_id == "7D.J28A..BDH"
    # A real 00 channel in the inventory keeps its own identity.
    assert set(apply_location_aliases(scan_files(tmp_path), {"7D.J28A.00.BDH"})) == {
        "7D.J28A.00.BDH"
    }
    trace.stats.location = "01"
    trace.write(str(path), format="MSEED")
    with pytest.raises(ValueError, match="Filename/header SEED ID mismatch"):
        index_waveform_files([path])


def test_fast_index_matches_full_header_read(tmp_path):
    from obspy import Stream, read

    from compy_pressure.archive import _record_span

    start = UTCDateTime("2020-01-01")
    traces = []
    for offset, length in [(0, 20000), (3000, 20000)]:  # second trace after a gap
        trace = Trace(np.arange(length, dtype=np.int32) % 977)
        trace.stats.update(
            dict(network="XX", station="T", channel="BDH", sampling_rate=10)
        )
        trace.stats.starttime = start + offset
        traces.append(trace)
    path = tmp_path / "XX.T..BDH.D.2020.001"
    Stream(traces).write(str(path), format="MSEED", reclen=512, encoding="STEIM2")
    span = _record_span(path, "XX.T..BDH")
    full = read(str(path), headonly=True)
    assert span.seed_id == "XX.T..BDH" and span.sample_rate == 10
    assert span.start == min(t.stats.starttime for t in full)
    assert span.end == max(t.stats.endtime for t in full)
    path.write_bytes(path.read_bytes() + b"\0" * 100)  # irregular: full read
    assert _record_span(path, "XX.T..BDH") is None


def test_header_times_cross_file_boundaries_and_gaps_stay_missing(tmp_path):
    start = UTCDateTime("2020-01-01")
    for day, offset, length in [(1, -2, 12), (2, 10, 10), (3, 25, 10)]:
        trace = Trace(np.arange(length, dtype=np.int32))
        trace.stats.update(
            dict(
                network="XX",
                station="TEST",
                location="",
                channel="BHZ",
                sampling_rate=1,
                starttime=start + offset,
            )
        )
        trace.write(str(tmp_path / f"XX.TEST..BHZ.D.2020.{day:03}"), format="MSEED")
    records = discover_waveforms(tmp_path)
    assert len(records) == 3
    assert records[0].start == start - 2
    stream = read_window(records, start, start + 30)
    assert [(trace.stats.starttime - start, trace.stats.npts) for trace in stream] == [
        (0, 20),
        (25, 5),
    ]
    assert all(not np.ma.isMaskedArray(trace.data) for trace in stream)


def test_strict_grid_refuses_rounding_of_overlapping_records(tmp_path):
    start = UTCDateTime("2020-01-01")
    for i, (offset, values) in enumerate(
        [(0, np.arange(1500)), (29.9984, np.arange(1500, 3000))]
    ):
        trace = Trace(values.astype(np.int32))
        trace.stats.update(
            dict(
                network="XX",
                station="TEST",
                channel="BHZ",
                sampling_rate=50,
                starttime=start + offset,
            )
        )
        trace.write(str(tmp_path / f"XX.TEST..BHZ.D.2020.{i + 1:03}"), format="MSEED")
    records = discover_waveforms(tmp_path)
    with pytest.raises(ValueError, match="sample-grid offset"):
        read_window(records, start, start + 59, strict_sample_grid=True)


def test_strict_grid_keeps_exact_duplicate_overlaps_and_continuity(tmp_path):
    start = UTCDateTime("2020-01-01")
    for i, (offset, values) in enumerate(
        [(0, np.arange(2000)), (30, np.arange(1500, 3500))]
    ):
        trace = Trace(values.astype(np.int32))
        trace.stats.update(
            dict(
                network="XX",
                station="TEST",
                channel="BHZ",
                sampling_rate=50,
                starttime=start + offset,
            )
        )
        trace.write(str(tmp_path / f"XX.TEST..BHZ.D.2020.{i + 1:03}"), format="MSEED")
    stream = read_window(
        discover_waveforms(tmp_path), start, start + 60, strict_sample_grid=True
    )
    assert len(stream) == 1
    np.testing.assert_array_equal(stream[0].data, np.arange(3000))


def test_conflicting_overlap_becomes_missing_not_averaged_or_chosen(tmp_path):
    start = UTCDateTime("2020-01-01")
    for day, offset, values in [(1, 0, np.arange(40)), (2, 30, np.arange(40) + 1000)]:
        trace = Trace(values.astype(np.int32))
        trace.stats.update(
            dict(
                network="XX",
                station="TEST",
                channel="BHZ",
                sampling_rate=1,
                starttime=start + offset,
            )
        )
        trace.write(str(tmp_path / f"XX.TEST..BHZ.D.2020.{day:03}"), format="MSEED")
    stream = read_window(discover_waveforms(tmp_path), start, start + 70)
    # Samples 30-39 disagree between files; neither version is retained.
    assert [(tr.stats.starttime - start, tr.stats.npts) for tr in stream] == [
        (0, 30),
        (40, 30),
    ]
    np.testing.assert_array_equal(stream[0].data, np.arange(30))
    np.testing.assert_array_equal(stream[1].data, np.arange(10, 40) + 1000)


def test_empty_daily_file_is_skipped_but_corrupt_file_fails(tmp_path):
    from compy_pressure.archive import index_waveform_files

    trace = Trace(np.arange(60, dtype=np.int32))
    trace.stats.update(dict(network="IM", station="CY301", channel="EDH"))
    good = tmp_path / "IM.CY301..EDH.D.2025.160"
    trace.write(str(good), format="MSEED")
    empty = tmp_path / "IM.CY301..EDH.D.2025.162"
    empty.write_bytes(b"")
    assert [r.path for r in index_waveform_files([good, empty])] == [good]
    empty.write_bytes(b"not miniSEED")
    with pytest.raises(RuntimeError, match="Cannot read MiniSEED header"):
        index_waveform_files([good, empty])
