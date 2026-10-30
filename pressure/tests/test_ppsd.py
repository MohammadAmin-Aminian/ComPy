"""Validate physical response correction, epoch isolation, and restart behavior."""

import json
from pathlib import Path

import numpy as np
import pytest
from obspy import Trace, UTCDateTime, read_inventory
from obspy.core.inventory import Channel, Inventory, Network, Site, Station
from obspy.core.inventory.response import Response
from obspy.signal import PPSD

from compy_pressure import ppsd as processing
from compy_pressure.archive import discover_waveforms
from compy_pressure.inventory import (
    build_master_dataset,
    load_master_dataset,
    save_master_dataset,
)
from compy_pressure.ppsd import process_epoch, response_handling, run_processing


@pytest.fixture
def dataset(tmp_path):
    start = UTCDateTime("2020-01-01")
    inventory_dir = tmp_path / "inventories"
    inventory_dir.mkdir()
    channel = Channel(
        code="BDH",
        location_code="00",
        latitude=1,
        longitude=2,
        elevation=-1000,
        depth=0,
        sample_rate=8,
        start_date=start,
        end_date=start + 4096,
        response=Response.from_paz(
            zeros=[],
            poles=[],
            stage_gain=4,
            input_units="PA",
            output_units="COUNTS",
        ),
    )
    station = Station(
        code="TEST",
        latitude=1,
        longitude=2,
        elevation=-1000,
        site=Site(name="Test"),
        channels=[channel],
    )
    Inventory([Network(code="XX", stations=[station])], source="Test").write(
        str(inventory_dir / "test.xml"),
        format="STATIONXML",
    )
    waveform_dir = tmp_path / "waveforms"
    waveform_dir.mkdir()
    data = np.random.default_rng(7).integers(-10000, 10000, 16384, dtype=np.int32)
    trace = Trace(data)
    trace.stats.update(
        dict(
            network="XX",
            station="TEST",
            location="00",
            channel="BDH",
            starttime=start,
            sampling_rate=8,
        )
    )
    # Split continuous data at a file boundary to test chunk overlap as well.
    trace.slice(start, start + 1024 - 1 / 8).write(
        str(waveform_dir / "XX.TEST.00.BDH.D.2020.001"),
        format="MSEED",
    )
    trace.slice(start + 1024, trace.stats.endtime).write(
        str(waveform_dir / "XX.TEST.00.BDH.D.2020.002"),
        format="MSEED",
    )
    master = build_master_dataset(inventory_dir).astype(object).fillna("")
    row = master.iloc[0]
    records = discover_waveforms(waveform_dir)
    return row, records, trace, tmp_path / "results"


def test_pressure_response_is_applied_once_and_npz_reloads(dataset):
    row, records, trace, output = dataset
    result = process_epoch(
        row,
        records,
        output,
        ppsd_length=512,
        minimum_frequency=0.02,
        chunk_seconds=1024,
    )
    corrected = PPSD.load_npz(result["ppsd_npz_file"])
    raw = PPSD(
        trace.stats,
        {"sensitivity": 1, "poles": [], "zeros": [], "gain": 1},
        special_handling="ringlaser",
        ppsd_length=512,
        period_limits=(0.25, 50),
        db_bins=(-200, 100, 1),
    )
    raw.add(trace)
    assert result["ppsd_status"] == "complete"
    assert corrected.special_handling == "hydrophone"
    assert result["ppsd_units"] == "Pa^2/Hz"
    assert result["ppsd_segments"] == len(raw.times_processed) == 7
    np.testing.assert_allclose(
        np.array(corrected.psd_values) - np.array(raw.psd_values),
        -20 * np.log10(4),
        atol=1e-4,
    )
    provenance = json.loads(open(result["ppsd_provenance_file"]).read())
    assert len(provenance["source_files"]) == 2
    assert read_inventory(
        str(
            output
            / "XX"
            / "TEST"
            / next((output / "XX" / "TEST").iterdir()).name
            / "inventory.xml"
        )
    ).get_response(row.Seed_id, trace.stats.starttime)


def test_resume_extends_pilot_without_duplicate_segments(dataset):
    row, records, trace, output = dataset
    options = dict(ppsd_length=512, minimum_frequency=0.02, chunk_seconds=1024)
    pilot = process_epoch(
        row, records, output, end=trace.stats.starttime + 1024, **options
    )
    assert pilot["ppsd_segments"] == 3
    full = process_epoch(row, records, output, **options)
    repeated = process_epoch(row, records, output, **options)
    assert full["ppsd_segments"] == repeated["ppsd_segments"] == 7
    times = PPSD.load_npz(full["ppsd_npz_file"]).times_processed
    assert len({str(time) for time in times}) == 7
    with pytest.raises(ValueError, match="Saved settings/sources differ"):
        process_epoch(row, records, output, **(options | {"overlap": 0.25}))


def test_signed_pressure_response_preserves_psd_and_complex_polarity(dataset):
    from compy_pressure.paired_spectra import checked_response

    row, records, _, output = dataset
    settings = dict(ppsd_length=512, minimum_frequency=0.02, chunk_seconds=1024)
    positive = process_epoch(row, records, output / "positive", **settings)
    inventory = read_inventory(row.inventory_file)
    channel = inventory[0][0][0]
    f = np.array([0.03, 0.05, 0.07])
    original_response = channel.response.get_evalresp_response_for_frequencies(
        f, output="DEF"
    )
    channel.response.instrument_sensitivity.value *= -1
    channel.response.response_stages[0].stage_gain *= -1
    inventory.write(row.inventory_file, format="STATIONXML")
    assert response_handling(channel) == ("hydrophone", "Pa^2/Hz")
    negative = process_epoch(row, records, output / "negative", **settings)
    assert negative["ppsd_status"] == "complete"
    a = PPSD.load_npz(positive["ppsd_npz_file"])
    b = PPSD.load_npz(negative["ppsd_npz_file"])
    np.testing.assert_allclose(a.psd_values, b.psd_values, atol=1e-5)
    signed_row = row.copy()
    signed_row["ppsd_npz_file"] = negative["ppsd_npz_file"]
    np.testing.assert_allclose(
        checked_response(signed_row, channel, f, "DEF"), -original_response, rtol=1e-12
    )
    saved = read_inventory(
        str(Path(negative["ppsd_npz_file"]).parent / "inventory.xml")
    )
    assert saved[0][0][0].response.instrument_sensitivity.value < 0


@pytest.mark.parametrize("value", [0.0, float("nan"), float("inf"), -float("inf")])
def test_invalid_sensitivity_remains_rejected(dataset, value):
    row, _, _, _ = dataset
    channel = read_inventory(row.inventory_file)[0][0][0]
    channel.response.instrument_sensitivity.value = value
    with pytest.raises(ValueError, match="finite and nonzero"):
        response_handling(channel)


def test_short_fft_and_unknown_units_are_rejected(dataset):
    row, records, _, output = dataset
    with pytest.raises(ValueError, match="too short"):
        process_epoch(row, records, output, ppsd_length=100, minimum_frequency=0.001)
    inventory = read_inventory(row.inventory_file)
    channel = inventory[0][0][0]
    channel.response.instrument_sensitivity.input_units = "HPA"
    with pytest.raises(ValueError, match="Unsupported response input units"):
        response_handling(channel)


def test_seismic_handling_requires_full_response(dataset):
    row, _, _, _ = dataset
    channel = read_inventory(row.inventory_file)[0][0][0]
    channel.response.instrument_sensitivity.input_units = "M/S"
    assert response_handling(channel) == (None, "m^2/s^4/Hz")
    channel.response.response_stages = []
    with pytest.raises(ValueError, match="Full response"):
        response_handling(channel)


def test_read_failure_is_reported_and_retried_on_resume(dataset, monkeypatch):
    row, records, _, output = dataset
    options = dict(ppsd_length=512, minimum_frequency=0.02, chunk_seconds=1024)
    real_reader = processing.read_window
    calls = 0

    def fail_first(*args):
        nonlocal calls
        calls += 1
        if calls == 1:
            raise OSError("Read failed")
        return real_reader(*args)

    monkeypatch.setattr(processing, "read_window", fail_first)
    partial = process_epoch(row, records, output, **options)
    assert partial["ppsd_status"] == "partial"
    assert "Read failed" in partial["ppsd_error"]
    repaired = process_epoch(row, records, output, **options)
    assert repaired["ppsd_status"] == "complete"
    assert repaired["ppsd_segments"] == 7


def test_no_psd_window_crosses_the_response_epoch_end(dataset):
    row, records, trace, output = dataset
    inventory = read_inventory(row.inventory_file)
    inventory[0][0][0].end_date = trace.stats.starttime + 1024
    inventory.write(row.inventory_file, format="STATIONXML")
    row = row.copy()
    row["channel_end_time"] = str(trace.stats.starttime + 1024)
    result = process_epoch(
        row,
        records,
        output,
        ppsd_length=512,
        minimum_frequency=0.02,
        chunk_seconds=1024,
    )
    assert result["ppsd_segments"] == 3
    assert UTCDateTime(result["ppsd_end_time"]) <= inventory[0][0][0].end_date


def test_fractional_file_timestamps_do_not_drop_interior_chunk_windows(dataset):
    row, records, original, output = dataset
    directory = records[0].path.parent
    trace = Trace(np.tile(original.data, 2), header=original.stats.copy())
    start = original.stats.starttime
    for index, phase in enumerate([0, 0.01, 0.003, 0.009]):
        piece = trace.slice(start + index * 1024, start + (index + 1) * 1024 - 1 / 8)
        piece.stats.starttime += phase
        piece.write(
            str(directory / f"XX.TEST.00.BDH.D.2020.{index + 1:03}"),
            format="MSEED",
        )
    result = process_epoch(
        row,
        discover_waveforms(directory),
        output,
        ppsd_length=512,
        minimum_frequency=0.02,
        chunk_seconds=1024,
    )
    times = PPSD.load_npz(result["ppsd_npz_file"]).times_processed
    assert result["ppsd_status"] == "complete"
    assert len(times) >= 14
    assert max(right - left for left, right in zip(times, times[1:])) <= 256 + 1 / 8


def test_parallel_channels_update_the_master_after_workers_finish(dataset):
    row, records, trace, output = dataset
    inventory = read_inventory(row.inventory_file)
    motion = inventory[0][0][0].copy()
    motion.code = "BHZ"
    motion.response = Response.from_paz(
        zeros=[],
        poles=[],
        stage_gain=4,
        input_units="M/S",
        output_units="COUNTS",
    )
    inventory[0][0].channels.append(motion)
    inventory.write(row.inventory_file, format="STATIONXML")
    trace.stats.channel = "BHZ"
    trace.write(
        str(records[0].path.parent / "XX.TEST.00.BHZ.D.2020.001"), format="MSEED"
    )
    master = build_master_dataset(records[0].path.parent.parent / "inventories")
    score = "test_pressure_spectrum_vs_expected_noise_levels_score"
    master.loc[0, score] = 0.5
    exports = output.parent / "exports"
    save_master_dataset(master, exports)
    updates = run_processing(
        records[0].path.parent,
        exports,
        output,
        workers=2,
        ppsd_length=512,
        minimum_frequency=0.02,
        chunk_seconds=1024,
    )
    saved = load_master_dataset(exports)
    assert len(updates) == 2
    assert saved.ppsd_status.eq("complete").all()
    assert saved.ppsd_segments.eq("7").all()
    assert set(saved.ppsd_units) == {"Pa^2/Hz", "m^2/s^4/Hz"}
    assert saved.loc[0, score] == "0.5"


def test_seismic_velocity_response_yields_acceleration_density(dataset):
    row, records, trace, output = dataset
    inventory = read_inventory(row.inventory_file)
    inventory[0][0][0].response = Response.from_paz(
        zeros=[], poles=[], stage_gain=4, input_units="M/S", output_units="COUNTS"
    )
    inventory.write(row.inventory_file, format="STATIONXML")
    result = process_epoch(
        row,
        records,
        output,
        ppsd_length=512,
        minimum_frequency=0.02,
        chunk_seconds=1024,
    )
    corrected = PPSD.load_npz(result["ppsd_npz_file"])
    raw = PPSD(
        trace.stats,
        {"sensitivity": 1, "poles": [], "zeros": [], "gain": 1},
        special_handling="ringlaser",
        ppsd_length=512,
        period_limits=(0.25, 50),
        db_bins=(-200, 100, 1),
    )
    raw.add(trace)
    assert result["ppsd_units"] == "m^2/s^4/Hz"
    assert corrected.special_handling is None
    # Counts/(m/s) removal followed by one time derivative: |i*omega/4|^2,
    # averaged in dB over each period bin's FFT frequencies as ObsPy does.
    expected = []
    for left, right in zip(
        corrected.period_bin_left_edges, corrected.period_bin_right_edges, strict=True
    ):
        periods = corrected.psd_periods
        omega = 2 * np.pi / periods[(periods >= left) & (periods <= right)]
        expected.append(np.mean(20 * np.log10(omega / 4)))
    np.testing.assert_allclose(
        np.array(corrected.psd_values) - np.array(raw.psd_values),
        np.broadcast_to(expected, np.shape(raw.psd_values)),
        atol=1e-4,
    )


def test_identical_repeated_station_blocks_select_one_epoch(tmp_path):
    def station(gain):
        response = Response.from_paz(
            zeros=[], poles=[], stage_gain=gain, input_units="PA", output_units="COUNTS"
        )
        channel = Channel(
            code="CDH",
            location_code="00",
            latitude=1,
            longitude=2,
            elevation=-2000,
            depth=0,
            sample_rate=250,
            start_date=UTCDateTime("2017-06-15"),
            end_date=UTCDateTime("2018-02-16"),
            response=response,
        )
        return Station(
            code="A415A",
            latitude=1,
            longitude=2,
            elevation=-2000,
            site=Site(name="S"),
            channels=[channel],
        )

    row = {
        "Seed_id": "Z3.A415A.00.CDH",
        "channel_start_time": "2017-06-15T00:00:00.000000Z",
        "channel_end_time": "2018-02-16T00:00:00.000000Z",
    }
    path = tmp_path / "Z3.xml"
    Inventory([Network(code="Z3", stations=[station(4), station(4)])]).write(
        str(path), format="STATIONXML"
    )
    inventory, channel = processing.select_epoch(row | {"inventory_file": str(path)})
    assert channel.code == "CDH"
    assert sum(len(sta) for net in inventory for sta in net) == 1
    Inventory([Network(code="Z3", stations=[station(4), station(5)])]).write(
        str(path), format="STATIONXML"
    )
    with pytest.raises(ValueError, match="found 2"):
        processing.select_epoch(row | {"inventory_file": str(path)})
