"""Verify paired days, distinct stations and known spectral offsets."""

import csv
import hashlib
import json
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pytest
from obspy import Trace, UTCDateTime
from obspy.core.inventory import Channel, Equipment, Inventory, Network, Site, Station
from obspy.core.inventory.response import Response
from obspy.signal import PPSD
from openpyxl import load_workbook

from compy_pressure import nearby_test as nearby
from compy_pressure import network
from compy_pressure.inventory import (
    RESULT_COLUMNS,
    TEST_COLUMNS,
    build_master_dataset,
    load_master_dataset,
    save_master_dataset,
    update_nearby_test_results,
)


def row(station="A", net="NN", location="00", lon=0, depth=1000, sensor="Hydrophone"):
    return {
        "Seed_id": f"{net}.{station}.{location}.BDH",
        "Lat": "0",
        "Lon": str(lon),
        "Profondeur": str(depth),
        "sensor_type": sensor,
        "sensor_model": "",
        "channel_start_time": "2020-01-01T00:00:00.000000Z",
        "channel_end_time": "",
        "inventory_file": "fixture.xml",
    }


def daily(record, days=None, frequencies=None, offset=0):
    frequencies = np.array(
        frequencies if frequencies is not None else [0.01, 0.02, 0.04, 0.08]
    )
    days = np.array(
        days if days is not None else [f"2020-01-{i:02}" for i in range(1, 31)]
    )
    values = np.tile(3 * np.log(frequencies) + offset, (len(days), 1))
    return nearby.DailyPressure(
        record,
        frequencies,
        frequencies / np.sqrt(2),
        frequencies * np.sqrt(2),
        days,
        values,
        {"minimum_hz": 0.001, "maximum_hz": 1},
        (10800, 0.5, 1),
        {},
        [14] * len(days),
        [22] * len(days),
    )


def test_geodesic_depth_sensor_tiers_and_self_station():
    criteria = nearby.Criteria()
    target = row()
    preferred = nearby.candidate(target, row("B", lon=0.1, depth=1020), criteria)
    assert preferred["tier"] == "preferred"
    assert preferred["distance_km"] == pytest.approx(11.13195, abs=1e-4)
    assert preferred["depth_difference_fraction"] == pytest.approx(20 / 1010)
    assert (
        nearby.candidate(target, row("B", lon=0.6), criteria)["tier"] == "exploratory"
    )
    assert (
        nearby.candidate(target, row("B", lon=2), criteria)["reason"]
        == "outside_maximum_distance"
    )
    assert (
        nearby.candidate(target, row("B", depth=1300), criteria)["reason"]
        == "outside_maximum_depth_difference"
    )
    assert (
        nearby.candidate(target, row("A", location="10"), criteria)["reason"]
        == "same_station"
    )
    assert (
        nearby.candidate(target, row("B", sensor=""), criteria)["tier"] == "exploratory"
    )
    failed = row("B") | {TEST_COLUMNS[0]: "0"}
    assert (
        nearby.candidate(target, failed, criteria)["reason"]
        == "reference_has_recorded_failed_check"
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"maximum_distance_km": np.nan},
        {"preferred_distance_km": 200},
        {"preferred_depth_fraction": 0.2},
        {"minimum_common_days": 0},
        {"minimum_reference_stations": 1},
        {"minimum_daily_coverage_hours": 25},
    ],
)
def test_invalid_criteria(kwargs):
    with pytest.raises(ValueError):
        nearby.Criteria(**kwargs)


def test_daily_coverage_is_union_and_midnight_windows_are_excluded():
    start = UTCDateTime("2020-01-01")
    times = [start + h * 3600 for h in [0, 1.5, 3, 6, 22]]
    values = np.array([[i, i + 1] for i in range(5)])
    days, spectra, counts, coverage = nearby.daily_spectra(times, values, 3 * 3600, 8)
    assert list(days) == ["2020-01-01"]
    assert counts == [4]
    assert coverage == [9]
    np.testing.assert_array_equal(spectra[0], [1.5, 2.5])
    assert nearby.daily_spectra(times, values, 3 * 3600, 10)[0].size == 0


def test_log_frequency_interpolation_and_no_extrapolation():
    target = daily(row(), offset=6)
    peer = daily(row("B"), frequencies=[0.012, 0.024, 0.048, 0.096])
    aligned, reason = nearby.align_pair(target, peer, nearby.Criteria())
    assert not reason
    assert aligned["frequency_indices"].tolist() == [1, 2, 3]
    np.testing.assert_allclose(aligned["difference"], 6, atol=1e-12)
    too_short = daily(row("B"), days=["2020-01-01"])
    assert (
        nearby.align_pair(target, too_short, nearby.Criteria())[1]
        == "insufficient_common_days"
    )
    incompatible = daily(row("B"))
    incompatible.processing = (3600, 0.5, 1)
    assert (
        nearby.align_pair(target, incompatible, nearby.Criteria())[1]
        == "incompatible_processing_settings"
    )


def test_reference_group_does_not_count_two_locations_as_two_stations():
    target = daily(row())
    pairs = []
    for record in [row("B"), row("B", location="10")]:
        peer = daily(record)
        aligned, _ = nearby.align_pair(target, peer, nearby.Criteria())
        audit = nearby.candidate(target.row, record, nearby.Criteria())
        pairs.append({"data": peer, "aligned": aligned, "audit": audit})
    group = nearby.reference_group(target, pairs, nearby.Criteria())
    assert len(group["pairs"]) == 1
    assert not group["frequency_mask"].any()
    peer = daily(row("C", net="ZZ"))
    aligned, _ = nearby.align_pair(target, peer, nearby.Criteria())
    pairs.append(
        {
            "data": peer,
            "aligned": aligned,
            "audit": nearby.candidate(target.row, peer.row, nearby.Criteria()),
        }
    )
    group = nearby.reference_group(target, pairs, nearby.Criteria())
    assert len(group["pairs"]) == 2
    assert group["frequency_mask"].all()


def test_joint_coverage_requires_two_references_on_the_same_days():
    days = [str((UTCDateTime("2020-01-01") + i * 86400).date) for i in range(60)]
    target = daily(row(), days=days)
    pairs = []
    for record, selected in [(row("B"), days[:30]), (row("C"), days[30:])]:
        peer = daily(record, days=selected)
        aligned, _ = nearby.align_pair(target, peer, nearby.Criteria())
        pairs.append(
            {
                "data": peer,
                "aligned": aligned,
                "audit": nearby.candidate(target.row, peer.row, nearby.Criteria()),
            }
        )
    assert not nearby.reference_group(target, pairs, nearby.Criteria())[
        "frequency_mask"
    ].any()


@pytest.fixture
def saved_neighbors(tmp_path):
    start = UTCDateTime("2020-01-01")
    inventory_dir = tmp_path / "inventories"
    inventory_dir.mkdir()
    stations, networks = [], []
    specs = [
        ("NN", "A", 0, 1000, 2),
        ("NN", "B", 0.1, 1010, 1),
        ("ZZ", "C", 0.2, 990, 1),
    ]
    for net, name, lon, depth, amplitude in specs:
        response = Response.from_paz(
            zeros=[], poles=[], stage_gain=4, input_units="PA", output_units="COUNTS"
        )
        response.instrument_sensitivity.frequency = 1
        channel = Channel(
            code="BDH",
            location_code="00",
            latitude=0,
            longitude=lon,
            elevation=-depth,
            depth=0,
            sample_rate=8,
            start_date=start,
            end_date=start + 31 * 86400,
            response=response,
            sensor=Equipment(type="Hydrophone"),
        )
        sta = Station(
            code=name,
            latitude=0,
            longitude=lon,
            elevation=-depth,
            site=Site(name=name),
            channels=[channel],
        )
        networks.append(Network(code=net, stations=[sta]))
        stations.append((net, sta, amplitude))
    inv = Inventory(networks, source="Test")
    inv.write(str(inventory_dir / "all.xml"), format="STATIONXML")
    master = build_master_dataset(inventory_dir).astype(object).fillna("")
    base = np.random.default_rng(5).integers(-1000, 1000, 4096, dtype=np.int32)
    for index, record in master.iterrows():
        net, sta, _, _ = record.Seed_id.split(".")
        amplitude = next(a for n, s, a in stations if n == net and s.code == sta)
        trace = Trace(base * amplitude)
        trace.stats.update(
            dict(
                network=net,
                station=sta,
                location="00",
                channel="BDH",
                sampling_rate=8,
                starttime=start,
            )
        )
        ppsd = PPSD(
            trace.stats,
            inv,
            special_handling="hydrophone",
            ppsd_length=512,
            period_limits=(0.25, 50),
        )
        for day in range(30):
            for hour in [2, 10, 18]:
                trace.stats.starttime = start + day * 86400 + hour * 3600
                ppsd.add(trace)
        directory = tmp_path / "results" / net / sta / f"{record.Seed_id}__epoch"
        directory.mkdir(parents=True)
        ppsd.save_npz(str(directory / "ppsd.npz"))
        inv.select(network=net, station=sta).write(
            str(directory / "inventory.xml"), format="STATIONXML"
        )
        master.loc[index, "ppsd_units"] = "Pa^2/Hz"
        master.loc[index, "ppsd_status"] = "complete"
        master.loc[index, "ppsd_npz_file"] = str(directory / "ppsd.npz")
    master.loc[2, TEST_COLUMNS[1]] = 0.5
    master.loc[2, nearby.STATUS_COLUMN] = "manual_review"
    exports = tmp_path / "exports"
    save_master_dataset(master, exports)
    book = load_workbook(exports / "master_inventory.xlsx")
    book.create_sheet("Manual notes")["A1"] = "preserve me"
    book.save(exports / "master_inventory.xlsx")
    book.close()
    return exports


def test_real_ppsd_offset_recovery_outputs_and_master_preservation(saved_neighbors):
    exports = saved_neighbors
    before = load_master_dataset(exports)
    hashes = {
        p: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in map(Path, before.ppsd_npz_file)
    }
    results = nearby.run_nearby_test(
        exports,
        criteria=replace(nearby.Criteria(), minimum_daily_coverage_hours=0.25),
        show=False,
    )
    target = next(result for result in results if result["seed_id"] == "NN.A.00.BDH")
    assert target["status"] == "comparison_only"
    assert target["reference_stations"] == 2
    assert target["joint_common_days"] == 30
    assert target["median_difference_db"] == pytest.approx(20 * np.log10(2), abs=2e-5)
    assert target["score"] is None
    assert all(Path(p).is_file() for p in target["files"].values())
    after = load_master_dataset(exports)
    assert after[RESULT_COLUMNS].equals(before[RESULT_COLUMNS])
    assert after.iloc[2][nearby.STATUS_COLUMN] == "manual_review"
    preserved = [c for c in before if c != nearby.STATUS_COLUMN]
    assert after[preserved].equals(before[preserved])
    for path, digest in hashes.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    book = load_workbook(exports / "master_inventory.xlsx")
    assert book["Manual notes"]["A1"].value == "preserve me"
    headers = [c.value for c in book["Master inventory"][1]]
    assert (
        book["Master inventory"]
        .cell(2, headers.index("nearby_reference_stations") + 1)
        .value
        == 2
    )
    book.close()


def test_no_peer_and_exploratory_results_are_unresolved_without_scores(saved_neighbors):
    master = load_master_dataset(saved_neighbors)
    criteria = replace(nearby.Criteria(), minimum_daily_coverage_hours=0.25)
    target = master.iloc[[0]]
    no_peers = nearby.run_nearby_test(
        saved_neighbors,
        selected_rows=target,
        reference_rows=master.iloc[:0],
        criteria=criteria,
        update_master=False,
        show=False,
    )[0]
    assert no_peers["status"] == "no_suitable_neighbors"
    assert no_peers["median_difference_db"] is None
    assert Path(no_peers["files"]["candidates"]).is_file()
    exploratory = nearby.run_nearby_test(
        saved_neighbors,
        selected_rows=target,
        reference_rows=master,
        criteria=replace(criteria, preferred_distance_km=1),
        update_master=False,
        show=False,
    )[0]
    assert exploratory["status"] == "exploratory_only"
    assert exploratory["score"] is None


def test_reference_change_invalidates_network_cache_and_network_view(
    saved_neighbors, tmp_path
):
    exports = saved_neighbors
    master = load_master_dataset(exports)
    files = {seed: [] for seed in master.Seed_id}
    options = {
        "criteria": asdict(
            replace(nearby.Criteria(), minimum_daily_coverage_hours=0.25)
        )
    }
    root = tmp_path / "workflow"
    report = network.run_stage("nearby", exports, root, files, options=options)
    assert report["status_counts"] == {"comparison_only": 3}
    before = next(
        e["fingerprint"]
        for e in report["entries"].values()
        if e["seed_id"] == "NN.A.00.BDH"
    )
    changed = load_master_dataset(exports)
    changed.loc[1, "Profondeur"] = "1400"
    changed.to_csv(exports / "master_inventory.csv", index=False)
    second = network.run_stage("nearby", exports, root, files, options=options)
    entry = next(e for e in second["entries"].values() if e["seed_id"] == "NN.A.00.BDH")
    assert entry["fingerprint"] != before
    assert entry["status"] == "insufficient_reference_stations"
    assert (
        len(network.view_plots(exports, {"NN.A.00.BDH": []}, "nearby", show=False)) == 3
    )


def test_rms_is_reported_with_medians(saved_neighbors, tmp_path):
    assert nearby.finite_rms([3.0, -4.0, np.nan]) == pytest.approx(np.sqrt(12.5))
    assert nearby.finite_rms([np.nan]) is None
    exports = saved_neighbors
    files = {seed: [] for seed in load_master_dataset(exports).Seed_id}
    options = {
        "criteria": asdict(
            replace(nearby.Criteria(), minimum_daily_coverage_hours=0.25)
        )
    }
    report = network.run_stage("nearby", exports, tmp_path, files, options=options)
    tables = {
        Path(p).name: Path(p)
        for entry in report["entries"].values()
        for p in entry["files"]
        if Path(p).suffix == ".csv"
    }
    with tables["frequency_comparison.csv"].open() as stream:
        rows = list(csv.DictReader(stream))
    assert rows
    for row in rows:
        # RMS includes the median offset, so it can never be smaller.
        assert (
            float(row["difference_rms_db"])
            >= abs(float(row["difference_median_db"])) - 1e-9
        )
    with tables["daily_comparison.csv"].open() as stream:
        assert "rms_difference_db" in next(csv.reader(stream))


def test_guard_and_update_scope(saved_neighbors, tmp_path):
    report = tmp_path / "_workflow/ppsd.json"
    report.parent.mkdir()
    report.write_text(json.dumps({"selected_epoch_keys": ["a", "b"], "entries": {}}))
    with pytest.raises(ValueError, match="PPSD batch is unfinished"):
        nearby.ensure_master_idle(tmp_path)
    report.write_text(
        json.dumps(
            {"selected_epoch_keys": ["a", "b"], "status_counts": {"complete": 2}}
        )
    )
    nearby.ensure_master_idle(tmp_path)
    master = load_master_dataset(saved_neighbors)
    with pytest.raises(ValueError, match="nearby_"):
        update_nearby_test_results(
            saved_neighbors,
            [master.iloc[0][nearby.KEYS].to_dict() | {TEST_COLUMNS[1]: 1.0}],
        )
