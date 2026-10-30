"""Full-archive coverage planning reuses valid work and computes only gaps."""

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from obspy import Trace, UTCDateTime
from obspy.core.inventory import (
    Channel,
    Equipment,
    Inventory,
    Network,
    Site,
    Station,
)
from obspy.core.inventory.response import Response

from compy_pressure import coverage, coverage_report, coverage_run, network
from compy_pressure.archive import index_waveform_files
from compy_pressure.ppsd import process_epoch

START = UTCDateTime("2020-01-01")
END = START + 6 * 86400
DAY = 86400
PARAMS = {
    "ppsd_length": 512.0,
    "overlap": 0.5,
    "minimum_frequency": 0.02,
    "chunk_seconds": 4096.0,
}


def channel(code, units, *, location="00", gain=4.0, start=START, end=END, sensor=None):
    response = Response.from_paz(
        zeros=[], poles=[], stage_gain=gain, input_units=units, output_units="COUNTS"
    )
    response.instrument_sensitivity.frequency = 0.1
    return Channel(
        code=code,
        location_code=location,
        latitude=1,
        longitude=2,
        elevation=-1000,
        depth=0,
        sample_rate=1.0,
        start_date=start,
        end_date=end,
        response=response,
        sensor=Equipment(description=sensor) if sensor else None,
    )


def write_day(folder, seed, day, *, filename_location=None, layout="sds", zeros=False):
    net, sta, loc, cha = seed.split(".")
    start = START + day * DAY
    rng = np.random.default_rng(day + len(seed))
    data = (
        np.zeros(DAY, dtype=np.int32)
        if zeros
        else rng.integers(-1000, 1000, DAY, dtype=np.int32)
    )
    trace = Trace(data)
    trace.stats.update(
        dict(
            network=net,
            station=sta,
            location=loc,
            channel=cha,
            sampling_rate=1.0,
            starttime=start,
        )
    )
    name_loc = loc if filename_location is None else filename_location
    if layout == "sds":
        name = f"{net}.{sta}.{name_loc}.{cha}.D.{start.year}.{start.julday:03d}"
        path = folder / f"{cha}.D" / name
    else:
        name = f"{net}.{sta}.{name_loc or '--'}.{cha}.{start.strftime('%Y%m%d')}.mseed"
        path = folder / "waveforms" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    trace.write(str(path), format="MSEED")
    return path


def write_csv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(rows[0])
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def master_row(seed, inventory, **fields):
    return {
        "Seed_id": seed,
        "channel_start_time": str(START),
        "channel_end_time": str(END),
        "inventory_file": str(inventory),
        "sample_rate_hz": "1.0",
        "input_units": "PA" if seed.endswith("DH") else "M/S",
        "Lat": "1.0",
        "Lon": "2.0",
        "Profondeur": "1000.0",
    } | fields


def build_env(tmp_path, *, product_days=2, local_p_row=True, c_end=END):
    root = tmp_path / "source"
    data = root / "Data/XX"
    inventory = root / "INVENTORIES/XX.xml"
    inventory.parent.mkdir(parents=True)
    stations = {
        "A": [
            channel("BDH", "PA", sensor="Differential Pressure Gauge"),
            channel("BHZ", "M/S"),
            channel("LCE", "V"),
        ],
        "B": [channel("BDH", "PA", sensor="Differential Pressure Gauge")],
        "C": [channel("BDH", "PA", location="", end=c_end, sensor="Hydrophone")],
    }
    Inventory(
        [
            Network(
                code="XX",
                stations=[
                    Station(
                        code=code,
                        latitude=1,
                        longitude=2,
                        elevation=-1000,
                        site=Site(name=code),
                        channels=channels,
                    )
                    for code, channels in stations.items()
                ],
            )
        ],
        source="test",
    ).write(str(inventory), format="STATIONXML")
    for day in range(4):
        for seed in ("XX.A.00.BDH", "XX.A.00.BHZ", "XX.A.00.LCE", "XX.B.00.BDH"):
            write_day(data / seed.split(".")[1], seed, day)
        # Filename location 00 for a blank header location (authorized alias).
        write_day(data / "C", "XX.C..BDH", day, filename_location="00", zeros=day == 3)
    # A second copy of one day in another layout: excluded, never merged.
    write_day(data / "B", "XX.B.00.BDH", 1, layout="daily")
    # Historical product for A.BDH covering only the first days.
    product_root = root / "results/NETWORKS/HIST"
    row = master_row("XX.A.00.BDH", inventory)
    files = sorted((data / "A/BDH.D").iterdir())
    summary = process_epoch(
        row,
        index_waveform_files(files),
        product_root,
        end=str(START + product_days * DAY),
        **PARAMS,
    )
    product = summary["ppsd_npz_file"]
    registry = root / "exports/master_inventory.csv"
    registered = master_row(
        "XX.A.00.BDH",
        inventory,
        ppsd_status="complete",
        ppsd_npz_file=product,
        ppsd_units="Pa^2/Hz",
    )
    status, updates, _ = network.diagnostic_job(
        "pressure", registered, [], tmp_path / "no_master", {}
    )
    registered |= {
        "pressure_expected_status": status,
        "pressure_expected_result_file": updates[0]["pressure_expected_result_file"],
    }
    write_csv(registry, [registered, master_row("XX.A.00.BHZ", inventory)])
    docs = root / "docs"
    docs.mkdir()
    (docs / "rayleigh.json").write_text(json.dumps({"networks": []}))
    write_csv(
        docs / "microseism.csv",
        [],
        ["Seed_id", "channel_start_time", "channel_end_time", "inventory_file"],
    )
    local_p = master_row(
        "XX.A.00.BDH",
        inventory,
        master_epoch_match="exact",
        local_p_status="diagnostic_only",
        local_p_review_flags="legacy_noise_guard_unverified",
    )
    not_applicable = {k: "" for k in local_p} | {
        "Seed_id": "XX.C..BDH",
        "master_epoch_match": "no_master_epoch",
        "local_p_status": "not_applicable_sample_rate",
        "local_p_note": "usable band below local-P frequencies",
    }
    write_csv(
        docs / "local_p.csv",
        [local_p, not_applicable] if local_p_row else [],
        list(local_p) + ["local_p_note"],
    )
    config = {
        "schema": coverage.SCHEMA,
        "networks": ["XX"],
        "data_directory": "Data",
        "inventory_sources": {"XX": ["INVENTORIES/XX.xml"]},
        "layout_preference": ["sds", "sds_mseed", "daily"],
        "allowed_location_aliases": [["00", ""]],
        "gap_tolerance_seconds": 1.0,
        "ppsd_parameters": PARAMS,
        "product_registries": ["exports/master_inventory.csv"],
        "product_roots": ["results/NETWORKS/HIST"],
        "evidence": {
            "masters": ["exports/master_inventory.csv"],
            "rayleigh": "docs/rayleigh.json",
            "microseism": "docs/microseism.csv",
            "local_p": "docs/local_p.csv",
        },
    }
    path = docs / "coverage_sources.json"
    path.write_text(json.dumps(config))
    return SimpleNamespace(
        root=root,
        config=coverage.load_config(path, root),
        campaign=tmp_path / "campaign",
        product=Path(product),
        inventory=inventory,
        registry=registry,
    )


def plan_for(env, **kwargs):
    return coverage.build_plan(env.config, env.campaign, progress=False, **kwargs)


def channel_of(plan, seed):
    (found,) = [c for c in plan["channels"] if c["Seed_id"] == seed]
    return found


def tree_hashes(root):
    return {
        str(p): coverage.digest(p) for p in sorted(Path(root).rglob("*")) if p.is_file()
    }


@pytest.fixture
def env(tmp_path):
    return build_env(tmp_path)


def test_interval_arithmetic():
    assert coverage.merge_intervals([(0, 1), (1.5, 2), (2, 3)], 0.6) == [(0, 3)]
    assert coverage.subtract_intervals([(0, 10)], [(2, 3), (5, 20)]) == [
        (0, 2),
        (3, 5),
    ]
    assert coverage.clip_intervals([(0, 10)], 4, None) == [(4, 10)]


def test_filename_layouts_and_location_alias():
    assert coverage.filename_identity("YV.RR01.00.BDH.D.2013.026") == (
        "sds",
        "YV.RR01.00.BDH",
    )
    assert coverage.filename_identity("XO.LA33.00.HDH.D.2018.152.mseed") == (
        "sds_mseed",
        "XO.LA33.00.HDH",
    )
    assert coverage.filename_identity("YL.A02W.--.LDH.20100206.mseed") == (
        "daily",
        "YL.A02W..LDH",
    )
    assert coverage.filename_identity("inventory.xml") == (None, None)
    aliases = {("00", "")}
    assert coverage.identity_issue("XO.LA33.00.HDH", "XO.LA33..HDH", aliases) is None
    assert coverage.identity_issue("XO.LA33.00.RKO", "XO.LA33.AN.R", aliases)


def test_plan_reuses_saved_work_and_plans_only_missing_intervals(env):
    plan = plan_for(env)
    a = channel_of(plan, "XX.A.00.BDH")
    # The saved product covers two of four available days: supplement the rest.
    assert a["ppsd_decision"] == "supplement"
    supplements = [
        t for t in plan["tasks"] if t["kind"] == "ppsd" and t["seed_id"] == a["Seed_id"]
    ]
    assert len(supplements) == 1
    start, end = map(UTCDateTime, supplements[0]["interval"])
    assert abs(start - (START + 2 * DAY)) < PARAMS["ppsd_length"]
    assert abs(end - (START + 4 * DAY)) < 2
    merge = [t for t in plan["tasks"] if t["kind"] == "ppsd_merge"]
    assert merge[0]["dependencies"] == [supplements[0]["task_id"]]
    assert merge[0]["base_product_sha256"] == coverage.digest(env.product)
    # Stations without products are computed over their full available span.
    b = channel_of(plan, "XX.B.00.BDH")
    assert b["ppsd_decision"] == "compute"
    # Header location blank, filename 00: same stream, no mismatch issue.
    assert channel_of(plan, "XX.C..BDH")["ppsd_decision"] == "compute"
    # Non-pressure channels are not planned by default; non-physical are blocked.
    assert channel_of(plan, "XX.A.00.BHZ")["ppsd_decision"] == (
        "not_planned_out_of_scope"
    )
    assert channel_of(plan, "XX.A.00.LCE")["ppsd_decision"] == "blocked"
    assert plan["archive"]["excluded_duplicate_copies"][0]["layout"] == "daily"
    for task in plan["tasks"]:
        for field in (
            "task_id",
            "key",
            "prerequisites",
            "intended_coverage",
            "parameters",
            "output_path",
            "reuse_decision",
            "dependencies",
            "reason",
        ):
            assert field in task, (task["kind"], field)
        assert Path(task["output_path"]).is_relative_to(env.campaign)


def test_all_channel_scope_is_explicit(env):
    plan = plan_for(env, scope="all-channels")
    z = channel_of(plan, "XX.A.00.BHZ")
    assert z["ppsd_decision"] == "compute"
    assert any(t["seed_id"] == "XX.A.00.BHZ" for t in plan["tasks"])
    assert channel_of(plan, "XX.A.00.LCE")["ppsd_decision"] == "blocked"


def test_full_product_is_reused_and_its_saved_result_kept(tmp_path):
    env = build_env(tmp_path, product_days=6)
    plan = plan_for(env)
    a = channel_of(plan, "XX.A.00.BDH")
    assert a["ppsd_decision"] == "reuse_existing"
    assert not [
        t
        for t in plan["tasks"]
        if t["seed_id"] == "XX.A.00.BDH" and t["kind"].startswith("ppsd")
    ]
    expected = a["checks"]["expected_pressure"]
    assert expected["decision"] == "reuse_existing"
    assert expected["result"]["binding"] == "registry_row_and_npz_hash"
    # Checks without saved results are planned against a view, not the source.
    infragravity = [
        t
        for t in plan["tasks"]
        if t["kind"] == "infragravity" and t["seed_id"] == "XX.A.00.BDH"
    ]
    assert infragravity[0]["output_path"].endswith("views")


def test_response_or_parameter_difference_requires_review(tmp_path):
    env = build_env(tmp_path, product_days=6)
    snapshot = env.product.with_name("inventory.xml")
    text = snapshot.read_text()
    snapshot.write_text(text.replace("<Value>4", "<Value>5", 1))
    plan = plan_for(env)
    a = channel_of(plan, "XX.A.00.BDH")
    assert a["ppsd_decision"] == "review_required"
    assert "response_differs" in a["ppsd_reason"]
    assert not [
        t
        for t in plan["tasks"]
        if t["seed_id"] == "XX.A.00.BDH" and t["kind"].startswith("ppsd")
    ]
    snapshot.write_text(text)
    manifest = env.product.with_name("processing.json")
    data = json.loads(manifest.read_text())
    data["config"]["ppsd_length_seconds"] = 1024.0
    manifest.write_text(json.dumps(data))
    plan = plan_for(env)
    assert "parameters_differ" in channel_of(plan, "XX.A.00.BDH")["ppsd_reason"]


def test_epoch_identity_splits_coverage_and_reports_unmatched_data(tmp_path):
    # The C epoch ends after three days: its last day is outside every epoch.
    env = build_env(tmp_path, c_end=START + 3 * DAY)
    plan = plan_for(env)
    c = channel_of(plan, "XX.C..BDH")
    assert c["available_end"].startswith("2020-01-04T00:00:00")
    (unmatched,) = plan["archive"]["unmatched_streams"]
    assert unmatched["seed_id"] == "XX.C..BDH"
    assert unmatched["reason"] == "data_outside_inventory_epochs"


def test_saved_product_from_different_epoch_requires_review(tmp_path):
    from obspy import read_inventory

    env = build_env(tmp_path, product_days=6)
    inventory = read_inventory(str(env.inventory))
    inventory[0][0][0].end_date = START + 3 * DAY
    inventory.write(str(env.inventory), format="STATIONXML")
    a = channel_of(plan_for(env), "XX.A.00.BDH")
    assert a["ppsd_decision"] == "review_required"
    assert "epoch_bounds_differ" in a["ppsd_reason"]
    assert a["ppsd_tasks"] == []


@pytest.mark.parametrize("status", ["no_data", "insufficient_data"])
def test_finished_empty_task_is_not_reported_as_covered(env, status):
    plan = plan_for(env)
    channel = channel_of(plan, "XX.B.00.BDH")
    task = next(t for t in plan["tasks"] if t["task_id"] == channel["ppsd_target"])
    coverage_run.write_state(
        env.campaign, task, status="complete", result_status=status, outputs=[]
    )
    rows, _ = coverage_report.coverage_tables(plan, env.campaign, [])
    row = next(r for r in rows if r["Seed_id"] == channel["Seed_id"])
    assert row["product_covered_days"] == 0
    assert row["available_days"] > 0


def test_calibration_evidence_is_reused_or_deferred_never_run(env):
    plan = plan_for(env)
    checks = channel_of(plan, "XX.A.00.BDH")["checks"]
    assert checks["local_p"]["decision"] == "reuse_historical"
    assert "legacy" in checks["local_p"]["reason"]
    assert checks["rayleigh"]["decision"] == "deferred"
    assert channel_of(plan, "XX.B.00.BDH")["checks"]["local_p"]["decision"] == (
        "deferred"
    )
    # An unjoined historical row saying "not applicable" is not a missing task.
    c_local_p = channel_of(plan, "XX.C..BDH")["checks"]["local_p"]
    assert c_local_p["decision"] == "not_applicable"
    assert "not_applicable_sample_rate" in c_local_p["reason"]
    assert not {t["kind"] for t in plan["tasks"]} & set(coverage.CALIBRATION_CHECKS)


def test_missing_configured_sources_are_refused(env):
    config = json.loads((env.root / "docs/coverage_sources.json").read_text())
    config["product_roots"].append("results/NETWORKS/ABSENT")
    path = env.root / "docs/broken.json"
    path.write_text(json.dumps(config))
    with pytest.raises(FileNotFoundError, match="ABSENT"):
        coverage.load_config(path, env.root)
    with pytest.raises(ValueError, match="protected"):
        coverage.guard_campaign_dir(env.root / "Data/campaign", env.config)


def test_execute_resume_and_no_duplicate_work(env, monkeypatch):
    protected = tree_hashes(env.root)
    dry = coverage_run.run_campaign(env.config, env.campaign, workers=1)
    assert dry["executed"] is False and dry["pending"] > 0
    assert not (env.campaign / "tasks").exists()

    calls = []
    original = network.ppsd_job

    def interrupt_once(row, *args):
        calls.append(row["Seed_id"])
        if len(calls) == 1:
            raise KeyboardInterrupt
        return original(row, *args)

    monkeypatch.setattr(network, "ppsd_job", interrupt_once)
    with pytest.raises(KeyboardInterrupt):
        coverage_run.run_campaign(env.config, env.campaign, execute=True, workers=1)
    states = [
        json.loads(p.read_text())["status"]
        for p in (env.campaign / "tasks").rglob("state.json")
    ]
    assert "running" in states or "interrupted" in states
    # The lock is released after an interruption; the next run resumes.
    result = coverage_run.run_campaign(
        env.config, env.campaign, execute=True, workers=1
    )
    assert result["ok"], result
    assert tree_hashes(env.root) == protected
    run = json.loads((Path(result["run_dir"]) / "run.json").read_text())
    assert run["protected_inputs_unchanged"] is True
    # Merged A product keeps the base windows and adds only new ones.
    merges = list((env.campaign / "products").rglob("merge.json"))
    merge = json.loads(merges[0].read_text())
    assert merge["supplements"][0]["added_windows"] > 0
    assert merge["total_windows"] == merge["base_windows"] + sum(
        s["added_windows"] for s in merge["supplements"]
    )
    # Checks on a reused saved product write into the campaign, never the source.
    assert not (env.product.parent / "infragravity").exists()

    calls.clear()
    monkeypatch.setattr(
        network,
        "diagnostic_job",
        lambda *a, **k: pytest.fail("completed check must not rerun"),
    )
    again = coverage_run.run_campaign(env.config, env.campaign, execute=True, workers=1)
    assert again["ok"] and calls == []
    assert set(again["counts"]) == {"already_complete"}
    assert tree_hashes(env.root) == protected


def test_failure_isolation_and_dependency_blocking(env, monkeypatch):
    original = network.ppsd_job

    def fail_b(row, *args):
        if row["Seed_id"] == "XX.B.00.BDH":
            raise RuntimeError("synthetic read failure")
        return original(row, *args)

    monkeypatch.setattr(network, "ppsd_job", fail_b)
    result = coverage_run.run_campaign(
        env.config, env.campaign, execute=True, workers=1
    )
    assert not result["ok"]
    run = json.loads((Path(result["run_dir"]) / "run.json").read_text())
    plan = json.loads((Path(result["run_dir"]) / "plan/plan.json").read_text())
    kinds = {t["task_id"]: (t["kind"], t["seed_id"]) for t in plan["tasks"]}
    outcome = {kinds[k]: v for k, v in run["tasks"].items()}
    assert outcome[("ppsd", "XX.B.00.BDH")] == "failed"
    assert outcome[("expected_pressure", "XX.B.00.BDH")] == "blocked_by_dependency"
    assert outcome[("ppsd", "XX.C..BDH")] == "complete"
    assert outcome[("expected_pressure", "XX.C..BDH")] == "complete"
    # The reference pool is a soft dependency: the failure is recorded, not hidden.
    assert outcome[("nearby", "XX.C..BDH")] == "complete"
    nearby_c = [k for k, v in kinds.items() if v == ("nearby", "XX.C..BDH")][0]
    state = json.loads((env.campaign / "tasks" / nearby_c / "state.json").read_text())
    assert state["pool_excluded"] == ["XX.B.00.BDH"]
    log = next((env.campaign / "tasks").glob("ppsd-*/log.txt"), None)
    assert log is not None
    # Once the failed PPSD succeeds, nearby results with the old pool are stale.
    monkeypatch.setattr(network, "ppsd_job", original)
    coverage_run.run_campaign(
        env.config, env.campaign, execute=True, workers=1, only_kinds=["ppsd"]
    )
    report = coverage_report.write_report(plan_for(env), env.campaign)
    rows = list(csv.DictReader((report / "tasks_status.csv").open()))
    states = {(r["kind"], r["seed_id"]): r["state"] for r in rows}
    assert states[("nearby", "XX.C..BDH")] != "complete"
    assert states[("expected_pressure", "XX.C..BDH")] == "complete"
    rerun = coverage_run.run_campaign(env.config, env.campaign, execute=True, workers=1)
    assert rerun["ok"], rerun
    state = json.loads((env.campaign / "tasks" / nearby_c / "state.json").read_text())
    assert state["pool_excluded"] == [] and state["attempts"] == 2
    run = json.loads((Path(rerun["run_dir"]) / "run.json").read_text())
    assert run["tasks"][nearby_c] == "complete"
    expected_c = [
        k for k, v in kinds.items() if v == ("expected_pressure", "XX.C..BDH")
    ][0]
    assert run["tasks"][expected_c] == "already_complete"


def test_checks_on_reused_product_write_only_campaign_views(tmp_path):
    env = build_env(tmp_path, product_days=6)
    protected = tree_hashes(env.root)
    result = coverage_run.run_campaign(
        env.config,
        env.campaign,
        execute=True,
        workers=1,
        only_kinds=["infragravity"],
        only_seeds=["XX.A.00.BDH"],
    )
    assert result["ok"], result
    assert tree_hashes(env.root) == protected
    views = list((env.campaign / "views").rglob("infragravity/comparison.json"))
    assert len(views) == 1
    view = json.loads((views[0].parents[1] / "view.json").read_text())
    assert view["sources"]["ppsd.npz"]["sha256"] == coverage.digest(env.product)


def test_campaign_lock_refuses_concurrent_execution(env):
    with coverage_run.campaign_lock(env.campaign):
        with pytest.raises(coverage_run.CampaignLocked):
            coverage_run.run_campaign(env.config, env.campaign, execute=True, workers=1)


def test_report_heatmaps_keep_floors_and_separate_weightings(env):
    result = coverage_run.run_campaign(
        env.config, env.campaign, execute=True, workers=1, only_kinds=["ppsd"]
    )
    assert result["executed"]
    plan = plan_for(env)
    output = coverage_report.write_report(plan, env.campaign)
    manifest = json.loads((output / "heatmaps/heatmap_manifest.json").read_text())
    # One channel code, two declared sensor families: never pooled.
    groups = {g["sensor_family"]: g for g in manifest["groups"]}
    assert set(groups) == {"dpg", "hydrophone"}
    assert groups["dpg"]["stations"] == 2 and groups["hydrophone"]["stations"] == 1
    assert groups["hydrophone"]["declared_sensor_labels"] == {"Hydrophone": 1}
    # Station C's all-zero day gives numerical floors; they are counted, kept.
    group = groups["hydrophone"]
    assert group["per_station"]["C"]["zero_power_floor_values"] > 0
    assert group["db_range_retained"][0] < -3000
    for family, group in groups.items():
        data = np.load(output / f"heatmaps/XX_BDH_{family}_heatmap.npz")
        balanced = data["station_balanced_probability"]
        covered = data["stations_per_frequency"] > 0
        assert np.allclose(balanced[:, covered].sum(axis=0), 1.0)
        assert data["pooled_counts"].sum(axis=0).max() == group["windows"]
    assert manifest["summary_curves"].startswith("none")
    rows = list(csv.DictReader((output / "coverage_channels.csv").open()))
    assert {r["Seed_id"] for r in rows} >= {"XX.A.00.BDH", "XX.B.00.BDH"}


@pytest.mark.parametrize(
    ("label", "family"),
    [
        ("Parosci APG Pressure/LDEO_APG_logger", "apg"),
        ("SIO Differential Pressure Gauge/Quanterra Q330 Lin", "dpg"),
        ("DPG5004 with logger 21 (card #?),Cox-Webb Differen", "dpg"),
        ("DPG", "dpg"),
        ("HTI-04-PCA/ULF,Hydrophone", "hydrophone"),
        ("LDEO HTI90U/LDEO OBS 2011 logger", "hydrophone"),
        ("", "unclassified"),
    ],
)
def test_declared_sensor_families_of_real_labels(label, family):
    assert coverage_report.sensor_family(label) == family


def test_parallel_workers_run_ppsd_tasks_in_separate_processes(env):
    protected = tree_hashes(env.root)
    result = coverage_run.run_campaign(
        env.config, env.campaign, execute=True, workers=2, only_kinds=["ppsd"]
    )
    assert result["ok"], result
    run = json.loads((Path(result["run_dir"]) / "run.json").read_text())
    plan = json.loads((Path(result["run_dir"]) / "plan/plan.json").read_text())
    ppsd = [t["task_id"] for t in plan["tasks"] if t["kind"] == "ppsd"]
    assert len(ppsd) == 3 and all(run["tasks"][t] == "complete" for t in ppsd)
    for task_id in ppsd:
        log = (env.campaign / "tasks" / task_id / "log.txt").read_text()
        assert "PPSD" in log and "segments saved" in log
    assert tree_hashes(env.root) == protected
