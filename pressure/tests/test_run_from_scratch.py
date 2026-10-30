"""From-scratch pipeline: configuration, workspaces, commands and records."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from obspy import UTCDateTime
from obspy.core.inventory import Channel, Inventory, Network, Site, Station
from obspy.core.inventory.response import Response

from scripts import run_from_scratch as rfs

CONFIG = Path(__file__).resolve().parents[1] / "docs/from_scratch_config.json"
START = UTCDateTime("2012-10-01")
END = UTCDateTime("2013-10-01")


def channel(code, units, response=True):
    kwargs = {}
    if response:
        kwargs["response"] = Response.from_paz(
            zeros=[0j],
            poles=[-0.1 + 0j],
            stage_gain=2000.0,
            input_units=units,
            output_units="COUNTS",
        )
    return Channel(
        code=code,
        location_code="00",
        latitude=-21.0,
        longitude=55.0,
        elevation=-2918,
        depth=0,
        sample_rate=50,
        start_date=START,
        end_date=END,
        dip=-90 if code.endswith("Z") else 0,
        azimuth=0,
        **kwargs,
    )


def write_inventory(path):
    channels = [channel("BDH", "PA"), channel("BHZ", "M/S"), channel("LOG", "", False)]
    station = Station("RR99", -21.0, 55.0, -2918, channels=channels, site=Site("RR99"))
    Inventory([Network("YV", stations=[station])]).write(str(path), "STATIONXML")
    return path


def test_config_loads_and_rejects_bad_settings(tmp_path):
    config = rfs.load_config(CONFIG)
    assert set(config["networks"]) == {"YV", "YL", "7D", "XO"}
    assert config["identity_control"]["g"] == "0.6787278140712512"
    assert config["networks"]["YV"]["pressure_sign"] == -1
    bad = json.loads(CONFIG.read_text()) | {"event_policies": ["sometimes"]}
    (tmp_path / "bad.json").write_text(json.dumps(bad))
    with pytest.raises(ValueError, match="event_policies"):
        rfs.load_config(tmp_path / "bad.json")
    (tmp_path / "old.json").write_text(json.dumps({"schema": "old"}))
    with pytest.raises(ValueError, match="schema"):
        rfs.load_config(tmp_path / "old.json")


def test_required_inputs_follow_the_selected_networks():
    config = rfs.load_config(CONFIG)
    yv = {str(p) for p in rfs.required_inputs(config, ["YV"])}
    assert any(p.endswith("INVENTORIES/YV.xml") for p in yv)
    assert not any("7D_ALL_NO_TILT" in p for p in yv)
    every = {str(p) for p in rfs.required_inputs(config, ["7D", "XO"])}
    assert any(
        p.endswith("7D_ALL_NO_TILT/_inputs/catalog_usgs_m5p5.json") for p in every
    )
    assert any(p.endswith("INVENTORIES/XO.xml") for p in every)


def test_snapshots_hold_one_exact_epoch_and_skip_unusable(tmp_path):
    from compy_pressure.ppsd import select_epoch

    xml = write_inventory(tmp_path / "YV.xml")
    master = rfs.inventory_master([xml], tmp_path / "links")
    assert len(master) == 3
    assert set(master.inventory_file) == {str(xml.resolve())}
    updated, counts = rfs.write_snapshots(master, tmp_path / "workspace")
    assert counts == {"written": 2, "unusable:ValueError": 1}
    row = updated[updated.Seed_id == "YV.RR99.00.BDH"].iloc[0]
    snapshot = Path(row[rfs.SNAPSHOT_COLUMN])
    from obspy import read_inventory

    saved = read_inventory(str(snapshot))
    channels = [c for n in saved for s in n for c in s]
    _, expected = select_epoch(row)
    assert len(channels) == 1 and channels[0].code == "BDH"
    assert channels[0].start_date == expected.start_date
    assert (
        channels[0].response.instrument_sensitivity.value
        == expected.response.instrument_sensitivity.value
    )
    log = updated[updated.Seed_id == "YV.RR99.00.LOG"].iloc[0]
    assert log[rfs.SNAPSHOT_COLUMN] == ""
    # An existing snapshot (e.g. from the XO preparation) is kept.
    again, counts = rfs.write_snapshots(updated, tmp_path / "workspace")
    assert counts["kept"] == 2


def test_rayleigh_commands_reproduce_each_frozen_setting(tmp_path):
    config = rfs.load_config(CONFIG)

    def args(name, policy, stations=()):
        command = [
            str(c)
            for c in rfs.rayleigh_command(
                config, name, policy, tmp_path / name, stations, "python"
            )
        ]
        return command, dict(zip(command, command[1:], strict=False))

    yv, options = args("YV", "bounded", ["RR52"])
    assert options["--rayleigh-pressure-sign"] == "-1"
    assert options["--rayleigh-event-sampling"] == "bounded"
    assert options["--rayleigh-tilt-method"] == "none"
    assert options["--master-dir"].endswith("YV/_calibration_master")
    assert options["--station"] == "RR52" and "--rayleigh-no-master-update" in yv
    assert "--rayleigh-event-resample-rate" not in yv
    xo, options = args("XO", "all")
    assert options["--rayleigh-event-resample-rate"] == "1.0"
    assert options["--rayleigh-native-rate-tolerance-ppm"] == "0.1"
    assert options["--data-dir"].endswith("XO/_inputs/waveforms")
    assert options["--rayleigh-catalog"].endswith("XO/_inputs/catalog_usgs_m5p5.json")
    assert options["--rayleigh-pressure-sign"] == "1"


def test_coverage_config_disables_every_kind_of_reuse(tmp_path):
    config = rfs.load_config(CONFIG)
    path = rfs.coverage_config(
        rfs.resolve(config["coverage_template"]), ["YV", "XO"], tmp_path / "c.json"
    )
    written = json.loads(path.read_text())
    assert written["product_registries"] == []
    assert written["product_roots"] == []
    assert written["evidence"] == {}
    assert written["networks"] == ["YV", "XO"]
    assert written["inventory_sources"]["YV"]


def workflow(path, entries):
    path.write_text(
        json.dumps(
            {
                "entries": {
                    seed: {
                        "seed_id": seed,
                        "status": status,
                        "update": {"rayleigh_pressure_multiplier": g},
                    }
                    for seed, (status, g) in entries.items()
                }
            }
        )
    )
    return path


def test_identity_and_comparison_with_the_frozen_workflow(tmp_path):
    frozen = workflow(
        tmp_path / "frozen.json",
        {
            "YV.RR52.00.BDH": ("candidate", 0.6787278140712512),
            "YV.RR03.00.BDH": ("candidate", 1.1),
            "YV.RR12.00.BDH": ("unresolved", ""),
        },
    )
    new = workflow(
        tmp_path / "new.json",
        {
            "YV.RR52.00.BDH": ("candidate", 0.6787278140712512),
            "YV.RR03.00.BDH": ("candidate", 1.21),
            "YV.RR12.00.BDH": ("unresolved", ""),
        },
    )
    assert rfs.identity_value(new, "YV.RR52.00.BDH") == "0.6787278140712512"
    assert rfs.identity_value(new, "YV.RR99.00.BDH") is None
    result = rfs.compare_workflows(new, frozen)
    assert result["stations_compared"] == 3 and result["identical"] == 2
    assert result["estimates_in_both"] == 2
    assert result["max_abs_gain_change"] == pytest.approx(0.1)
    assert result["differences"][0][0] == "YV.RR03.00.BDH"


def task(campaign, name, kind, status, outputs):
    folder = campaign / "tasks" / name
    folder.mkdir(parents=True)
    (folder / "state.json").write_text(
        json.dumps(
            {
                "kind": kind,
                "status": status,
                "outputs": [{"path": str(p), "sha256": "x"} for p in outputs],
            }
        )
    )


def product(folder, seed, units="Pa^2/Hz"):
    folder.mkdir(parents=True)
    (folder / "ppsd.npz").write_bytes(b"npz")
    (folder / "processing.json").write_text(
        json.dumps(
            {
                "config": {
                    "seed_id": seed,
                    "epoch_start": str(START),
                    "epoch_end": str(END),
                    "psd_units": units,
                },
                "source_inventory": "/inv/YV.xml",
            }
        )
    )
    return folder / "ppsd.npz"


def test_registry_lists_completed_ppsds_with_metadata(tmp_path):
    campaign = tmp_path / "campaign"
    base = product(tmp_path / "p/base", "YV.RR99.00.BDH")
    merged = product(tmp_path / "p/merged", "YV.RR99.00.BDH")
    vertical = product(tmp_path / "p/z", "YV.RR99.00.BHZ", "m^2/s^4/Hz")
    failed = product(tmp_path / "p/failed", "YV.RR98.00.BDH")
    task(campaign, "a", "ppsd", "complete", [base])
    task(campaign, "b", "ppsd_merge", "complete", [merged])
    task(campaign, "c", "ppsd", "complete", [vertical])
    task(campaign, "d", "ppsd", "failed", [failed])
    task(campaign, "e", "expected_pressure", "complete", [tmp_path / "x.json"])
    master = tmp_path / "master.csv"
    master.write_text(
        "Seed_id,channel_start_time,channel_end_time,sensor_type\n"
        f"YV.RR99.00.BDH,{START},{END},DPG\n"
    )
    rows = rfs.registry_rows(campaign, [master])
    assert [r["Seed_id"] for r in rows] == ["YV.RR99.00.BDH", "YV.RR99.00.BHZ"]
    assert rows[0]["ppsd_npz_file"] == str(merged)
    assert rows[0]["sensor_type"] == "DPG" and rows[0]["ppsd_status"] == "complete"
    assert rows[1]["ppsd_units"] == "m^2/s^4/Hz"
    path = rfs.write_registry(rows, tmp_path / "registry/master_inventory.csv")
    assert path.read_text().count("\n") == 3


def test_seeded_prepare_inputs_must_match_their_sources(tmp_path):
    response = tmp_path / "7D.J67A.xml"
    response.write_text("<xml/>")
    catalog = tmp_path / "catalog.json"
    catalog.write_text("{}")
    network = {"local_responses": [str(response)], "local_catalog": str(catalog)}
    workspace = tmp_path / "workspace"
    rfs.seed_prepare_inputs(network, workspace)
    rfs.seed_prepare_inputs(network, workspace)  # idempotent
    assert (workspace / "_inputs/responses/7D.J67A.xml").read_text() == "<xml/>"
    (workspace / "_inputs/catalog_usgs_m5p5.json").write_text("changed")
    with pytest.raises(ValueError, match="catalogue differs"):
        rfs.seed_prepare_inputs(network, workspace)


def test_run_record_survives_reload(tmp_path):
    run = rfs.Run(tmp_path / "run")
    run.data.update(code={"state": "abc"})
    run.mark("workspaces", "complete", YV={"written": 2})
    assert run.done("workspaces") and not run.done("identity")
    again = rfs.Run(tmp_path / "run")
    assert again.stage("workspaces")["YV"] == {"written": 2}


def test_stage_selection_respects_only_and_skip():
    config = rfs.load_config(CONFIG)
    args = SimpleNamespace(networks=None, stations=None, only=None, skip=["coverage"])
    assert "coverage" not in rfs.Pipeline(args, config).selected_stages()
    args = SimpleNamespace(networks=["YV"], stations=None, only=["summary"], skip=None)
    assert rfs.Pipeline(args, config).selected_stages() == ["summary"]
    args = SimpleNamespace(networks=["ZZ"], stations=None, only=None, skip=None)
    with pytest.raises(SystemExit):
        rfs.Pipeline(args, config)


def test_prepare_scripts_never_download_offline(monkeypatch):
    from scripts import prepare_7d_run, prepare_xo_run

    for module in (prepare_7d_run, prepare_xo_run):
        with pytest.raises(SystemExit, match="--offline"):
            module.download("https://example.invalid", Path("t.xml"), True)

        class Reply:
            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

            def read(self):
                return b"data"

        monkeypatch.setattr(module, "urlopen", lambda url, timeout: Reply())
        assert (
            module.download("https://example.invalid", Path("t.xml"), False) == b"data"
        )


def test_microseism_sources_override_and_refuse_bad_input(tmp_path):
    from scripts import run_microseism_check as micro

    original = {k: dict(v) for k, v in micro.NETWORKS.items()}
    try:
        master, data = tmp_path / "master", tmp_path / "data"
        master.mkdir()
        data.mkdir()
        recorded = micro.apply_sources({"XO": {"master": master, "data": data}})
        assert micro.NETWORKS["XO"]["master"] == master.resolve()
        assert micro.NETWORKS["XO"]["resample_hz"] == 1.0  # completed-run setting
        assert recorded["XO"]["data"] == str(data.resolve())
        with pytest.raises(FileNotFoundError):
            micro.apply_sources({"YV": {"master": tmp_path / "no", "data": data}})
        with pytest.raises(ValueError):
            micro.apply_sources({"ZZ": {"master": master, "data": data}})
    finally:
        micro.NETWORKS.clear()
        micro.NETWORKS.update(original)
