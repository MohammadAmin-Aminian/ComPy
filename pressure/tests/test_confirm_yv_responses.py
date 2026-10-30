"""Official-versus-saved response comparison: provenance, roots and refusals."""

import io
import json
import urllib.error
from types import SimpleNamespace

import numpy as np
import pytest
from obspy import UTCDateTime
from obspy.core.inventory import Channel, Inventory, Network, Site, Station
from obspy.core.inventory.response import Response

from scripts.confirm_yv_responses import (
    STATIONS_2016,
    channel_summary,
    compare_summaries,
    fetch,
    hydrophone_context,
    literature_check,
    load_official,
    match_roots,
    normalization_check,
    planned_queries,
    pz_a0,
    run_compare,
    run_fetch,
)
from scripts.diagnose_yv_rayleigh import sha256

P = -0.074016 + 0.07347j
Q = -502.65 + 596.9j
GOOD_Z = [P, P.conjugate(), Q, Q.conjugate()]
BAD_Z = [P, P, Q, Q]
START = UTCDateTime("2012-10-20T18:27:45")
END = UTCDateTime("2013-06-29T20:55:33")


def make_channel(code, poles, zeros, a0=None, units="M/S", end=END, gain=1963.0):
    if a0 is None:
        a0 = pz_a0(GOOD_Z if poles == BAD_Z else poles, zeros, 1.0)
    response = Response.from_paz(
        zeros=zeros,
        poles=poles,
        stage_gain=gain,
        input_units=units,
        output_units="V",
        normalization_factor=a0,
    )
    # Stated metadata sensitivity, not ObsPy's recalculation from the roots.
    response.instrument_sensitivity.value = gain
    return Channel(
        code=code,
        location_code="00",
        latitude=0,
        longitude=0,
        elevation=-4000,
        depth=0,
        sample_rate=50,
        start_date=START,
        end_date=end,
        response=response,
    )


def write_inventory(path, station, channels):
    site = Station(station, 0, 0, -4000, channels=channels, site=Site(station))
    Inventory([Network("YV", stations=[site])]).write(str(path), "STATIONXML")


def test_match_roots_respects_sign_count_and_tolerance():
    assert match_roots(GOOD_Z, list(reversed(GOOD_Z))) == (True, 0.0)
    assert not match_roots(BAD_Z, GOOD_Z)[0]
    assert match_roots(GOOD_Z[:3], GOOD_Z) == (False, None)
    nudged = [r * (1 + 1e-6) for r in GOOD_Z]
    assert not match_roots(nudged, GOOD_Z)[0]
    assert match_roots(nudged, GOOD_Z, rtol=1e-4)[0]


def test_normalization_check_identifies_the_paired_root_set():
    channel = make_channel("BHZ", BAD_Z, [0j, 0j])
    stage = channel.response.response_stages[0]
    check = normalization_check(stage)
    assert abs(check["paired_relative_difference"]) < 1e-12
    assert check["stored_relative_difference"] < -0.03
    good = normalization_check(
        make_channel("BHZ", GOOD_Z, [0j, 0j]).response.response_stages[0]
    )
    assert "from_paired_roots" not in good
    assert abs(good["stored_relative_difference"]) < 1e-12


def test_compare_summaries_reports_only_real_differences():
    base = channel_summary(make_channel("BHZ", GOOD_Z, [0j, 0j]))
    assert (
        compare_summaries(base, channel_summary(make_channel("BHZ", GOOD_Z, [0j, 0j])))
        == []
    )
    flipped = channel_summary(make_channel("BHZ", BAD_Z, [0j, 0j]))
    assert [d["field"] for d in compare_summaries(base, flipped)] == ["stage1.poles"]
    other_units = channel_summary(make_channel("BHZ", GOOD_Z, [0j, 0j], units="m/s"))
    fields = [d["field"] for d in compare_summaries(base, other_units)]
    assert "sensitivity.input_units" in fields and "stage1.input_units" in fields
    shorter = channel_summary(make_channel("BHZ", GOOD_Z, [0j, 0j], end=END - 1))
    assert [d["field"] for d in compare_summaries(base, shorter)] == ["end"]


def test_literature_check_separates_stored_and_paired_roots():
    rr56 = literature_check("RR56", make_channel("BHZ", BAD_Z, [0j, 0j]))
    assert rr56["table"] == "A1a"
    assert rr56["poles_match_as_stored"] is False
    assert rr56["poles_match_after_pairing"] is True
    rr17 = literature_check("RR17", make_channel("BHZ", GOOD_Z, [0j, 0j]))
    assert rr17["poles_match_as_stored"] is True
    assert "poles_match_after_pairing" not in rr17
    p = -0.074016 + 0.074016j
    q = -395.1 + 850.69j
    rr22 = [p, p.conjugate(), -471.24 + 0j, q.conjugate(), q.conjugate(), -2199.1 + 0j]
    check = literature_check("RR22", make_channel("BHZ", rr22, [0j, 0j]))
    assert check["poles_match_as_stored"] is False
    assert check["poles_match_after_pairing"] is True
    assert (
        literature_check("RR22", make_channel("BDH", [0.099 + 0j], [0j], units="PA"))
        is None
    )
    assert literature_check("RR99", make_channel("BHZ", GOOD_Z, [0j, 0j])) is None


def test_hydrophone_context_reports_implied_time_constant():
    context = hydrophone_context(
        "RR17", make_channel("BDH", [0.747005 + 0j], [0j], units="PA", gain=0.000191)
    )
    assert context["family"] == STATIONS_2016["RR17"]["hydrophone_family"]
    assert context["implied_R_C_total_s"] == pytest.approx(1.33868, rel=1e-5)
    assert context["implied_C_total_nF_at_R_210_Mohm"] == pytest.approx(
        6.3747, rel=1e-4
    )
    assert context["stage_sensitivity_uV_per_Pa"] == pytest.approx(191.0)


def test_planned_queries_are_limited_to_requested_streams():
    queries = planned_queries(["RR17", "RR56"], ["BDH", "BHZ"])
    resif = [q for q in queries if q["name"].startswith("resif_")]
    assert [q["station"] for q in resif] == ["RR17", "RR56"]
    for query in resif:
        assert query["url"].startswith("https://ws.resif.fr/fdsnws/station/1/query?")
        assert "level=response" in query["url"]
        assert "channel=BDH%2CBHZ" in query["url"]
        assert "%2A" not in query["url"] and "*" not in query["url"]
    assert {q["name"] for q in queries} == {
        "eida_routing",
        "fedcatalog",
        "resif_RR17",
        "resif_RR56",
    }


class FakeResponse(io.BytesIO):
    def __init__(self, body, status=200):
        super().__init__(body)
        self.status = status
        self.headers = {"Content-Type": "application/xml", "Date": "today"}


def opener_for(outcomes):
    def opener(request, timeout):
        outcome = outcomes[request.full_url.split("?")[0]]
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    return opener


def test_fetch_records_success_no_data_and_failures():
    query = planned_queries(["RR17"], ["BDH"])[-1]
    url = query["url"].split("?")[0]
    record, body = fetch(query, opener_for({url: FakeResponse(b"<xml/>")}))
    assert record["outcome"] == "saved" and body == b"<xml/>"
    assert record["sha256"] == sha256_bytes(b"<xml/>")
    assert record["bytes"] == 6 and record["headers"]["Date"] == "today"
    record, body = fetch(query, opener_for({url: FakeResponse(b"", status=204)}))
    assert record["outcome"] == "no_data" and body is None
    error = urllib.error.HTTPError(url, 503, "busy", {}, None)
    record, body = fetch(query, opener_for({url: error}))
    assert record["outcome"] == "failed" and record["http_status"] == 503
    record, body = fetch(query, opener_for({url: urllib.error.URLError("offline")}))
    assert record["outcome"] == "failed" and record["http_status"] is None


def sha256_bytes(data):
    import hashlib

    return hashlib.sha256(data).hexdigest()


def test_run_fetch_is_exclusive_and_keeps_failures(tmp_path):
    queries = planned_queries(["RR17"], ["BDH"])
    outcomes = {q["url"].split("?")[0]: FakeResponse(b"x") for q in queries}
    outcomes["https://ws.resif.fr/fdsnws/station/1/query"] = urllib.error.URLError(
        "down"
    )
    args = SimpleNamespace(
        out=tmp_path / "official", stations=["RR17"], channels=["BDH"], timeout=1
    )
    manifest = run_fetch(args, opener_for(outcomes))
    names = {q["name"]: q for q in manifest["queries"]}
    assert names["resif_RR17"]["outcome"] == "failed"
    assert not (tmp_path / "official" / "YV.RR17.response.xml").exists()
    assert (tmp_path / "official" / "eida_routing.json").read_bytes() == b"x"
    with pytest.raises(FileExistsError):
        run_fetch(args, opener_for(outcomes))


def test_load_official_refuses_changed_download(tmp_path):
    (tmp_path / "a.xml").write_bytes(b"one")
    record = {"name": "resif_RR17", "outcome": "saved", "file": "a.xml"}
    record["sha256"] = sha256_bytes(b"one")
    (tmp_path / "manifest.json").write_text(json.dumps({"queries": [record]}))
    assert load_official(tmp_path)[1]["resif_RR17"]["sha256"] == record["sha256"]
    (tmp_path / "a.xml").write_bytes(b"two")
    with pytest.raises(ValueError, match="changed since retrieval"):
        load_official(tmp_path)


def synthetic_station(tmp_path, official_z):
    """Saved RR56-like products plus an official download."""
    pressure = make_channel("BDH", [0.102613 + 0j], [0j], units="PA", gain=0.0002)
    vertical = make_channel("BHZ", BAD_Z, [0j, 0j])
    source = tmp_path / "YV.xml"
    write_inventory(source, "RR56", [pressure, vertical])
    epoch_keys, snapshots = [], {}
    for channel in (pressure, vertical):
        seed_id = f"YV.RR56.00.{channel.code}"
        folder = tmp_path / "saved" / f"{seed_id}__20121020T182745__abc"
        folder.mkdir(parents=True)
        write_inventory(folder / "inventory.xml", "RR56", [channel])
        snapshots[str(folder / "inventory.xml")] = sha256(folder / "inventory.xml")
        epoch_keys.append(
            {
                "Seed_id": seed_id,
                "channel_start_time": str(START),
                "channel_end_time": str(END),
                "inventory_file": str(source),
            }
        )
    comparison = tmp_path / "comparison.json"
    comparison.write_text(
        json.dumps(
            {
                "epoch_keys": epoch_keys,
                "source_inventory_sha256": {str(source): sha256(source)},
                "response_snapshot_sha256": snapshots,
            }
        )
    )
    workflow = tmp_path / "rayleigh.json"
    entry = {
        "seed_id": "YV.RR56.00.BDH",
        "update": {"rayleigh_result_file": str(comparison)},
    }
    workflow.write_text(json.dumps({"entries": {"x": entry}}))
    official = tmp_path / "official"
    official.mkdir()
    official_vertical = make_channel(
        "BHZ", official_z, [0j, 0j], a0=pz_a0(GOOD_Z, [0j, 0j], 1)
    )
    write_inventory(
        official / "YV.RR56.response.xml", "RR56", [pressure, official_vertical]
    )
    record = {
        "name": "resif_RR56",
        "outcome": "saved",
        "file": "YV.RR56.response.xml",
        "sha256": sha256(official / "YV.RR56.response.xml"),
    }
    (official / "manifest.json").write_text(json.dumps({"queries": [record]}))
    return SimpleNamespace(
        workflow=str(workflow),
        workflow_sha256=sha256(workflow),
        official=str(official),
        out=str(tmp_path / "out"),
        stations=["RR56"],
    ), snapshots


def test_run_compare_binds_hashes_and_reports_root_defects(tmp_path):
    args, _ = synthetic_station(tmp_path, BAD_Z)
    evidence = run_compare(args)
    channels = evidence["stations"]["RR56"]["channels"]
    vertical = channels["YV.RR56.00.BHZ"]
    assert vertical["official_differences"] == []
    assert vertical["source_inventory_differences"] == []
    assert len(vertical["local_physical"][0]["unpaired_poles"]) == 4
    assert vertical["literature"]["poles_match_after_pairing"] is True
    pressure = channels["YV.RR56.00.BDH"]
    assert pressure["local_physical"][0]["positive_real_part_poles"] == [
        [0.102613, 0.0]
    ]
    terms = evidence["stations"]["RR56"]["hypothesis_band_terms"]
    assert terms["H0"]["phase_deg"] == [0.0, 0.0, 0.0]
    assert np.all(np.array(terms["H2"]["amplitude_factor"]) > 1.3)
    with pytest.raises(FileExistsError):
        run_compare(args)


def test_run_compare_reports_a_corrected_official_response(tmp_path):
    args, _ = synthetic_station(tmp_path, GOOD_Z)
    vertical = run_compare(args)["stations"]["RR56"]["channels"]["YV.RR56.00.BHZ"]
    assert [d["field"] for d in vertical["official_differences"]] == ["stage1.poles"]
    assert vertical["official"]["evalresp_over_sensitivity"] == pytest.approx(1.0)
    assert vertical["local"]["evalresp_over_sensitivity"] > 1.03
    assert vertical["official_physical"][0]["unpaired_poles"] == []


def test_run_compare_refuses_changed_snapshot_and_workflow(tmp_path):
    args, snapshots = synthetic_station(tmp_path, BAD_Z)
    first = next(iter(snapshots))
    with open(first, "a", encoding="utf-8") as handle:
        handle.write("\n")
    with pytest.raises(ValueError, match="recorded value"):
        run_compare(args)
    args.workflow_sha256 = "0" * 64
    with pytest.raises(ValueError, match="recorded value"):
        run_compare(args)
