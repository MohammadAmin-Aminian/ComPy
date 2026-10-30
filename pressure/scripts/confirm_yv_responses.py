"""Compare official FDSN response metadata with saved YV P/Z response snapshots.

Two read-only steps:

* ``fetch`` downloads public StationXML from fixed FDSN queries into a new,
  exclusive folder. It records provider, exact URL, retrieval UTC, HTTP status,
  selected headers, byte count, SHA-256 and any failure. It never writes over
  an existing file or inventory.
* ``compare`` binds each station's frozen no-tilt workflow entry to its saved
  P/Z response snapshots (recorded hashes re-verified), the source inventory
  channel and the official channel epoch with the same start and end. It
  compares epoch, sample rate, units, stage gains, normalization, transfer
  function type, poles/zeros, decimation and filter coefficients. It also
  reports physical-validity checks (conjugate pairs, positive real poles, A0
  consistency with the stored roots) and the primary-literature values.

Archive agreement is metadata provenance, not physical calibration. Nothing is
repaired, fitted, applied or written outside the new output folder.
"""

import argparse
import csv
import hashlib
import json
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from scripts.diagnose_yv_rayleigh import conjugate_repair, sha256, unpaired_roots

RESIF_STATION = "https://ws.resif.fr/fdsnws/station/1/query"
EIDA_ROUTING = "https://www.orfeus-eu.org/eidaws/routing/1/query"
FEDCATALOG = "https://service.iris.edu/irisws/fedcatalog/1/query"
# Same window as the ModuleURI recorded in the local RESIF export of YV.xml.
QUERY_START = "2012-01-01T00:00:00"
QUERY_END = "2014-01-01T00:00:00"
USER_AGENT = "bruit-fm-response-confirmation/1 (read-only metadata review)"
KEPT_HEADERS = ("Date", "Last-Modified", "Content-Type", "Content-Length", "ETag")
BAND_FREQUENCIES_HZ = (0.03, 0.05, 0.07)
ROOT_RTOL = 1e-9
VALUE_RTOL = 1e-9
LITERATURE_RTOL = 1e-4

STAEHLER_2016 = (
    "Stähler et al. (2016), Adv. Geosci. 41, 43-63, doi:10.5194/adgeo-41-43-2016"
)
P_A1A = -0.074016 + 0.07347j
Q_A1A = -502.65 + 596.9j
P_A1C = -0.074016 + 0.074016j
Q_A1C = -395.1 + 850.69j
TABLE_A1 = {
    "A1a": {
        "label": "60 s Güralp CMG-OBS40T, German LOBSTER OBS; "
        "caption: all 60 s stations but RR13 and RR22",
        "poles": [P_A1A, P_A1A.conjugate(), Q_A1A, Q_A1A.conjugate()],
        "zeros": [0j, 0j],
    },
    "A1c": {
        "label": "60 s Güralp CMG-OBS40T used at RR22",
        "poles": [
            P_A1C,
            P_A1C.conjugate(),
            -471.24 + 0j,
            Q_A1C,
            Q_A1C.conjugate(),
            -2199.1 + 0j,
        ],
        "zeros": [0j, 0j],
    },
}
# Appendix A2 and Table 3, transcribed from the primary PDF; values as printed.
HYDROPHONE_FAMILIES = {
    "DEPAS": {
        "sensitivity_uV_per_Pa": "about 185, sample standard deviation 8",
        "pole": "p1 = -1/(R C_total), C_total = C_amp C_hyd/(C_amp + C_hyd)",
        "C_hyd_nF": "manufacturer reference 45; "
        "calibrated mean 56.3, sd 3.5 (60 units)",
        "R_Mohm": "210 or 500, depending on instrument version",
        "C_amp": "not printed",
        "station_values": "supplied by DEPAS per instrument; not printed",
    },
    "Geomar HTI-01-PCA": {
        "sensitivity_uV_per_Pa": "199.5 (DEPAS HTI-01-PCA average; "
        "no individual calibration)",
        "pole": "p1 = 0.10774 rad/s as printed (C_hyd nominal 50 nF)",
        "stations": "RR33, RR39, RR53, RR56, RR45, RR55",
    },
}
STATIONS_2016 = {
    "RR17": {
        "vertical_table": "A1a",
        "hydrophone_family": "DEPAS",
        "table3": "2205 m; deployed 3 Oct 2012; recovered 23 Oct 2013; "
        "end of record 23 Oct 2013; seismometer good; hydrophone good; skew 1.82 s",
    },
    "RR18": {
        "vertical_table": "A1a",
        "hydrophone_family": "DEPAS",
        "table3": "4743 m; deployed 6 Oct 2012; recovered 29 Oct 2013; "
        "end of record 29 Oct 2013; seismometer good; hydrophone good; skew 0.36 s",
    },
    "RR22": {
        "vertical_table": "A1c",
        "hydrophone_family": "DEPAS",
        "table3": "4920 m; deployed 9 Oct 2012; recovered 1 Nov 2013; "
        "end of record 1 Nov 2013; seismometer good; hydrophone good; skew 0.89 s",
    },
    "RR56": {
        "vertical_table": "A1a",
        "hydrophone_family": "Geomar HTI-01-PCA",
        "table3": "4230 m; deployed 21 Oct 2012; recovered 29 Nov 2013; "
        "end of record 29 Jun 2013; seismometer good; hydrophone good; skew NA; Geomar",
    },
}


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def planned_queries(stations, channels, start=QUERY_START, end=QUERY_END):
    """Fixed official queries: routing lookups and per-station response XML."""
    station_list = ",".join(stations)
    channel_list = ",".join(channels)
    queries = [
        {
            "name": "eida_routing",
            "provider": "EIDA routing service (ORFEUS)",
            "purpose": "which EIDA data centre serves these YV streams",
            "url": EIDA_ROUTING
            + "?"
            + urllib.parse.urlencode(
                {
                    "network": "YV",
                    "station": station_list,
                    "channel": channel_list,
                    "starttime": start,
                    "endtime": end,
                    "service": "station",
                    "format": "json",
                }
            ),
            "filename": "eida_routing.json",
        },
        {
            "name": "fedcatalog",
            "provider": "EarthScope/IRIS FDSN federator catalogue",
            "purpose": "which FDSN data centres list these YV streams",
            "url": FEDCATALOG
            + "?"
            + urllib.parse.urlencode(
                {
                    "net": "YV",
                    "sta": station_list,
                    "cha": channel_list,
                    "starttime": start,
                    "endtime": end,
                    "format": "text",
                }
            ),
            "filename": "fedcatalog.txt",
        },
    ]
    for station in stations:
        params = {
            "network": "YV",
            "station": station,
            "channel": channel_list,
            "starttime": start,
            "endtime": end,
            "level": "response",
            "format": "xml",
        }
        queries.append(
            {
                "name": f"resif_{station}",
                "provider": "RESIF data centre FDSN station service",
                "purpose": f"official StationXML responses for YV.{station}",
                "url": RESIF_STATION + "?" + urllib.parse.urlencode(params),
                "filename": f"YV.{station}.response.xml",
                "station": station,
            }
        )
    return queries


def fetch(query, opener=urllib.request.urlopen, timeout=60):
    """Run one GET; return (provenance record, body bytes or None)."""
    record = {k: query[k] for k in query if k != "filename"}
    record["method"] = "GET"
    record["request_utc"] = utc_now()
    request = urllib.request.Request(query["url"], headers={"User-Agent": USER_AGENT})
    body = None
    try:
        with opener(request, timeout=timeout) as response:
            body = response.read()
            record["http_status"] = int(response.status)
            headers = response.headers
    except urllib.error.HTTPError as error:
        record["http_status"] = int(error.code)
        record["error"] = f"HTTPError {error.code}: {error.reason}"
        headers = error.headers
    except (urllib.error.URLError, TimeoutError, OSError) as error:
        record["http_status"] = None
        record["error"] = f"{type(error).__name__}: {error}"
        headers = None
    record["response_utc"] = utc_now()
    record["headers"] = {
        k: headers.get(k) for k in KEPT_HEADERS if headers and headers.get(k)
    }
    if record.get("http_status") == 204:
        record["outcome"] = "no_data"
        body = None
    elif body is not None and record.get("http_status") == 200:
        record["outcome"] = "saved"
        record["bytes"] = len(body)
        record["sha256"] = hashlib.sha256(body).hexdigest()
    else:
        record["outcome"] = "failed"
        body = None
    return record, body


def run_fetch(args, opener=urllib.request.urlopen):
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.mkdir(exist_ok=False)
    manifest = {
        "tool": "scripts/confirm_yv_responses.py fetch",
        "code_sha256": sha256(Path(__file__)),
        "created_utc": utc_now(),
        "note": "Fresh official downloads; never written over the local YV.xml.",
        "queries": [],
    }
    for query in planned_queries(args.stations, args.channels):
        record, body = fetch(query, opener, args.timeout)
        if body is not None:
            with open(out / query["filename"], "xb") as handle:
                handle.write(body)
            record["file"] = query["filename"]
        manifest["queries"].append(record)
        print(f"{record['name']}: {record['outcome']} {record.get('http_status')}")
    with open(out / "manifest.json", "x", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    return manifest


def roots(values):
    return [complex(v) for v in values or []]


def root_pairs(values):
    return [[v.real, v.imag] for v in values]


def match_roots(a, b, rtol=ROOT_RTOL):
    """Greedy nearest multiset match; returns (matched, max relative difference)."""
    a, b = roots(a), roots(b)
    if len(a) != len(b):
        return False, None
    free = list(b)
    worst = 0.0
    for r in sorted(a, key=lambda x: (x.real, x.imag)):
        j = int(np.argmin([abs(r - x) for x in free]))
        scale = max(abs(r), abs(free[j]))
        worst = max(worst, abs(r - free[j]) / scale if scale else 0.0)
        free.pop(j)
    return worst <= rtol, worst


def pz_a0(poles, zeros, frequency):
    """Normalization making |A0 * prod(s - z) / prod(s - p)| = 1 at s = 2*pi*i*f."""
    s = 2j * np.pi * float(frequency)
    value = np.prod([s - z for z in zeros]) / np.prod([s - p for p in poles])
    return float(1.0 / abs(value))


def normalization_check(stage):
    """Recorded A0 against A0 recomputed from the stored and paired roots."""
    if stage.pz_transfer_function_type != "LAPLACE (RADIANS/SECOND)":
        return {"checked": False, "reason": stage.pz_transfer_function_type}
    poles, zeros = roots(stage.poles), roots(stage.zeros)
    frequency = stage.normalization_frequency
    recorded = float(stage.normalization_factor)
    stored = pz_a0(poles, zeros, frequency)
    result = {
        "checked": True,
        "recorded": recorded,
        "from_stored_roots": stored,
        "stored_relative_difference": stored / recorded - 1.0,
    }
    if unpaired_roots(poles) or unpaired_roots(zeros):
        try:
            paired = pz_a0(conjugate_repair(poles), conjugate_repair(zeros), frequency)
        except ValueError as error:
            result["from_paired_roots"] = None
            result["pairing_error"] = str(error)
        else:
            result["from_paired_roots"] = paired
            result["paired_relative_difference"] = paired / recorded - 1.0
    return result


def coefficient_hash(values):
    data = np.asarray([float(v) for v in values or []], dtype="<f8")
    return hashlib.sha256(data.tobytes()).hexdigest(), int(data.size)


def stage_summary(stage):
    summary = {
        "sequence": stage.stage_sequence_number,
        "class": type(stage).__name__,
        "input_units": stage.input_units,
        "output_units": stage.output_units,
        "gain": stage.stage_gain,
        "gain_frequency": stage.stage_gain_frequency,
        "decimation": [
            getattr(stage, f"decimation_{k}", None)
            for k in ("input_sample_rate", "factor", "offset", "delay", "correction")
        ],
    }
    if hasattr(stage, "poles"):
        summary.update(
            {
                "transfer_function_type": stage.pz_transfer_function_type,
                "normalization_factor": stage.normalization_factor,
                "normalization_frequency": stage.normalization_frequency,
                "poles": root_pairs(roots(stage.poles)),
                "zeros": root_pairs(roots(stage.zeros)),
            }
        )
    for name in ("numerator", "denominator", "coefficients"):
        if hasattr(stage, name):
            summary[name] = coefficient_hash(getattr(stage, name))
    if hasattr(stage, "symmetry"):
        summary["symmetry"] = stage.symmetry
    if hasattr(stage, "cf_transfer_function_type"):
        summary["transfer_function_type"] = stage.cf_transfer_function_type
    return summary


def channel_summary(channel):
    response = channel.response
    sensitivity = response.instrument_sensitivity
    stages = [stage_summary(s) for s in response.response_stages]
    product = float(np.prod([s["gain"] for s in stages]))
    reference = sensitivity.frequency
    evalresp = response.get_evalresp_response_for_frequencies([reference], output="DEF")
    sensor = channel.sensor
    return {
        "start": str(channel.start_date),
        "end": str(channel.end_date) if channel.end_date else None,
        "location": channel.location_code,
        "code": channel.code,
        "sample_rate": channel.sample_rate,
        "sensor": {
            k: getattr(sensor, k, None)
            for k in ("type", "description", "manufacturer", "model", "serial_number")
        },
        "sensitivity": {
            "value": sensitivity.value,
            "frequency": reference,
            "input_units": sensitivity.input_units,
            "output_units": sensitivity.output_units,
        },
        "stage_gain_product": product,
        "sensitivity_over_stage_product": sensitivity.value / product,
        "evalresp_abs_at_reference": float(abs(evalresp[0])),
        "evalresp_over_sensitivity": float(abs(evalresp[0]) / sensitivity.value),
        "stages": stages,
    }


def physical_checks(channel):
    """Root validity of every analogue stage; no repair is applied."""
    checks = []
    for stage in channel.response.response_stages:
        if not hasattr(stage, "poles"):
            continue
        poles, zeros = roots(stage.poles), roots(stage.zeros)
        checks.append(
            {
                "sequence": stage.stage_sequence_number,
                "unpaired_poles": root_pairs(unpaired_roots(poles)),
                "unpaired_zeros": root_pairs(unpaired_roots(zeros)),
                "positive_real_part_poles": root_pairs(
                    [p for p in poles if p.real > 0]
                ),
                "normalization": normalization_check(stage),
            }
        )
    return checks


def _close(a, b, rtol=VALUE_RTOL):
    if a is None or b is None:
        return a is b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return bool(np.isclose(a, b, rtol=rtol, atol=0.0))
    return a == b


def compare_summaries(local, official, rtol=ROOT_RTOL):
    """Field-by-field differences between two channel summaries."""
    differences = []

    def note(field, a, b):
        differences.append({"field": field, "local": a, "official": b})

    for key in ("start", "end", "location", "code", "sample_rate"):
        if not _close(local[key], official[key]):
            note(key, local[key], official[key])
    for key, value in local["sensitivity"].items():
        if not _close(value, official["sensitivity"][key]):
            note(f"sensitivity.{key}", value, official["sensitivity"][key])
    if len(local["stages"]) != len(official["stages"]):
        note("stage_count", len(local["stages"]), len(official["stages"]))
        return differences
    for a, b in zip(local["stages"], official["stages"], strict=True):
        prefix = f"stage{a['sequence']}"
        for key in sorted(set(a) | set(b)):
            if key in ("poles", "zeros"):
                ra = [complex(*x) for x in a.get(key, [])]
                rb = [complex(*x) for x in b.get(key, [])]
                matched, worst = match_roots(ra, rb, rtol)
                if not matched:
                    note(f"{prefix}.{key}", a.get(key), b.get(key))
            elif key == "decimation":
                if not all(_close(x, y) for x, y in zip(a[key], b[key], strict=True)):
                    note(f"{prefix}.{key}", a[key], b[key])
            elif not _close(a.get(key), b.get(key)):
                note(f"{prefix}.{key}", a.get(key), b.get(key))
    return differences


def literature_check(station, channel):
    """Compare the analogue vertical stage to the station's Table A1 entry."""
    expected = STATIONS_2016.get(station)
    if expected is None or channel.code[-1] != "Z":
        return None
    table = TABLE_A1[expected["vertical_table"]]
    stage = next(s for s in channel.response.response_stages if hasattr(s, "poles"))
    poles, zeros = roots(stage.poles), roots(stage.zeros)
    as_stored, stored_diff = match_roots(poles, table["poles"], LITERATURE_RTOL)
    result = {
        "table": expected["vertical_table"],
        "label": table["label"],
        "poles_match_as_stored": as_stored,
        "max_relative_difference_as_stored": stored_diff,
        "zeros_match": match_roots(zeros, table["zeros"], LITERATURE_RTOL)[0],
    }
    if unpaired_roots(poles):
        try:
            paired = conjugate_repair(poles)
        except ValueError:
            result["poles_match_after_pairing"] = None
        else:
            matched, diff = match_roots(paired, table["poles"], LITERATURE_RTOL)
            result["poles_match_after_pairing"] = matched
            result["max_relative_difference_after_pairing"] = diff
    return result


def hydrophone_context(station, channel):
    expected = STATIONS_2016.get(station)
    if expected is None or channel.code[-1] != "H":
        return None
    stage = next(s for s in channel.response.response_stages if hasattr(s, "poles"))
    pole = max(roots(stage.poles), key=abs)
    rc = 1.0 / abs(pole)
    return {
        "family": expected["hydrophone_family"],
        "printed": HYDROPHONE_FAMILIES[expected["hydrophone_family"]],
        "stage_sensitivity_uV_per_Pa": stage.stage_gain * 1e6,
        "pole_rad_s": [pole.real, pole.imag],
        "implied_R_C_total_s": rc,
        "implied_C_total_nF_at_R_210_Mohm": rc / 210e6 * 1e9,
        "implied_C_total_nF_at_R_500_Mohm": rc / 500e6 * 1e9,
    }


def band_terms(terms_function, pressure, vertical):
    """Fixed-hypothesis phase and amplitude factor in the Rayleigh band.

    Deterministic metadata arithmetic from the existing diagnosis helper: the
    phase nominal processing would show, and the factor on the saved |P|/|A|
    ratio, if each hypothesis were true. No saved spectrum is read.
    """
    frequency = np.asarray(BAND_FREQUENCIES_HZ)
    try:
        terms = terms_function(frequency, pressure, vertical)
    except ValueError as error:
        return {"error": str(error)}
    result = {"frequency_hz": list(BAND_FREQUENCIES_HZ)}
    for name, (phase, amplitude) in terms.items():
        result[name] = {
            "phase_deg": [float(x) for x in phase],
            "amplitude_factor": [float(x) for x in amplitude],
        }
    return result


def select_channels(inventory, seed_id):
    network, station, location, code = seed_id.split(".")
    return [
        c
        for n in inventory
        if n.code == network
        for s in n
        if s.code == station
        for c in s
        if c.code == code and c.location_code == location
    ]


def exact_epoch(channels, start, end):
    matches = [
        c for c in channels if str(c.start_date) == start and str(c.end_date) == end
    ]
    return matches[0] if len(matches) == 1 else None


def load_official(folder):
    """Official files verified against their fetch manifest."""
    folder = Path(folder)
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    files = {}
    for record in manifest["queries"]:
        if record.get("outcome") != "saved":
            continue
        path = folder / record["file"]
        if sha256(path) != record["sha256"]:
            raise ValueError(f"Official file changed since retrieval: {path}")
        files[record["name"]] = {
            "path": path,
            "sha256": record["sha256"],
            "record": record,
        }
    return manifest, files


def run_compare(args):
    from obspy import read_inventory

    from scripts.diagnose_yv_rayleigh import (
        hypothesis_terms,
        laplace_stages,
        snapshot_channel,
    )

    out = Path(args.out)
    hashes = {}

    def bind(path, expected=None):
        path = Path(path)
        actual = sha256(path)
        if expected is not None and actual != expected:
            raise ValueError(f"Input hash differs from recorded value: {path}")
        hashes[str(path)] = actual
        return actual

    workflow_path = Path(args.workflow)
    bind(workflow_path, args.workflow_sha256)
    workflow = json.loads(workflow_path.read_text(encoding="utf-8"))
    manifest, official_files = load_official(args.official)
    for item in official_files.values():
        bind(item["path"], item["sha256"])
    inventories = {}
    stations = {}
    for station in args.stations:
        entries = [
            e
            for e in workflow["entries"].values()
            if e["seed_id"].split(".")[1] == station
        ]
        if len(entries) != 1:
            raise ValueError(f"Expected one workflow entry for {station}")
        comparison_path = Path(entries[0]["update"]["rayleigh_result_file"])
        bind(comparison_path)
        comparison = json.loads(comparison_path.read_text(encoding="utf-8"))
        for source, digest in comparison["source_inventory_sha256"].items():
            bind(source, digest)
            if source not in inventories:
                inventories[source] = read_inventory(source)
        official_item = official_files.get(f"resif_{station}")
        official_inventory = (
            read_inventory(str(official_item["path"])) if official_item else None
        )
        channels = {}
        loaded = {}
        for key in comparison["epoch_keys"]:
            seed_id = key["Seed_id"]
            snapshot = next(
                p
                for p in comparison["response_snapshot_sha256"]
                if Path(p).parent.name.startswith(seed_id + "__")
            )
            bind(snapshot, comparison["response_snapshot_sha256"][snapshot])
            local = snapshot_channel(snapshot, seed_id, key)
            loaded[local.code[-1]] = local
            local_summary = channel_summary(local)
            source_channel = exact_epoch(
                select_channels(inventories[key["inventory_file"]], seed_id),
                local_summary["start"],
                local_summary["end"],
            )
            record = {
                "seed_id": seed_id,
                "epoch_key": key,
                "snapshot": snapshot,
                "snapshot_sha256": hashes[snapshot],
                "local": local_summary,
                "local_physical": physical_checks(local),
                "source_inventory_differences": (
                    compare_summaries(local_summary, channel_summary(source_channel))
                    if source_channel
                    else "no exact source epoch"
                ),
                "literature": literature_check(station, local),
                "hydrophone": hydrophone_context(station, local),
            }
            if official_inventory is None:
                record["official"] = "no official file"
            else:
                epochs = select_channels(official_inventory, seed_id)
                match = exact_epoch(
                    epochs, local_summary["start"], local_summary["end"]
                )
                record["official_epochs"] = [
                    [str(c.start_date), str(c.end_date)] for c in epochs
                ]
                record["official_file_sha256"] = official_item["sha256"]
                if match is None:
                    record["official"] = "no exact official epoch"
                else:
                    official_summary = channel_summary(match)
                    record["official"] = official_summary
                    record["official_physical"] = physical_checks(match)
                    record["official_differences"] = compare_summaries(
                        local_summary, official_summary
                    )
            channels[seed_id] = record
        stations[station] = {
            "comparison": str(comparison_path),
            "channels": channels,
            "hypothesis_band_terms": band_terms(
                hypothesis_terms,
                laplace_stages(loaded["H"]),
                laplace_stages(loaded["Z"]),
            ),
        }
    for path, digest in hashes.items():
        if sha256(path) != digest:
            raise ValueError(f"Input changed during comparison: {path}")
    evidence = {
        "tool": "scripts/confirm_yv_responses.py compare",
        "code_sha256": sha256(Path(__file__)),
        "created_utc": utc_now(),
        "literature": STAEHLER_2016,
        "official_manifest": manifest,
        "tolerances": {
            "roots_and_values_relative": ROOT_RTOL,
            "literature_roots_relative": LITERATURE_RTOL,
        },
        "input_sha256": hashes,
        "stations": stations,
    }
    out.mkdir(parents=True, exist_ok=False)
    with open(out / "evidence.json", "x", encoding="utf-8") as handle:
        json.dump(evidence, handle, indent=2, default=str)
        handle.write("\n")
    with open(out / "channels.csv", "x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(csv_row("", {}).keys()), lineterminator="\n"
        )
        writer.writeheader()
        for station, data in stations.items():
            for record in data["channels"].values():
                writer.writerow(csv_row(station, record))
    return evidence


def csv_row(station, record):
    local = record.get("local") or {}
    physical = (record.get("local_physical") or [{}])[0]
    norm = physical.get("normalization") or {}
    official = record.get("official")
    differences = record.get("official_differences")
    literature = record.get("literature") or {}
    hydrophone = record.get("hydrophone") or {}
    return {
        "station": station,
        "seed_id": record.get("seed_id", ""),
        "start": local.get("start", ""),
        "end": local.get("end", ""),
        "snapshot_sha256": record.get("snapshot_sha256", ""),
        "source_inventory_identical": record.get("source_inventory_differences") == [],
        "official_epoch": "exact" if isinstance(official, dict) else (official or ""),
        "official_identical": differences == [] if differences is not None else "",
        "official_difference_fields": ";".join(d["field"] for d in differences or []),
        "local_unpaired_poles": len(physical.get("unpaired_poles", [])),
        "local_positive_poles": len(physical.get("positive_real_part_poles", [])),
        "a0_recorded": norm.get("recorded", ""),
        "a0_stored_roots_rel": norm.get("stored_relative_difference", ""),
        "a0_paired_roots_rel": norm.get("paired_relative_difference", ""),
        "sensitivity": (local.get("sensitivity") or {}).get("value", ""),
        "sensitivity_over_stage_product": local.get(
            "sensitivity_over_stage_product", ""
        ),
        "literature_table": literature.get("table", ""),
        "literature_match_as_stored": literature.get("poles_match_as_stored", ""),
        "literature_match_after_pairing": literature.get(
            "poles_match_after_pairing", ""
        ),
        "hydrophone_pole_rad_s": (hydrophone.get("pole_rad_s") or [""])[0],
        "hydrophone_R_C_total_s": hydrophone.get("implied_R_C_total_s", ""),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = parser.add_subparsers(dest="command", required=True)
    fetch_parser = sub.add_parser("fetch", help="download official StationXML")
    fetch_parser.add_argument("--out", required=True, help="new folder; must not exist")
    fetch_parser.add_argument("--stations", nargs="+", required=True)
    fetch_parser.add_argument("--channels", nargs="+", default=["BDH", "BHZ"])
    fetch_parser.add_argument("--timeout", type=float, default=60.0)
    compare_parser = sub.add_parser("compare", help="compare saved and official")
    compare_parser.add_argument("--workflow", required=True)
    compare_parser.add_argument("--workflow-sha256")
    compare_parser.add_argument("--official", required=True, help="fetch folder")
    compare_parser.add_argument(
        "--out", required=True, help="new folder; must not exist"
    )
    compare_parser.add_argument("--stations", nargs="+", required=True)
    args = parser.parse_args()
    if args.command == "fetch":
        run_fetch(args)
    else:
        run_compare(args)


if __name__ == "__main__":
    main()
