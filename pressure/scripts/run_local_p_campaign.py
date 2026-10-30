"""All-station local-P pressure/vertical check: plan, dry run, resumable run.

The plan is built from the waveform archive (channel folders and filenames),
the network StationXML and a fixed pairing rule, before any data are read:
each pressure channel with archive files is paired with the highest-rate
vertical (velocity, dip +-90) at the same location that also has files.
Pairs too slow for the local-P band, without full responses or without a
catalogue are listed with an explicit status and are not run.

Timing: per-channel decisions come from docs/local_p_timing.json, which records
the timing evidence (never gains): "drift_model" (lag predicted from a fitted
logger clock drift, refined by onsets), "search" (envelope lag search) or
"none". Unlisted channels use "none": zero lag with the independent
pressure/vertical onset check. A pressure sensor naming a separate APG logger
without a recorded decision is flagged `separate_logger_timing_unverified`.

Each pair runs in its own process into <campaign>/<NET.STA>/<P>__<Z>/ with a
log in <campaign>/logs/. A pair is skipped on resume when its summary is
complete and its fingerprint (code, options, inputs) matches; otherwise the
old folder is moved to <campaign>/attempts/. Inputs are read in place.
"""

import argparse
import csv
import fcntl
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from uuid import uuid4

import numpy as np
from obspy import read_inventory

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))

from compy_pressure.local_p_check import (  # noqa: E402
    DEFAULTS,
    DayFileIndex,
    file_state,
    sha256,
)
from compy_pressure.local_p_io import local_p_output  # noqa: E402
from compy_pressure.ppsd import atomic_json, response_handling  # noqa: E402

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

DATA = Path.home() / "Desktop/Data_Bruit_FM/Data"
INVENTORIES = Path.home() / "Desktop/Data_Bruit_FM/INVENTORIES"
MAIN = Path.home() / "Projects/bruit-fm"
CATALOGS = ROOT / "results/LOCAL_P/catalogs"
NETWORKS = {
    "YV": {
        "inventory": [INVENTORIES / "YV.xml"],
        "catalog": CATALOGS / "usgs_local_YV_20121231_20131127_m2p5.csv",
    },
    "XO": {
        "inventory": [INVENTORIES / "XO.xml"],
        "catalog": CATALOGS / "usgs_local_XO_20180520_20190101_m2p5.csv",
    },
    "7D": {
        # The INVENTORIES/7D.xml high-rate epochs have no response stages; the
        # Rayleigh 7D run fetched full responses for these stations.
        "inventory": [
            MAIN
            / "results/NETWORKS/7D_ALL_NO_TILT/_inputs/responses"
            / (
                "7D.J28A,J29A,J30A,J35A,J36A,J38A,J43A,J44A,J45A,"
                "J52A,J53A,J54A,J55A,J63A,J67A.xml"
            )
        ],
        "catalog": CATALOGS / "usgs_local_7D_20111019_20120101_m2p0.csv",
    },
    "YL": {"inventory": "per_station", "catalog": None},
}
CODE_FILES = [
    ROOT / "../compy_pressure" / name
    for name in (
        "local_p_check.py",
        "local_p_physics.py",
        "rayleigh_calibration.py",
        "archive.py",
        "ppsd.py",
        "local_p_io.py",
    )
] + [Path(__file__)]
TIMING_FILE = ROOT / "docs/local_p_timing.json"


def timing_options(seed_id, sensor, decisions):
    """Pipeline options and a review note for one pressure channel."""
    decision = decisions.get(seed_id)
    if decision is None:
        note = (
            "separate_logger_timing_unverified"
            if "APG_LOGGER" in sensor.upper()
            else ""
        )
        return {"timing_lag_search_s": 0.0, "timing_drift": None}, "none", note
    mode = decision["mode"]
    if mode == "drift_model":
        return {"timing_lag_search_s": 0.0, "timing_drift": decision["drift"]}, mode, ""
    if mode == "search":
        lag = float(decision.get("search_s", 150.0))
        return {"timing_lag_search_s": lag, "timing_drift": None}, mode, ""
    if mode == "none":
        return {"timing_lag_search_s": 0.0, "timing_drift": None}, mode, ""
    raise ValueError(f"Unknown timing mode for {seed_id}: {mode}")


MIN_BAND_HZ = DEFAULTS["bands_hz"][0][1]


def archive_channels(station_dir):
    """Channel codes that have at least one daily file (from names only)."""
    codes = set()
    for _, _, names in os.walk(station_dir):
        for name in names:
            parts = name.split(".")
            if len(parts) >= 5 and len(parts[3]) == 3:
                codes.add(parts[3])
    return codes


def inventory_channels(paths, network, station):
    found = []
    for path in paths:
        for net in read_inventory(str(path)).select(network=network, station=station):
            for sta in net:
                for cha in sta:
                    found.append((path, sta, cha))
    return found


def classify(cha):
    try:
        _, kind = response_handling(cha)
    except ValueError:
        return None
    units = cha.response.instrument_sensitivity.input_units.upper().replace(" ", "")
    if kind == "Pa^2/Hz":
        return "pressure"
    if (
        units == "M/S"
        and cha.dip is not None
        and np.isclose(abs(cha.dip), 90, atol=0.1)
    ):
        return "vertical"
    return None


def plan_network(network, data_root, stations=None, decisions=None):
    config = NETWORKS[network]
    decisions = decisions or {}
    rows = []
    for station_dir in sorted(p for p in (data_root / network).iterdir() if p.is_dir()):
        station = station_dir.name
        if stations and station not in stations:
            continue
        if config["inventory"] == "per_station":
            inventories = [station_dir / "inventory.xml"]
        else:
            inventories = config["inventory"]
        present = archive_channels(station_dir)
        channels = inventory_channels(inventories, network, station)
        kinds = {}
        for path, sta, cha in channels:
            kind = classify(cha)
            if kind and cha.code in present:
                key = (kind, cha.location_code, cha.code)
                kinds.setdefault(key, []).append((path, sta, cha))
        pressures = sorted(k for k in kinds if k[0] == "pressure")
        archive_pressure = sorted(c for c in present if c[1:] in {"DH", "XH"})
        covered = {k[2] for k in pressures}
        # Every archive pressure channel is listed, with or without a response.
        missing = [c for c in archive_pressure if c not in covered]
        for code in missing or ([] if pressures else ["?"]):
            rows.append(
                {
                    "network": network,
                    "station": station,
                    "pressure_seed_id": f"{network}.{station}..{code}",
                    "vertical_seed_id": "",
                    "plan_status": "no_full_pressure_response",
                    "note": "pressure files exist; no epoch with a full Pa response",
                }
            )
        for _, loc, code in pressures:
            entries = kinds[("pressure", loc, code)]
            p_rate = max(float(c.sample_rate) for _, _, c in entries)
            verticals = [k for k in kinds if k[0] == "vertical" and k[1] == loc]
            row = {
                "network": network,
                "station": station,
                "pressure_seed_id": f"{network}.{station}.{loc}.{code}",
                "inventory": str(entries[0][0]),
                "pressure_rate_hz": p_rate,
                "pressure_sensor": entries[0][2].sensor.description
                if entries[0][2].sensor
                else "",
            }
            if not verticals:
                rows.append(
                    {**row, "vertical_seed_id": "", "plan_status": "no_vertical_pair"}
                )
                continue
            best = max(
                verticals,
                key=lambda k: max(float(c.sample_rate) for _, _, c in kinds[k]),
            )
            v_rate = max(float(c.sample_rate) for _, _, c in kinds[best])
            row.update(
                {
                    "vertical_seed_id": f"{network}.{station}.{loc}.{best[2]}",
                    "vertical_rate_hz": v_rate,
                }
            )
            timing, mode, flag = timing_options(
                row["pressure_seed_id"], row["pressure_sensor"], decisions
            )
            row.update({"timing": timing, "timing_mode": mode})
            if flag:
                row["note"] = flag
            nyquist_band = min(p_rate, v_rate) / 2 * DEFAULTS["nyquist_fraction"]
            if nyquist_band < MIN_BAND_HZ:
                row["plan_status"] = "not_applicable_sample_rate"
                row["note"] = (
                    f"usable band <= {nyquist_band:.2f} Hz is below local-P frequencies"
                )
            elif not config["catalog"] or not Path(config["catalog"]).is_file():
                row["plan_status"] = "no_catalog"
            else:
                row["plan_status"] = "eligible"
                row["catalog"] = str(config["catalog"])
                row["waveform_dir"] = str(station_dir)
            rows.append(row)
    return rows


def fingerprint(row, options):
    manifests = {}
    for seed in (row["pressure_seed_id"], row["vertical_seed_id"]):
        index = DayFileIndex([row["waveform_dir"]], seed)
        paths = {p.resolve() for group in index.by_day.values() for p in group}
        manifests[seed] = [file_state(p) for p in sorted(paths)]
    content = {
        "pair": [row["pressure_seed_id"], row["vertical_seed_id"]],
        "code": {str(p.relative_to(ROOT)): sha256(p) for p in CODE_FILES},
        "options": options,
        "inventory": sha256(row["inventory"]),
        "catalog": sha256(row["catalog"]),
        "waveform_root": str(Path(row["waveform_dir"]).resolve()),
        "waveforms": manifests,
        "dependencies": {p: version(p) for p in ("numpy", "scipy", "obspy")},
    }
    return hashlib.sha256(json.dumps(content, sort_keys=True).encode()).hexdigest()[:16]


def run_pair(row, campaign, options, python):
    folder = (
        campaign
        / f"{row['network']}.{row['station']}"
        / (f"{row['pressure_seed_id']}__{row['vertical_seed_id']}")
    )
    log = (
        campaign / "logs" / f"{row['pressure_seed_id']}__{row['vertical_seed_id']}.log"
    )
    pair_options = {**options, **row["timing"]}
    digest = fingerprint(row, pair_options)
    summary_path = folder / "summary.json"
    if summary_path.is_file():
        try:
            previous = json.loads(summary_path.read_text())
        except (OSError, ValueError):
            previous = {}  # preserve corrupt attempts, never reuse them
        if previous.get("execution_status") in {"completed", "skipped"} and (
            previous.get("campaign_fingerprint") == digest
        ):
            return previous, folder, log, "resumed"
    if folder.exists():
        stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ") + uuid4().hex[:6]
        target = campaign / "attempts" / stamp / folder.relative_to(campaign)
        target.parent.mkdir(parents=True, exist_ok=True)
        folder.rename(target)
        if log.is_file():
            log.rename(target / f"launcher_{uuid4().hex[:8]}.log")
    options_file = campaign / "options" / f"{row['pressure_seed_id']}.json"
    options_file.parent.mkdir(parents=True, exist_ok=True)
    options_file.write_text(json.dumps(pair_options, indent=2) + "\n")
    command = [
        python,
        "-m",
        "compy_pressure.local_p_check",
        "--pressure",
        row["pressure_seed_id"],
        "--vertical",
        row["vertical_seed_id"],
        "--inventory",
        row["inventory"],
        "--waveform-dir",
        row["waveform_dir"],
        "--catalog",
        row["catalog"],
        "--results-dir",
        str(folder),
        "--options-json",
        str(options_file),
    ]
    log.parent.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "PYTHONPATH": str(ROOT.parent)}
    with log.open("w") as stream:
        stream.write(" ".join(command) + "\n\n")
        stream.flush()
        result = subprocess.run(
            command, stdout=stream, stderr=subprocess.STDOUT, env=env
        )
    if not summary_path.is_file():
        return (
            {
                "execution_status": "failed",
                "error": f"exit {result.returncode}; no summary",
            },
            folder,
            log,
            "ran",
        )
    summary = json.loads(summary_path.read_text())
    if result.returncode:
        summary.update(
            execution_status="failed", error=f"Child exited {result.returncode}"
        )
    summary["campaign_fingerprint"] = digest
    summary["campaign_dependencies_unchanged"] = (
        fingerprint(row, pair_options) == digest
    )
    summary["dependency_versions"] = {
        p: version(p) for p in ("numpy", "scipy", "obspy")
    }
    if not summary["campaign_dependencies_unchanged"]:
        summary.update(
            execution_status="failed",
            error="Campaign dependencies changed during pair run",
        )
    atomic_json(summary_path, summary)
    return summary, folder, log, "ran"


def summary_row(row, summary=None, folder=None, log=None, action=""):
    out = {
        k: row.get(k, "")
        for k in (
            "network",
            "station",
            "pressure_seed_id",
            "vertical_seed_id",
            "plan_status",
            "pressure_rate_hz",
            "vertical_rate_hz",
            "timing_mode",
            "note",
            "inventory",
        )
    }
    if summary is None:
        return out
    out.update(
        {
            "action": action,
            "execution_status": summary.get("execution_status", ""),
            "outcome": summary.get("outcome", ""),
            "error": summary.get("error", ""),
            "water_depth_m": summary.get("water_depth_m", ""),
            "candidates_processed": summary.get("candidates_processed", ""),
            "accepted_events": (summary.get("event_reason_counts") or {}).get(
                "accepted", 0
            ),
            "regime_counts": json.dumps(summary.get("regime_counts", {})),
            "frequency_assessment": (summary.get("frequency_dependence") or {}).get(
                "assessment", ""
            ),
            "polarity_majority": summary.get("polarity_majority", ""),
            "reasons": json.dumps(summary.get("event_reason_counts", {})),
            "result_dir": str(folder or ""),
            "log": str(log or ""),
        }
    )
    for band in summary.get("bands", []):
        tag = f"{band['band_hz'][0]:g}_{band['band_hz'][1]:g}hz"
        out[f"{tag}_events"] = band["events"]
        out[f"{tag}_gain_median"] = band.get("gain_median", "")
        out[f"{tag}_gain_mad"] = band.get("gain_mad", "")
    return out


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--networks", nargs="*", default=list(NETWORKS))
    parser.add_argument("--stations", nargs="*")
    parser.add_argument("--pressure-ids", nargs="*", help="run only these channels")
    parser.add_argument("--timing-file", type=Path, default=TIMING_FILE)
    parser.add_argument("--data-root", type=Path, default=DATA)
    parser.add_argument("--max-events", type=int, default=DEFAULTS["max_events"])
    parser.add_argument("--dry-run", action="store_true", help="print the plan only")
    args = parser.parse_args()
    campaign = local_p_output(args.campaign_dir)
    decisions = json.loads(args.timing_file.read_text())["decisions"]
    plan = []
    for network in args.networks:
        plan.extend(
            plan_network(network, args.data_root, set(args.stations or []), decisions)
        )
    if args.pressure_ids:
        plan = [r for r in plan if r["pressure_seed_id"] in set(args.pressure_ids)]
    counts = {}
    for row in plan:
        counts[row["plan_status"]] = counts.get(row["plan_status"], 0) + 1
    for row in plan:
        print(
            f"{row['plan_status']:28s} {row['pressure_seed_id']:18s} "
            f"{row.get('vertical_seed_id', ''):18s} "
            f"timing={row.get('timing_mode', '')} "
            f"{row.get('note', '')}"
        )
    print(json.dumps(counts, indent=2))
    if args.dry_run:
        return
    campaign.mkdir(parents=True, exist_ok=True)
    with (campaign / ".campaign.lock").open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("Another local-P campaign writer is active") from exc
        execute_plan(campaign, plan, args, counts)


def execute_plan(campaign, plan, args, counts):
    """Execute with the campaign lock held, isolating each pair's failures."""
    options = {"max_events": args.max_events}
    inputs = sorted(
        {r["inventory"] for r in plan if r.get("inventory")}
        | {r["catalog"] for r in plan if r.get("catalog")}
    )
    before = {p: sha256(p) for p in inputs}
    rows = []
    for row in plan:
        if row["plan_status"] != "eligible":
            rows.append(summary_row(row))
            continue
        print(
            f"running {row['pressure_seed_id']} / {row['vertical_seed_id']}", flush=True
        )
        try:
            summary, folder, log, action = run_pair(
                row, campaign, options, sys.executable
            )
        except Exception as exc:
            # Isolate parent-side IO, fingerprint and corrupt-summary failures.
            folder, log, action = None, None, "parent_failed"
            summary = {
                "execution_status": "failed",
                "error": f"{type(exc).__name__}: {exc}",
            }
        rows.append(summary_row(row, summary, folder, log, action))
        print(
            f"  -> {summary.get('execution_status')} {summary.get('outcome', '')}",
            flush=True,
        )
    after = {p: sha256(p) for p in inputs}
    fields = []
    for row in rows:
        fields.extend(k for k in row if k not in fields)
    with (campaign / "campaign_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "finished_utc": datetime.now(UTC).isoformat(),
        "plan_counts": counts,
        "execution_counts": {},
        "outcome_counts": {},
        "options": options,
        "inputs_unchanged": before == after,
        "inputs": before,
        "rows": rows,
    }
    for row in rows:
        for key, field in (
            ("execution_counts", "execution_status"),
            ("outcome_counts", "outcome"),
        ):
            value = row.get(field) or "not_run"
            report[key][value] = report[key].get(value, 0) + 1
    atomic_json(campaign / "campaign_summary.json", report)
    print(
        json.dumps(
            {
                k: report[k]
                for k in ("execution_counts", "outcome_counts", "inputs_unchanged")
            },
            indent=2,
        )
    )
    if before != after:
        raise SystemExit("Protected inputs changed during the campaign")


if __name__ == "__main__":
    main()
