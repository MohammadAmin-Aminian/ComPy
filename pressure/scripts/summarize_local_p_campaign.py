"""Tabulate a completed local-P campaign beside saved Rayleigh/microseism values.

Read-only. Pairs rerun in a later campaign (--supersede) replace the base
campaign's row and are marked. Writes a CSV and a Markdown table; the
microseism value is the frozen-rule conditional g when one exists. Legacy
best-band/weekly medians are never promoted. Bands and channels differ between methods;
the table is context, not validation.
"""

import argparse
import csv
import json
from pathlib import Path

from compy_pressure.local_p_check import PROTOCOL_VERSION
from compy_pressure.local_p_io import fresh_file

ROOT = Path(__file__).resolve().parents[1]
RAYLEIGH = ROOT / "docs/calibration_network_summary.json"
FROZEN = ROOT / "docs/microseism_table_staging.csv"
BANDS = ["0.5-1", "1-2", "2-4", "4-8", "8-16"]


def markdown_rows(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) > 3 and "." in cells[0] and not set(cells[0]) <= {"-"}:
            rows.append(cells)
    return rows


def microseism_values():
    return {
        ".".join(seed.split(".")[:2]): (
            float(row["microseism_conditional_g"]),
            "frozen rule (conditional_unvalidated)",
        )
        for seed, row in microseism_records().items()
    }


def microseism_records():
    with FROZEN.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    result = {}
    for row in rows:
        if (
            row["calibration_check_pressure_vertical_microseisms_status"]
            != "conditional_repeatable_unvalidated"
            or not row["microseism_conditional_g"]
        ):
            continue
        if row["Seed_id"] in result:
            raise ValueError("Ambiguous microseism pressure epoch")
        result[row["Seed_id"]] = row
    return result


def rayleigh_values():
    data = json.loads(RAYLEIGH.read_text())
    return {
        f"{n['network']}.{s['station']}": (s["pressure_multiplier"], s["seed_id"])
        for n in data["networks"]
        for s in n["stations"]
    }


def pair_rows(campaign):
    report = json.loads((campaign / "campaign_summary.json").read_text())
    return {r["pressure_seed_id"]: r for r in report["rows"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--supersede", type=Path, nargs="*", default=[])
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    args.csv, args.markdown = fresh_file(args.csv), fresh_file(args.markdown)
    if args.csv == args.markdown:
        raise ValueError("CSV and Markdown outputs must differ")
    rows = pair_rows(args.campaign_dir)
    source = {k: args.campaign_dir.name for k in rows}
    for later in args.supersede:
        for key, row in pair_rows(later).items():
            rows[key] = row
            source[key] = later.name
    rayleigh, micro = rayleigh_values(), microseism_values()
    micro_records = {
        ".".join(k.split(".")[:2]): v for k, v in microseism_records().items()
    }
    table = []
    for key in sorted(rows, key=lambda k: (k.split(".")[0], k)):
        row = rows[key]
        summary = {}
        if row.get("result_dir"):
            path = Path(row["result_dir"]) / "summary.json"
            summary = json.loads(path.read_text()) if path.is_file() else {}
        station = ".".join(key.split(".")[:2])
        bands = {
            f"{b['band_hz'][0]:g}-{b['band_hz'][1]:g}": b
            for b in summary.get("bands", [])
            if b.get("events")
        }
        reasons = {
            k: v
            for k, v in (summary.get("event_reason_counts") or {}).items()
            if k not in {"accepted", "beyond_max_events_priority"}
        }
        ray = rayleigh.get(station, (None, ""))
        mic = micro.get(station, (None, ""))
        mic_record = micro_records.get(station, {})
        status = summary.get("outcome") or row.get("plan_status", "")
        legacy = bool(summary) and summary.get("protocol_version") != PROTOCOL_VERSION
        if legacy and status == "conditional_gain_candidate":
            status = "legacy_protocol_unvalidated"
        item = {
            "pressure_seed_id": key,
            "source_campaign": source[key],
            "status": status,
            "protocol_version": summary.get(
                "protocol_version", "legacy_v1" if summary else ""
            ),
            "legacy_noise_guard_unverified": legacy,
            "timing_mode": summary.get("timing_mode", ""),
            "water_depth_m": summary.get("water_depth_m", ""),
            "accepted": (summary.get("event_reason_counts") or {}).get("accepted", 0),
            "processed": summary.get("candidates_processed", ""),
            "regime": ",".join(summary.get("regime_counts", {}) or []),
            "polarity_majority": summary.get("polarity_majority", ""),
            "frequency_assessment": (summary.get("frequency_dependence") or {}).get(
                "assessment", ""
            ),
            "main_rejection": max(reasons, key=reasons.get) if reasons else "",
            "rayleigh_g": ray[0] if ray[0] is not None else "",
            "rayleigh_channel": ray[1],
            "microseism_g": mic[0] if mic[0] is not None else "",
            "microseism_basis": mic[1],
            "microseism_channel": mic_record.get("Seed_id", ""),
            "microseism_channel_start_time": mic_record.get("channel_start_time", ""),
            "microseism_channel_end_time": mic_record.get("channel_end_time", ""),
            "microseism_inventory_file": mic_record.get("inventory_file", ""),
            "same_pressure_channel_as_microseism": mic_record.get("Seed_id") == key
            if mic_record
            else "",
        }
        for band in BANDS:
            b = bands.get(band)
            item[f"g_{band}_hz"] = round(b["gain_median"], 3) if b else ""
            item[f"n_{band}_hz"] = b["events"] if b else ""
        table.append(item)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("x", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=list(table[0]) if table else ["pressure_seed_id"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(table)
    run = [t for t in table if t["timing_mode"] or t["processed"] != ""]
    lines = [
        "| Pressure channel | Status | Depth (m) | Accepted/processed | "
        + " | ".join(f"g {b} Hz (n)" for b in BANDS)
        + " | Polarity | Rayleigh g | Microseism g | Main rejection |",
        "|" + "---|" * (9 + len(BANDS) - 1),
    ]
    for t in run:
        cells = [
            t["pressure_seed_id"]
            + (" ‡" if t["source_campaign"] != args.campaign_dir.name else ""),
            t["status"]
            + (" (legacy guard)" if t["legacy_noise_guard_unverified"] else ""),
            f"{float(t['water_depth_m']):.0f}" if t["water_depth_m"] != "" else "",
            f"{t['accepted']}/{t['processed']}",
        ]
        for band in BANDS:
            g = t[f"g_{band}_hz"]
            cells.append(f"{g:.2f} ({t[f'n_{band}_hz']})" if g != "" else "")
        rg, mg = t["rayleigh_g"], t["microseism_g"]
        cells += [
            {1: "+1", -1: "−1"}.get(t["polarity_majority"], ""),
            f"{rg:.3f}"
            + ("†" if t["rayleigh_channel"] not in ("", t["pressure_seed_id"]) else "")
            if rg != ""
            else "",
            f"{mg:.3f}"
            + ("*" if "frozen" in t["microseism_basis"] else "")
            + ("†" if t["same_pressure_channel_as_microseism"] is False else "")
            if mg != ""
            else "",
            t["main_rejection"].replace("_", " "),
        ]
        lines.append("| " + " | ".join(cells) + " |")
    lines += [
        "",
        "* Frozen-rule model-conditional microseism estimate; "
        "not accepted calibration.",
        "† Different pressure channel. Frequencies also differ between methods.",
        "Rows labelled legacy guard have not been rerun with the v2 noise safeguard.",
    ]
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    with args.markdown.open("x") as stream:
        stream.write("\n".join(lines) + "\n")
    counts = {}
    for t in table:
        counts[t["status"]] = counts.get(t["status"], 0) + 1
    print(
        json.dumps(
            {"rows": len(table), "run": len(run), "status_counts": counts}, indent=2
        )
    )


if __name__ == "__main__":
    main()
