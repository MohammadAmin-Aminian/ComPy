"""Compose explicit saved campaigns and export a fresh, unscored phase review.

No waveform computation or master update. Later campaigns replace only their
explicit pressure IDs. Original summaries/spectra stay in place; historical
noise-guard results remain labelled legacy. No newest-directory guessing.
"""

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

from compy_pressure.local_p_check import PROTOCOL_VERSION, assert_fresh_output
from compy_pressure.ppsd import atomic_json

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def review(base, supersede, output, master):
    out = assert_fresh_output(output, [])
    if out.exists():
        raise FileExistsError("Phase review requires a new output directory")
    protected, rows, sources = {}, {}, {}

    def protect(path):
        path = Path(path).resolve()
        value = digest(path)
        if str(path) in protected and protected[str(path)] != value:
            raise RuntimeError(f"Protected input changed between reads: {path}")
        protected[str(path)] = value
        return path

    for path in [
        master,
        master.with_suffix(".xlsx"),
        ROOT / "docs/microseism_table_staging.csv",
        ROOT / "docs/calibration_network_summary.json",
    ]:
        protect(path)
    for campaign in [base, *supersede]:
        report = json.loads(protect(campaign / "campaign_summary.json").read_text())
        keys = [r["pressure_seed_id"] for r in report["rows"]]
        if len(keys) != len(set(keys)):
            raise ValueError("Duplicate pressure ID in campaign")
        for row in report["rows"]:
            seed = row["pressure_seed_id"]
            if campaign != base and seed not in rows:
                raise ValueError(
                    "Superseding campaign introduces an unplanned pressure ID"
                )
            rows[seed] = row
            sources[seed] = str(campaign.resolve())
            if row.get("result_dir"):
                folder = Path(row["result_dir"])
                summary = json.loads(protect(folder / "summary.json").read_text())
                for name in ["events.csv", "spectra.npz"]:
                    if (folder / name).is_file():
                        protect(folder / name)
                for name, state in summary.get("inputs", {}).items():
                    protect(name)
                    if protected[str(Path(name).resolve())] != state["sha256"]:
                        raise ValueError(f"Saved input hash changed: {name}")
    ordered = [rows[k] for k in sorted(rows)]
    out.mkdir(parents=True)
    atomic_json(
        out / "campaign_summary.json",
        {
            "derived_only": True,
            "sources": sources,
            "rows": ordered,
            "plan_counts": dict(Counter(r["plan_status"] for r in ordered)),
        },
    )
    env = {**os.environ, "PYTHONPATH": f"{ROOT.parent}:{ROOT}"}
    for script, args in [
        (
            "export_local_p_results.py",
            [
                "--campaign-dir",
                str(out),
                "--master",
                str(master),
                "--output-dir",
                str(out / "master_export"),
            ],
        ),
        (
            "summarize_local_p_campaign.py",
            [
                "--campaign-dir",
                str(out),
                "--csv",
                str(out / "station_table.csv"),
                "--markdown",
                str(out / "station_table.md"),
            ],
        ),
        (
            "apply_local_p_results.py",
            [
                "--export",
                str(out / "master_export/local_p_master_export.csv"),
                "--master-dir",
                str(master.parent),
            ],
        ),
    ]:
        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / script), *args], check=True, env=env
        )
    with (out / "station_table.csv").open(newline="") as stream:
        table = list(csv.DictReader(stream))
    if any(digest(path) != value for path, value in protected.items()):
        raise RuntimeError("Protected phase-review input changed")
    report = {
        "protocol_version": PROTOCOL_VERSION,
        "base_campaign": str(base.resolve()),
        "superseding_campaigns": [str(p.resolve()) for p in supersede],
        "pairs": len(table),
        "stations": len(
            {".".join(r["pressure_seed_id"].split(".")[:2]) for r in table}
        ),
        "outcomes": dict(Counter(r["status"] for r in table)),
        "legacy_protocol_rows": sum(
            r["legacy_noise_guard_unverified"] == "True" for r in table
        ),
        "protected_inputs_unchanged": True,
        "protected_sha256": protected,
        "master_updated": False,
        "score_assigned": False,
        "gain_applied": False,
    }
    atomic_json(out / "phase_report.json", report)
    print(
        json.dumps(
            {k: v for k, v in report.items() if k != "protected_sha256"}, indent=2
        )
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-campaign", type=Path, required=True)
    parser.add_argument("--supersede", type=Path, nargs="*", default=[])
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--master",
        type=Path,
        default=Path.home() / "Projects/bruit-fm/exports/master_inventory.csv",
    )
    args = parser.parse_args()
    review(args.base_campaign, args.supersede, args.output_dir, args.master)


if __name__ == "__main__":
    main()
