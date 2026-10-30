"""Run one baseline station's bounded pressure/vertical microseism check."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
from uuid import uuid4

from compy_pressure.inventory import load_master_dataset
from compy_pressure.microseism_diagnostic import run_diagnostic
from compy_pressure.paired_spectra import KEYS, run_paired_spectra

ROOT = Path(__file__).resolve().parents[1] / "pressure"
DEFAULT_MANIFEST = (
    Path(__file__).with_name("reference_data") / "microseism_baseline_stations.csv"
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_manifest(path=DEFAULT_MANIFEST):
    """Read and validate the frozen baseline-station comparison manifest."""
    path = Path(path)
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "network",
        "station",
        "pressure_seed_id",
        "vertical_seed_id",
        "fixed_pressure_multiplier",
        "comparison_direction",
        "master_dir",
        "data_dir",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Microseism manifest is empty or lacks required columns")
    names = [f"{row['network']}.{row['station']}" for row in rows]
    if len(names) != len(set(names)):
        raise ValueError("Microseism manifest has duplicate network/station rows")
    for row, name in zip(rows, names, strict=True):
        gain = float(row["fixed_pressure_multiplier"])
        if gain <= 0:
            raise ValueError(f"Invalid fixed multiplier for {name}")
        if row["comparison_direction"] not in {"g", "1/g"}:
            raise ValueError(f"Invalid comparison direction for {name}")
        if row["pressure_seed_id"].split(".")[:2] != name.split("."):
            raise ValueError(f"Pressure channel does not match {name}")
        if row["vertical_seed_id"].split(".")[:2] != name.split("."):
            raise ValueError(f"Vertical channel does not match {name}")
    return rows


def resolve_project_path(value):
    """Resolve manifest paths without changing their recorded text."""
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def select_station(rows, station):
    """Require a fully qualified and unique NETWORK.STATION selector."""
    matches = [row for row in rows if f"{row['network']}.{row['station']}" == station]
    if len(matches) != 1:
        raise ValueError(f"Unknown baseline station {station!r}; use --list")
    return matches[0]


def exact_pressure_rows(master, row):
    """Select only the manifest's exact pressure ID and refuse epoch ambiguity."""
    selected = master[master.Seed_id.eq(row["pressure_seed_id"])].drop_duplicates(KEYS)
    if len(selected) != 1:
        raise ValueError(
            f"Expected one epoch for {row['pressure_seed_id']}, found {len(selected)}"
        )
    vertical = master[master.Seed_id.eq(row["vertical_seed_id"])].drop_duplicates(KEYS)
    if len(vertical) != 1:
        raise ValueError(
            f"Expected one epoch for {row['vertical_seed_id']}, found {len(vertical)}"
        )
    return selected


def validate_comparison_binding(path, row):
    """Refuse to display one station's frozen gain on another P/Z product."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Paired comparison does not exist: {path}")
    summary = json.loads(path.read_text())
    pressure = summary.get("pressure_epoch", {}).get("Seed_id")
    vertical = summary.get("vertical_seed_id")
    if pressure != row["pressure_seed_id"] or vertical != row["vertical_seed_id"]:
        raise ValueError(
            "Paired comparison channels do not match the station manifest: "
            f"{pressure} / {vertical}"
        )
    return summary


def comparison_context(row, manifest_path):
    """Record the frozen first-method and literature references beside the result."""
    gain = float(row["fixed_pressure_multiplier"])
    direction = row["comparison_direction"]
    coefficient = gain if direction == "g" else 1 / gain
    published = row.get("published_coefficient", "")
    published_value = float(published) if published else None
    return {
        "manifest": {
            "path": str(Path(manifest_path).resolve()),
            "sha256": sha256(manifest_path),
        },
        "first_method": {
            "status": "fixed_provisional_reference",
            "pressure_multiplier_g": gain,
            "comparison_direction": direction,
            "comparison_coefficient": coefficient,
            "source": "docs/socrates_all_station_comparison.csv: previous_g",
        },
        "publication": {
            "reference": row.get("published_reference", ""),
            "coefficient": published_value,
            "status": row.get("published_comparison_status", ""),
            "first_method_difference_percent": (
                100 * (coefficient / published_value - 1)
                if published_value is not None
                else None
            ),
        },
        "microseism_gain_estimated": False,
    }


def run_station(
    station,
    *,
    manifest_path=DEFAULT_MANIFEST,
    start=None,
    days=7,
    comparison=None,
    fresh_pair=False,
    master_dir=None,
    data_dir=None,
):
    """Run one station only; never update a master, score, response, or waveform."""
    rows = load_manifest(manifest_path)
    row = select_station(rows, station)
    fixed_gain = float(row["fixed_pressure_multiplier"])
    saved = row.get("saved_comparison", "")
    if comparison is not None:
        comparison_path = Path(comparison).expanduser().resolve()
    elif saved and not fresh_pair:
        comparison_path = resolve_project_path(saved)
    else:
        master_path = (
            Path(master_dir).expanduser().resolve()
            if master_dir is not None
            else resolve_project_path(row["master_dir"])
        )
        data_path = (
            Path(data_dir).expanduser().resolve()
            if data_dir is not None
            else resolve_project_path(row["data_dir"])
        )
        station_data = data_path / row["network"] / row["station"]
        if station_data.is_dir():
            data_path = station_data
        master = load_master_dataset(master_path)
        selected = exact_pressure_rows(master, row)
        paired_output = (
            ROOT
            / "results/DIAGNOSTICS/MICROSEISM_PAIRED"
            / row["network"]
            / row["station"]
            / f"run_{uuid4().hex[:12]}"
        )
        resample = row.get("event_resample_rate_hz", "")
        tolerance = row.get("native_rate_tolerance_ppm", "")
        paired = run_paired_spectra(
            master_path,
            data_dir=data_path,
            network=row["network"],
            station=row["station"],
            start=start,
            days=days,
            frequency_max=0.4,
            vertical_id=row["vertical_seed_id"],
            pressure_id=row["pressure_seed_id"],
            output_directory=paired_output,
            event_resample_rate=float(resample) if resample else None,
            native_rate_tolerance_ppm=float(tolerance) if tolerance else 0,
            update_master=False,
            selected_rows=selected,
            reference_rows=master,
        )
        comparison_path = Path(paired["files"]["result"])
        if paired["status"] != "comparison_only":
            raise RuntimeError(
                f"{station} pairing ended as {paired['status']}: {comparison_path}"
            )

    validate_comparison_binding(comparison_path, row)
    diagnostic_output = (
        ROOT
        / "results/DIAGNOSTICS/MICROSEISM_BASELINE"
        / row["network"]
        / row["station"]
        / f"run_{uuid4().hex[:12]}"
    )
    result = run_diagnostic(
        comparison_path,
        diagnostic_output,
        pressure_multiplier=fixed_gain,
        comparison_context=comparison_context(row, manifest_path),
    )
    return {
        "station": station,
        "pressure_seed_id": row["pressure_seed_id"],
        "vertical_seed_id": row["vertical_seed_id"],
        "fixed_pressure_multiplier": fixed_gain,
        "paired_comparison": str(comparison_path),
        "diagnostic": str(result),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--station", help="Fully qualified NETWORK.STATION")
    parser.add_argument("--list", action="store_true", help="List supported stations")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--start", help="Optional UTC start; otherwise choose midpoint")
    parser.add_argument("--days", type=int, default=7)
    parser.add_argument("--comparison", type=Path, help="Reuse an explicit saved run")
    parser.add_argument(
        "--fresh-pair",
        action="store_true",
        help="Create fresh paired spectra even where a designated saved run exists",
    )
    parser.add_argument("--master-dir", type=Path, help="Explicit master override")
    parser.add_argument("--data-dir", type=Path, help="Explicit waveform-root override")
    args = parser.parse_args()
    rows = load_manifest(args.manifest)
    if args.list:
        for row in rows:
            print(
                f"{row['network']}.{row['station']}  "
                f"{row['pressure_seed_id']} / {row['vertical_seed_id']}  "
                f"g={row['fixed_pressure_multiplier']}"
            )
        return
    if not args.station:
        parser.error("--station is required unless --list is used")
    result = run_station(
        args.station,
        manifest_path=args.manifest,
        start=args.start,
        days=args.days,
        comparison=args.comparison,
        fresh_pair=args.fresh_pair,
        master_dir=args.master_dir,
        data_dir=args.data_dir,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
