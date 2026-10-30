"""Register five pilot stations or all available YL XMLs and cache the catalogue."""

import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

from obspy import UTCDateTime, read_inventory

from compy_pressure.archive import index_waveform_files
from compy_pressure.inventory import load_master_dataset, synchronize_inventory
from compy_pressure.network import scan_files

STATIONS = ("A02W", "A06W", "A08W", "A10W", "A12W")


def unchanged_complete_ppsd(row, files):
    """Skip completed spectra only when saved inventory and raw inputs still match."""
    npz = Path(row.get("ppsd_npz_file", ""))
    manifest = Path(row.get("ppsd_provenance_file", ""))
    if (
        row.get("ppsd_status") != "complete"
        or not npz.is_file()
        or not manifest.is_file()
    ):
        return False
    if not (npz.parent / "inventory.xml").is_file():
        return False
    saved = json.loads(manifest.read_text())
    source = [
        dict(path=str(p), size=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns)
        for p in sorted(map(Path, files.get(row.Seed_id, [])))
    ]
    unchanged = (
        saved["source_files"] == source
        and saved["config"]["inventory_sha256"]
        == hashlib.sha256(Path(row.inventory_file).read_bytes()).hexdigest()
    )
    if not unchanged:
        return False
    summary = saved.get("summary", {})
    if summary.get("ppsd_requested_start_time") != saved["config"]["origin"]:
        return False
    requested_end = summary.get("ppsd_requested_end_time")
    if not requested_end or summary.get("ppsd_status") != "complete":
        return False
    # 'Complete' describes the requested interval, which can be a short pilot.
    # Read headers only to establish the full available end for unchanged inputs.
    records = index_waveform_files([item["path"] for item in source])
    full_end = max(r.end + 1 / r.sample_rate for r in records)
    epoch_end = saved["config"].get("epoch_end")
    if epoch_end:
        full_end = min(full_end, UTCDateTime(epoch_end))
    return UTCDateTime(requested_end) == full_end


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path.home() / "Desktop/Data_Bruit_FM/Data",
    )
    parser.add_argument("--master-dir", type=Path, default=Path("exports"))
    parser.add_argument("--catalog-only", action="store_true")
    parser.add_argument(
        "--all-stations",
        action="store_true",
        help="Discover all station folders with inventory.xml below Data/YL",
    )
    args = parser.parse_args()
    root = args.data_dir / "YL"
    stations = (
        tuple(sorted(p.parent.name for p in root.glob("*/inventory.xml")))
        if args.all_stations
        else STATIONS
    )
    if not stations:
        raise FileNotFoundError(f"No YL station XMLs found below {root}")
    starts, ends, seeds = [], [], set()
    for station in stations:
        path = args.data_dir / "YL" / station / "inventory.xml"
        inv = read_inventory(str(path))
        channels = [
            c
            for n in inv
            if n.code == "YL"
            for s in n
            if s.code == station
            for c in s
            if c.code in {"LDH", "LHZ"}
        ]
        if {c.code for c in channels} != {"LDH", "LHZ"}:
            raise ValueError(f"Missing pressure/vertical metadata: {path}")
        if any(not c.response or not c.start_date or not c.end_date for c in channels):
            raise ValueError(f"Missing response/deployment bounds: {path}")
        starts.extend(c.start_date for c in channels)
        ends.extend(c.end_date for c in channels)
        seeds.update(f"YL.{station}.{c.location_code}.{c.code}" for c in channels)
    params = dict(
        format="geojson",
        starttime=str(min(starts)),
        endtime=str(max(ends) + 1),
        minmagnitude=5.5,
        orderby="time-asc",
    )
    url = "https://earthquake.usgs.gov/fdsnws/event/1/query?" + urlencode(params)
    prepared = Path("results/RAYLEIGH_REVIEW")
    if args.all_stations:
        dates = f"{min(starts).strftime('%Y%m%d')}_{(max(ends) + 1).strftime('%Y%m%d')}"
        catalog = prepared / f"catalog_usgs_m5p5_YL_ALL_{dates}.json"
    else:
        catalog = prepared / "catalog_usgs_m5p5_YL_2009_2010.json"
    if not catalog.exists():
        with urlopen(url, timeout=60) as response:
            content = response.read()
        doc = json.loads(content)
        if doc.get("type") != "FeatureCollection" or not doc.get("features"):
            raise ValueError("USGS returned no catalogue events")
        catalog.parent.mkdir(parents=True, exist_ok=True)
        with catalog.open("xb") as handle:
            handle.write(content)
        print(f"Saved {len(doc['features'])} events: {catalog}")
    else:
        print(f"Preserved existing catalogue: {catalog}")
    catalog_pointer = prepared / (
        "yl_all_catalog_path.txt" if args.all_stations else "yl_five_catalog_path.txt"
    )
    catalog_pointer.write_text(str(catalog) + "\n")
    seed_file = prepared / (
        "yl_all_pressure_vertical_seeds.txt"
        if args.all_stations
        else "yl_five_pressure_vertical_seeds.txt"
    )
    seed_text = "\n".join(sorted(seeds)) + "\n"
    if seed_file.exists() and seed_file.read_text() != seed_text:
        raise ValueError(
            f"Existing seed list differs; retain it before preparing again: {seed_file}"
        )
    if not seed_file.exists():
        seed_file.write_text(seed_text)
    print(f"Prepared {len(stations)} stations / {len(seeds)} P/Z channels: {seed_file}")
    if not args.catalog_only:
        if args.all_stations:
            print("YL", synchronize_inventory(root, args.master_dir))
        else:
            for station in stations:
                print(
                    station,
                    synchronize_inventory(root / station, args.master_dir),
                )
        files = scan_files(root, networks=["YL"], seed_ids=seeds)
        master = load_master_dataset(args.master_dir)
        needed = sorted(
            {
                row.Seed_id
                for _, row in master[master.Seed_id.isin(seeds)].iterrows()
                if not unchanged_complete_ppsd(row, files)
            }
        )
        plan = prepared / (
            "yl_all_ppsd_needed_seeds.txt"
            if args.all_stations
            else "yl_five_ppsd_needed_seeds.txt"
        )
        # This is a regenerated execution selector, not a scientific result.
        plan.write_text("\n".join(needed) + ("\n" if needed else ""))
        print(
            f"PPSD needed for {len(needed)} channels; "
            f"completed unchanged products retained: {plan}"
        )
    print(
        "Use original signed YL responses with explicit pressure-sign +1 for the pilot."
    )
    print("RR52 dispersion is a provisional borrowed model, not a YL reproduction.")


if __name__ == "__main__":
    main()
