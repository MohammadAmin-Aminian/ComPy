"""Prepare isolated 7D P/Z processing using header-identified read-only symlinks."""

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import pandas as pd
from obspy import UTCDateTime, read

from compy_pressure.archive import waveform_seed_id
from compy_pressure.inventory import (
    EPOCH_COLUMNS,
    build_master_dataset,
    load_master_dataset,
    save_master_dataset,
)
from compy_pressure.ppsd import is_pressure_unit, response_handling, select_epoch

KEYS = EPOCH_COLUMNS + ["inventory_file"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def download(url, target, offline):
    """Fetch a URL for a new local file, or refuse when working offline."""
    if offline:
        raise SystemExit(
            f"--offline: local copy missing, refusing to download it: {target}"
        )
    with urlopen(url, timeout=60) as response:
        return response.read()


def source_header(path):
    filename_id = waveform_seed_id(path)
    traces = read(str(path), format="MSEED", headonly=True)
    ids = {t.id for t in traces}
    rates = {float(t.stats.sampling_rate) for t in traces}
    if len(ids) != 1 or len(rates) != 1:
        raise ValueError(f"Mixed identifiers or sample rates: {path}")
    seed = ids.pop()
    if filename_id is None:
        raise ValueError(f"Unrecognized filename: {path}")
    old, actual = filename_id.split("."), seed.split(".")
    if old[:2] != actual[:2] or old[3] != actual[3]:
        raise ValueError(f"Network/station/channel mismatch: {path} versus {seed}")
    # Only the observed 00/blank location discrepancy is authorized here.
    if old[2] != actual[2] and not (old[2] == "00" and actual[2] == ""):
        raise ValueError(f"Unexpected location mismatch: {path} versus {seed}")
    return dict(
        source=str(path.resolve()),
        filename_seed_id=filename_id,
        seed_id=seed,
        sample_rate_hz=rates.pop(),
        start_utc=str(min(t.stats.starttime for t in traces)),
        end_utc_exclusive=str(max(t.stats.endtime + t.stats.delta for t in traces)),
        size=path.stat().st_size,
        mtime_ns=path.stat().st_mtime_ns,
    )


def choose_pairs(master, available):
    """Prefer native LDH/LHZ at 1 Hz; retain exact epochs and equal native rates."""
    rows, plans = [], []
    stations = sorted({s.split(".")[1] for s in available})
    for station in stations:
        candidates = []
        scope = master[master.Seed_id.str.startswith(f"7D.{station}.")]
        for _, p in scope.iterrows():
            if not is_pressure_unit(p.input_units) or p.Seed_id not in available:
                continue
            if p.Seed_id.split(".")[-1] not in {"LDH", "BDH"}:
                continue
            prefix = p.Seed_id.rsplit(".", 1)[0]
            for _, z in scope.iterrows():
                if z.Seed_id.rsplit(".", 1)[0] != prefix or z.Seed_id not in available:
                    continue
                if z.Seed_id.split(".")[-1] not in {"LHZ", "BHZ"}:
                    continue
                if float(p.sample_rate_hz) != float(z.sample_rate_hz):
                    continue
                if not z.dip_deg or abs(abs(float(z.dip_deg)) - 90) > 0.1:
                    continue
                bounds = []
                for row in (p, z):
                    records = available[row.Seed_id]
                    if any(
                        r["sample_rate_hz"] != float(row.sample_rate_hz)
                        for r in records
                    ):
                        raise ValueError(f"Header/metadata rates differ: {row.Seed_id}")
                    left = min(UTCDateTime(r["start_utc"]) for r in records)
                    right = max(UTCDateTime(r["end_utc_exclusive"]) for r in records)
                    if row.channel_start_time:
                        left = max(left, UTCDateTime(row.channel_start_time))
                    if row.channel_end_time:
                        right = min(right, UTCDateTime(row.channel_end_time))
                    bounds.append((left, right))
                left, right = max(b[0] for b in bounds), min(b[1] for b in bounds)
                if right - left >= 10800:
                    candidates.append((float(p.sample_rate_hz), p, z, left, right))
        if not candidates:
            raise ValueError(f"No equal-rate, epoch-matched P/Z pair for {station}")
        rate = min(c[0] for c in candidates)
        choices = [c for c in candidates if c[0] == rate]
        if len(choices) != 1:
            raise ValueError(f"Ambiguous pressure/vertical epochs for {station}")
        _, p, z, left, right = choices[0]
        rows.extend([p, z])
        plans.append(
            dict(
                station=station,
                pressure=p.Seed_id,
                vertical=z.Seed_id,
                sample_rate_hz=rate,
                common_start_utc=str(left),
                common_end_utc=str(right),
                depth_m=float(p.Profondeur),
                sensor=p.sensor_description,
            )
        )
    return pd.DataFrame(rows).reset_index(drop=True), plans


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path.home() / "Desktop/Data_Bruit_FM/Data/7D",
    )
    parser.add_argument("--source-master", type=Path, default=Path("exports"))
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results/NETWORKS/7D_NO_TILT"),
    )
    parser.add_argument("--station", action="append")
    parser.add_argument("--low-rate-only", action="store_true")
    parser.add_argument("--fetch-responses", action="store_true")
    parser.add_argument("--exclude-incompatible-files", action="store_true")
    parser.add_argument(
        "--no-catalog", action="store_true", help="Preparation check only"
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Use only local response/catalogue copies already in --results-dir",
    )
    args = parser.parse_args()
    root = args.results_dir.resolve()
    protected = {str(p): sha(p) for p in args.source_master.glob("master_inventory.*")}
    master = load_master_dataset(args.source_master)
    master = master[master.Seed_id.str.startswith("7D.")].drop_duplicates(KEYS)
    if args.station:
        master = master[master.Seed_id.str.split(".").str[1].isin(args.station)]
    if master.empty:
        raise ValueError("No requested 7D metadata in source master")
    codes = {"LDH", "LHZ"} if args.low_rate_only else {"LDH", "LHZ", "BDH", "BHZ"}
    available = defaultdict(list)
    exclusions = []
    for path in sorted(args.data_dir.rglob("*")):
        seed = waveform_seed_id(path)
        if not path.is_file() or seed is None:
            continue
        parts = seed.split(".")
        if parts[0] != "7D" or parts[3] not in codes:
            continue
        if args.station and parts[1] not in args.station:
            continue
        try:
            item = source_header(path)
        except ValueError as exc:
            if not args.exclude_incompatible_files:
                raise
            exclusions.append(dict(source=str(path), reason=str(exc)))
            continue
        available[item["seed_id"]].append(item)
    if args.fetch_responses:
        root.mkdir(parents=True, exist_ok=True)
        response_dir = root / "_inputs/responses"
        response_dir.mkdir(parents=True, exist_ok=True)
        names = sorted({seed.split(".")[1] for seed in available})
        groups = names if args.station else [",".join(names)]
        for station in groups:
            target = response_dir / f"7D.{station}.xml"
            if not target.exists():
                records = [
                    r
                    for s, rs in available.items()
                    if s.split(".")[1] in station.split(",")
                    for r in rs
                ]
                query = dict(
                    network="7D",
                    station=station,
                    channel=",".join(sorted(codes)),
                    starttime=min(r["start_utc"] for r in records),
                    endtime=max(r["end_utc_exclusive"] for r in records),
                    level="response",
                    format="xml",
                )
                url = "https://service.iris.edu/fdsnws/station/1/query?" + urlencode(
                    query
                )
                content = download(url, target, args.offline)
                with target.open("xb") as handle:
                    handle.write(content)
        master = build_master_dataset(response_dir).astype(object).fillna("")
    if args.exclude_incompatible_files:
        for seed, items in list(available.items()):
            rates = set(
                master.loc[master.Seed_id.eq(seed), "sample_rate_hz"].map(float)
            )
            usable = []
            for item in items:
                if item["sample_rate_hz"] in rates:
                    usable.append(item)
                else:
                    exclusions.append(
                        dict(
                            source=item["source"],
                            reason=(
                                f"Native rate {item['sample_rate_hz']} "
                                f"differs from metadata {rates}"
                            ),
                        )
                    )
            if usable:
                available[seed] = usable
            else:
                del available[seed]
    selected, stations = choose_pairs(master, available)
    if selected.empty:
        raise ValueError("No usable downloaded P/Z pairs")
    for _, row in selected.iterrows():
        protected[row.inventory_file] = sha(row.inventory_file)
        _, channel = select_epoch(row)
        response_handling(channel)
    records = [r for seed in sorted(set(selected.Seed_id)) for r in available[seed]]
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "preparation.json"
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        if previous["source_files"] != records or previous["stations"] != stations:
            raise ValueError("Inputs/scope changed; use a fresh --results-dir")
        old_hashes = previous["source_inventory_sha256"]
        if any(sha(p) != h for p, h in old_hashes.items()):
            raise ValueError("Inventory changed; use a fresh --results-dir")
    view = root / "_inputs/waveforms"
    for item in records:
        source = Path(item["source"])
        tag = hashlib.sha256(str(source).encode()).hexdigest()[:16]
        stamp = UTCDateTime(item["start_utc"])
        link = view / tag / f"{item['seed_id']}.D.{stamp.year}.{stamp.julday:03d}"
        link.parent.mkdir(parents=True, exist_ok=True)
        if link.is_symlink():
            if link.resolve() != source:
                raise ValueError(f"Existing alias points elsewhere: {link}")
        elif link.exists():
            raise ValueError(f"Refusing to replace existing file: {link}")
        else:
            link.symlink_to(source)
    isolated = root / "_master"
    if not isolated.exists():
        save_master_dataset(selected, isolated)
    else:
        existing = load_master_dataset(isolated)
        if set(map(tuple, existing[KEYS].to_numpy())) != set(
            map(tuple, selected[KEYS].to_numpy())
        ):
            raise ValueError(
                "Isolated master scope differs; use a fresh results directory"
            )
    catalog = root / "_inputs/catalog_usgs_m5p5.json"
    if not args.no_catalog and not catalog.exists():
        query = dict(
            format="geojson",
            minmagnitude=5.5,
            orderby="time-asc",
            starttime=min(p["common_start_utc"] for p in stations),
            endtime=max(p["common_end_utc"] for p in stations),
        )
        url = "https://earthquake.usgs.gov/fdsnws/event/1/query?" + urlencode(query)
        content = download(url, catalog, args.offline)
        doc = json.loads(content)
        if doc.get("type") != "FeatureCollection" or not doc.get("features"):
            raise ValueError("USGS catalogue is empty or invalid")
        with catalog.open("xb") as handle:
            handle.write(content)
    for item in records:
        stat = Path(item["source"]).stat()
        if (stat.st_size, stat.st_mtime_ns) != (item["size"], item["mtime_ns"]):
            raise ValueError("Raw source changed during preparation")
    if any(sha(p) != h for p, h in protected.items()):
        raise ValueError("Original master/inventory changed during preparation")
    manifest = dict(
        stations=stations,
        source_files=records,
        source_inventory_sha256={
            p: h for p, h in protected.items() if p.endswith(".xml")
        },
        data_view=str(view),
        isolated_master=str(isolated),
        catalog=str(catalog),
        source_files_modified=False,
        excluded_files=exclusions,
        timestamps_modified=False,
        original_master_preserved=True,
        filename_view=(
            "Symlinks named by actual headers; "
            "only 00/blank location mismatches allowed"
        ),
    )
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Prepared {len(stations)} stations / {len(selected)} P/Z epochs")
    print(f"Excluded files: {len(exclusions)}; see preparation.json")
    print("Native rates:", {p["station"]: p["sample_rate_hz"] for p in stations})
    print(f"Isolated master: {isolated}\nHeader-based waveform view: {view}")
    if catalog.exists():
        count = len(json.loads(catalog.read_text())["features"])
        print(f"Cached catalogue events: {count}")


if __name__ == "__main__":
    main()
