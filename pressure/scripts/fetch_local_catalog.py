"""Fetch a USGS ComCat local-event catalogue (CSV) into a new file.

The query box is the station coordinates from a StationXML file plus a radius.
Existing files are never overwritten. A sidecar JSON records the exact query
URL, retrieval time, SHA-256 and row count. The CSV keeps ComCat's location,
depth and magnitude uncertainty columns.
"""

import argparse
import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from obspy import read_inventory

from compy_pressure.local_p_io import fresh_file

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

SERVICE = "https://earthquake.usgs.gov/fdsnws/event/1/query"


def station_box(inventory_path, stations, radius_km):
    inventory = read_inventory(str(inventory_path))
    coordinates = [
        (sta.latitude, sta.longitude)
        for net in inventory
        for sta in net
        if not stations or sta.code in stations
    ]
    if not coordinates:
        raise ValueError("No matching stations in the inventory")
    lat, lon = np.array(coordinates).T
    dlat = radius_km / 111.195
    dlon = radius_km / (111.195 * np.cos(np.radians(np.max(np.abs(lat)) + dlat)))
    return {
        "minlatitude": round(float(lat.min() - dlat), 3),
        "maxlatitude": round(float(lat.max() + dlat), 3),
        "minlongitude": round(float(lon.min() - dlon), 3),
        "maxlongitude": round(float(lon.max() + dlon), 3),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--stations", nargs="*", default=[])
    parser.add_argument("--radius-km", type=float, default=320.0)
    parser.add_argument("--start", required=True)
    parser.add_argument("--end", required=True)
    parser.add_argument("--min-magnitude", type=float, required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.suffix.lower() != ".csv":
        raise SystemExit("Catalogue output must end in .csv")
    output = fresh_file(args.output)
    sidecar_path = fresh_file(output.with_suffix(".json"))
    query = {
        "format": "csv",
        "starttime": args.start,
        "endtime": args.end,
        "minmagnitude": args.min_magnitude,
        "orderby": "time-asc",
        **station_box(args.inventory, set(args.stations), args.radius_km),
    }
    url = f"{SERVICE}?{urlencode(query)}"
    with urlopen(url, timeout=120) as response:
        content = response.read()
    rows = list(csv.DictReader(io.StringIO(content.decode())))
    output.parent.mkdir(parents=True, exist_ok=True)
    sidecar = {
        "url": url,
        "query": query,
        "retrieved_utc": datetime.now(UTC).isoformat(),
        "sha256": hashlib.sha256(content).hexdigest(),
        "rows": len(rows),
        "inventory": str(args.inventory.resolve()),
        "stations": args.stations,
        "radius_km": args.radius_km,
        "source": "USGS ComCat FDSN event service; contributing networks per row",
    }
    # No fixed sibling .tmp and no replacement of either existing destination.
    with TemporaryDirectory(prefix=".catalog-", dir=output.parent) as staging:
        import os

        staged = Path(staging) / "catalog.csv"
        staged.write_bytes(content)
        with sidecar_path.open("x") as stream:
            stream.write(json.dumps(sidecar, indent=2) + "\n")
        os.link(staged, output)  # atomic exclusive publication; never replace
    print(json.dumps({k: sidecar[k] for k in ("rows", "sha256", "url")}, indent=2))


if __name__ == "__main__":
    main()
