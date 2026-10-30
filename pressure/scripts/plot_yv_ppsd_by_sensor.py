"""One pooled PPSD probability figure per YV pressure-sensor type.

Groups (StationXML sensor fields plus Stähler et al. 2016, Table 3/App. A):
INSU DPG; Geomar HTI-01-PCA hydrophones (paper: RR45, RR53, RR55, RR56);
all other DEPAS HTI hydrophones. Each group reuses
``compy_pressure.plot_network_ppsd.plot_network`` on a filtered copy of the master
(written in the output folder; the real master is only read).
"""

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path

from compy_pressure.plot_network_ppsd import plot_network
from compy_pressure.plot_ppsd import MASTER_CSV

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

GEOMAR = {"RR45", "RR53", "RR55", "RR56"}  # Stähler et al. (2016), App. A2


def sensor_group(row):
    station = row["Seed_id"].split(".")[1]
    if "DPG" in row.get("sensor_type", "").upper():
        return "INSU_DPG"
    if station in GEOMAR:
        return "Geomar_HTI-01-PCA_hydrophone"
    return "DEPAS_HTI_hydrophone"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--master-csv", default=str(MASTER_CSV))
    parser.add_argument("--channel", default="BDH")
    parser.add_argument("--output-dir")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    master = Path(args.master_csv).resolve()
    out = Path(
        args.output_dir
        or master.parent.parent
        / f"results/REVIEWS/PPSD_BY_SENSOR/run_{datetime.now(UTC):%Y%m%dT%H%M%SZ}"
    ).resolve()
    out.mkdir(parents=True, exist_ok=False)
    with master.open(newline="") as stream:
        reader = csv.DictReader(stream)
        fields, rows = reader.fieldnames, list(reader)
    groups = {}
    for row in rows:
        parts = row["Seed_id"].split(".")
        saved = row.get("ppsd_status") in {"complete", "partial"}
        if parts[0] == "YV" and parts[-1] == args.channel and saved:
            groups.setdefault(sensor_group(row), []).append(row)
    for name, members in sorted(groups.items()):
        subset = out / f"{name}_master.csv"
        with subset.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(members)
        stations = sorted({r["Seed_id"].split(".")[1] for r in members})
        print(f"{name}: {len(stations)} stations: {' '.join(stations)}")
        plot_network(
            subset,
            network="YV",
            channels=(args.channel,),
            output=out / name,
            probability=True,
            show=args.show,
        )
    print(f"Figures: {out}")


if __name__ == "__main__":
    main()
