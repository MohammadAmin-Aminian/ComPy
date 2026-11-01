"""Stage A driver for slow archives (USB hard disks).

Many workers reading day files in parallel make a spinning disk seek between
files, and throughput collapses. This driver copies one station's remaining
raw day files (plus one day either side, for the padding) to a local cache in
a single sequential stream, then runs ``a_preprocess.py`` on the cache with
``--workers`` processes. The copy of the next station overlaps with the
processing of the current one. When a station is finished, the driver writes
work/<STA>/stageA.done (read by run_bdc.py) and deletes its cache.

Days already processed are those listed in existing work/<STA>/qc/*.json.

Usage: python a_fast.py [--workers N] STATION [STATION ...]
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

import a_preprocess as A
from config import CHANNELS, NETWORK, ROOT, SDS, SDS_LOCATION, WORK

CACHE = ROOT / "cache_sds"
HERE = Path(__file__).resolve().parent


def processed_days(station):
    days = set()
    for path in (WORK / station / "qc").glob("*.json"):
        days.update(json.loads(path.read_text())["days"])
    return days


def remaining(station):
    done = processed_days(station)  # str(UTCDateTime), as written by stage A
    return [d for d in A.station_days(station) if str(d) not in done]


def copy_station(station, days):
    """Sequentially copy the day files of ``days`` (and their neighbours)."""
    from obspy import UTCDateTime

    wanted = [
        UTCDateTime(t)
        for t in sorted({(d + k * A.DAY).timestamp for d in days for k in (-1, 0, 1)})
    ]
    n_bytes = 0
    for channel in CHANNELS:
        for day in wanted:
            rel = (
                Path(str(day.year))
                / NETWORK
                / station
                / f"{channel}.D"
                / (
                    f"{NETWORK}.{station}.{SDS_LOCATION}.{channel}.D.{day.year}.{day.julday:03d}"
                )
            )
            src, dst = SDS / rel, CACHE / rel
            if not src.exists() or (
                dst.exists() and dst.stat().st_size == src.stat().st_size
            ):
                continue
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
            n_bytes += dst.stat().st_size
    return n_bytes


def process(station, workers):
    env = dict(os.environ, COMPY_SDS=str(CACHE), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1")  # fmt: skip
    log = open(ROOT / "logs" / f"a_fast_{station}.log", "w")
    return subprocess.run([sys.executable, "a_preprocess.py", station, "--workers", str(workers)],
                          cwd=HERE, env=env, stdout=log, stderr=subprocess.STDOUT).returncode  # fmt: skip


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stations", nargs="+")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    copied = {}

    def copier():
        for station in args.stations:
            t0 = time.time()
            days = remaining(station)
            n = copy_station(station, days) if days else 0
            copied[station] = (len(days), n, time.time() - t0)
            print(f"{time.strftime('%H:%M:%S')} {station} copied {len(days)} days, "
                  f"{n / 1e9:.1f} GB in {time.time() - t0:.0f} s", flush=True)  # fmt: skip

    thread = threading.Thread(target=copier, daemon=True)
    thread.start()
    for station in args.stations:
        while station not in copied:
            time.sleep(5)
        if copied[station][0]:
            code = process(station, args.workers)
        else:
            code = 0
        print(
            f"{time.strftime('%H:%M:%S')} {station} processed, exit {code}", flush=True
        )
        if code == 0:
            (WORK / station).mkdir(parents=True, exist_ok=True)
            (WORK / station / "stageA.done").write_text(
                time.strftime("%Y-%m-%dT%H:%M:%S")
            )
            for year in A.YEARS:
                path = CACHE / str(year) / NETWORK / station
                if path.is_dir():
                    shutil.rmtree(path)
    thread.join()
    print("A_FAST_DONE", flush=True)


if __name__ == "__main__":
    main()
