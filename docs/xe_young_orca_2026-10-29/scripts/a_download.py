"""Stage A driver fetching raw day files from an FDSN dataselect service.

Same role as a_fast.py, for stations whose local archive is too slow to read:
the remaining days (those without a stage-A QC record, plus one day either
side for padding) are requested per channel and day from the FDSN service
with ``--connections`` parallel requests, written into the local cache in SDS
layout, decimated with ``a_preprocess.py``, and the station is marked done.
Downloading the next station overlaps with processing the current one.

The day list comes from the local SDS archive (directory listing only), so the
same days are processed whichever source supplies them.

Usage: python a_download.py [--connections 6] [--workers 8] STATION [STATION ...]
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import a_preprocess as A
from config import CHANNELS, NETWORK, ROOT, SDS_LOCATION, WORK

DATASELECT = "https://service.earthscope.org/fdsnws/dataselect/1/query"
CACHE = ROOT / "cache_sds"
HERE = Path(__file__).resolve().parent


def remaining(station):
    done = set()
    for path in (WORK / station / "qc").glob("*.json"):
        done.update(json.loads(path.read_text())["days"])
    return [d for d in A.station_days(station) if str(d) not in done]


def fetch(station, channel, day):
    rel = (
        Path(str(day.year))
        / NETWORK
        / station
        / f"{channel}.D"
        / (
            f"{NETWORK}.{station}.{SDS_LOCATION}.{channel}.D.{day.year}.{day.julday:03d}"
        )
    )
    dst = CACHE / rel
    if dst.exists() and dst.stat().st_size > 0:
        return 0
    url = (f"{DATASELECT}?net={NETWORK}&sta={station}&loc=*&cha={channel}"
           f"&start={(day - 180).strftime('%Y-%m-%dT%H:%M:%S')}"
           f"&end={(day + 86400 + 180).strftime('%Y-%m-%dT%H:%M:%S')}")  # fmt: skip
    dst.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=300) as response:
                if response.status == 204:
                    return 0
                tmp = dst.with_name(dst.name + ".part")  # unique per day and channel
                with open(tmp, "wb") as fh:
                    shutil.copyfileobj(response, fh, 1 << 20)
                tmp.rename(dst)
                return dst.stat().st_size
        except Exception as exc:  # retry transient network errors
            error = exc
            time.sleep(10 * (attempt + 1))
    print(f"FETCH-FAIL {station} {channel} {day.date}: {error!r}", flush=True)
    return 0


def download_station(station, days, connections):
    import obspy

    wanted = sorted({(d + k * A.DAY).timestamp for d in days for k in (-1, 0, 1)})
    jobs = [(station, c, obspy.UTCDateTime(t)) for c in CHANNELS for t in wanted]
    with ThreadPoolExecutor(connections) as pool:
        return sum(pool.map(lambda job: fetch(*job), jobs))


def process(station, workers):
    env = dict(os.environ, COMPY_SDS=str(CACHE), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1")  # fmt: skip
    with open(ROOT / "logs" / f"a_download_{station}.log", "w") as log:
        return subprocess.run([sys.executable, "a_preprocess.py", station, "--workers", str(workers)],
                              cwd=HERE, env=env, stdout=log, stderr=subprocess.STDOUT).returncode  # fmt: skip


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stations", nargs="+")
    parser.add_argument("--connections", type=int, default=6)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    fetched = {}

    def downloader():
        for station in args.stations:
            t0 = time.time()
            days = remaining(station)
            n = download_station(station, days, args.connections) if days else 0
            fetched[station] = len(days)
            print(f"{time.strftime('%H:%M:%S')} {station} downloaded {len(days)} days, "
                  f"{n / 1e9:.1f} GB in {time.time() - t0:.0f} s", flush=True)  # fmt: skip

    thread = threading.Thread(target=downloader, daemon=True)
    thread.start()
    for station in args.stations:
        while station not in fetched:
            time.sleep(5)
        code = process(station, args.workers) if fetched[station] else 0
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
    print("A_DOWNLOAD_DONE", flush=True)


if __name__ == "__main__":
    main()
