"""Run stages B (clean), D (calibrate) and C (compliance) per station as soon
as stage A has finished that station, with at most --jobs stations at once.

A station counts as finished in stage A when every expected 6-day chunk has
written its QC file and no output of that station changed for SETTLE_S
seconds. Logs go to logs/<stage>_<STA>.log; a station's chain stops at the
first failing stage.
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

from a_preprocess import CHUNK_DAYS, DAY, station_days
from config import ROOT, STATIONS, WORK

PY = sys.executable
HERE = Path(__file__).resolve().parent
SETTLE_S = 300
ENV = {
    "MPLBACKEND": "Agg", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1", "MKL_NUM_THREADS": "1",
}  # fmt: skip


def expected_chunks(station):
    days = station_days(station)
    runs, run = [], []
    for day in days:
        if run and day - run[-1] > 1.5 * DAY:
            runs.append(run)
            run = []
        run.append(day)
    runs += [run] if run else []
    return [r[i] for r in runs for i in range(0, len(r), CHUNK_DAYS)]


def stage_a_done(station, chunks):
    if (WORK / station / "stageA.done").exists():  # written by a_fast.py
        return True
    qc = WORK / station / "qc"
    if not all(
        (qc / f"{station}.{d.year}.{d.julday:03d}.json").exists() for d in chunks
    ):
        return False
    files = [
        p for sub in ("phys", "counts", "qc") for p in (WORK / station / sub).glob("*")
    ]
    newest = max((p.stat().st_mtime for p in files), default=0)
    return time.time() - newest > SETTLE_S


def chain(station):
    steps = [("b_clean", ["b_clean.py", station]), ("d_calibrate", ["d_calibrate.py", station]),
             ("c_compliance", ["c_compliance.py", station])]  # fmt: skip
    (ROOT / "logs").mkdir(exist_ok=True)
    script = " && ".join(
        f"'{PY}' {' '.join(args)} > '{ROOT}/logs/{name}_{station}.log' 2>&1"
        for name, args in steps
    )
    env = dict(os.environ, **ENV)
    return subprocess.Popen(["/bin/bash", "-c", script], cwd=HERE, env=env)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=6)
    parser.add_argument("stations", nargs="*", default=list(STATIONS))
    args = parser.parse_args()
    chunks = {s: expected_chunks(s) for s in args.stations}
    waiting, running, finished = list(args.stations), {}, {}
    while waiting or running:
        for station, proc in list(running.items()):
            if proc.poll() is not None:
                finished[station] = proc.returncode
                del running[station]
                print(
                    f"{time.strftime('%H:%M:%S')} {station} chain exit {proc.returncode}",
                    flush=True,
                )
        for station in list(waiting):
            if len(running) >= args.jobs:
                break
            if stage_a_done(station, chunks[station]):
                waiting.remove(station)
                running[station] = chain(station)
                print(f"{time.strftime('%H:%M:%S')} {station} started", flush=True)
        time.sleep(30)
    print("BDC_DONE", finished, flush=True)


if __name__ == "__main__":
    main()
