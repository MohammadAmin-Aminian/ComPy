"""Run ``r01_run_station.py`` over a list of stations, a few at a time.

Run: PYTHONPATH=<worktree> python r02_batch.py OUTDIR [--stations NET.STA ... | --tier 1|2|all]
     [--parallel 2] [--processes 6] [--tag auto] [extra arguments for r01 after --]

One log per station in OUTDIR/logs; a station whose JSON exists is skipped, so the
command can be repeated after an interruption. Failures are recorded in
OUTDIR/failed_<tag>.csv with the last line of the log and do not stop the run.
"""

import argparse
import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import archive  # noqa: E402


def one(sid, out, processes, tag, extra):
    net, sta = sid.split(".")
    done = Path(out) / net / sta / f"{sta}_{tag}.json"
    log = Path(out) / "logs" / f"{sid}_{tag}.log"
    if done.exists():
        return sid, "skipped", ""
    log.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(HERE / "r01_run_station.py"),
        sid,
        str(out),
        "--processes",
        str(processes),
        "--tag",
        tag,
        *extra,
    ]
    env = dict(
        os.environ,
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
        MPLBACKEND="Agg",
    )
    with open(log, "w") as fh:
        code = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env)
    lines = [ln for ln in log.read_text().splitlines() if ln.strip()]
    if code == 0 and done.exists():
        first = next((ln for ln in lines if ln.startswith(sid)), "")
        return sid, "done", first
    return sid, "failed", lines[-1] if lines else "no output"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out")
    parser.add_argument("--stations", nargs="*")
    parser.add_argument("--tier", default=None)
    parser.add_argument("--parallel", type=int, default=2)
    parser.add_argument("--processes", type=int, default=6)
    parser.add_argument("--tag", default="auto")
    args, extra = parser.parse_known_args()
    extra = [e for e in extra if e != "--"]
    plan = archive.plan()
    if args.stations:
        sids = args.stations
    else:
        sids = [
            s
            for s, r in plan.items()
            if args.tier in (None, "all") or r["tier"] == args.tier
        ]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    failed = []
    with ThreadPoolExecutor(args.parallel) as pool:
        for sid, status, line in pool.map(
            lambda s: one(s, out, args.processes, args.tag, extra), sids
        ):
            print(f"[{time.time() - t0:7.0f} s] {status:8s} {line or sid}", flush=True)
            if status == "failed":
                failed.append((sid, line))
    with open(out / f"failed_{args.tag}.csv", "w", newline="") as fh:
        csv.writer(fh).writerows([("sid", "last_log_line"), *failed])
    print(
        f"finished: {len(sids)} stations, {len(failed)} failed, {time.time() - t0:.0f} s",
        flush=True,
    )


if __name__ == "__main__":
    main()
