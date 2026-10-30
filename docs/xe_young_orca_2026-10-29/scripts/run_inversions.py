"""Run Stage E for every station whose compliance is ready, as results arrive.

A station is ready when results/<STA>/<STA>_summary.json exists and its
c_compliance log is complete; it is inverted once (variants final, main,
final_h2). Stations listed with --redo-dc first rerun stages D and C.
Stops when all stations in config.STATIONS are inverted or have failed, or
after --max-hours.
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

from config import RESULTS, ROOT, STATIONS

HERE = Path(__file__).resolve().parent
ENV = dict(os.environ, MPLBACKEND="Agg", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1")  # fmt: skip


def run(args, log):
    with open(ROOT / "logs" / log, "w") as fh:
        return subprocess.run([sys.executable, *args], cwd=HERE, env=ENV, stdout=fh,
                              stderr=subprocess.STDOUT).returncode  # fmt: skip


def ready(station):
    summary = RESULTS / station / f"{station}_summary.json"
    log = ROOT / "logs" / f"c_compliance_{station}.log"
    return (
        summary.exists() and log.exists() and time.time() - summary.stat().st_mtime > 30
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processes", type=int, default=6)
    parser.add_argument("--redo-dc", nargs="*", default=[])
    parser.add_argument("--max-hours", type=float, default=12)
    parser.add_argument("--variants", nargs="+", default=["final", "main", "final_h2"])
    parser.add_argument(
        "--reverse", action="store_true", help="take stations from the end"
    )
    args = parser.parse_args()
    for s in args.redo_dc:
        code = run(["d_calibrate.py", s], f"d_calibrate_{s}.log") or run(["c_compliance.py", s],
                                                                       f"c_compliance_{s}.log")  # fmt: skip
        print(f"{time.strftime('%H:%M:%S')} {s} redo D/C exit {code}", flush=True)
    done = {
        s
        for s in STATIONS
        if (RESULTS / s / "inversion" / f"{s}_{args.variants[0]}.json").exists()
    }
    t_end = time.time() + 3600 * args.max_hours
    order = list(reversed(STATIONS)) if args.reverse else list(STATIONS)
    while len(done) < len(STATIONS) and time.time() < t_end:
        todo = [s for s in order if s not in done and ready(s)]
        if not todo:
            time.sleep(60)
            continue
        s = todo[0]
        lock = RESULTS / s / "inversion" / f".{args.variants[0]}.lock"
        lock.parent.mkdir(parents=True, exist_ok=True)
        try:  # several queues may run; the lock gives each station to one
            lock.touch(exist_ok=False)
        except FileExistsError:
            done.add(s)
            continue
        code = run(["e_invert.py", s, *args.variants, "--processes", str(args.processes)],
                   f"e_invert_{s}.log")  # fmt: skip
        done.add(s)
        print(f"{time.strftime('%H:%M:%S')} {s} inversion exit {code}", flush=True)
    print("INVERSIONS_DONE", sorted(done), flush=True)


if __name__ == "__main__":
    main()
