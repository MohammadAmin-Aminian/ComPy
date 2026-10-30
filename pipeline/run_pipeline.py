"""Run A → B → D → C → E → figures → report, with logs and a run manifest.

Set COMPY_RUN_DIR and COMPY_PIPELINE_CONFIG as in pipeline/README.md.
The default runs every configured station; positional station codes narrow it.
Existing stage checkpoints are reused. No scheduler or background queue is needed.
"""

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import config

HERE = Path(__file__).resolve().parent
STAGES = ("A", "B", "D", "C", "E", "F", "G")


def positive_int(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stations", nargs="*", default=list(config.STATIONS))
    parser.add_argument("--workers", type=positive_int, default=4)
    parser.add_argument("--processes", type=positive_int, default=4)
    parser.add_argument("--start-at", choices=STAGES, default="A")
    parser.add_argument("--stop-after", choices=STAGES, default="G")
    parser.add_argument(
        "--fetch", action="store_true", help="fetch inventory and catalogue first"
    )
    args = parser.parse_args(argv)
    if not args.stations:
        parser.error("select at least one station")
    first, last = STAGES.index(args.start_at), STAGES.index(args.stop_after)
    if first > last:
        parser.error("--start-at must precede --stop-after")
    env = dict(os.environ, MPLBACKEND="Agg", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
               COMPY_STATIONS=",".join(args.stations))  # fmt: skip
    # Imports of config in child stages must resolve paths exactly as in this process.
    env["COMPY_RUN_DIR"] = str(config.ROOT)
    if env.get("COMPY_PIPELINE_CONFIG"):
        env["COMPY_PIPELINE_CONFIG"] = str(Path(env["COMPY_PIPELINE_CONFIG"]).resolve())
    log_dir = config.ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).isoformat()
    manifest = log_dir / f"pipeline_{uuid4().hex[:12]}.json"
    record = dict(started_utc=now, status="running", stations=args.stations,
                  stages=list(STAGES[first:last + 1]), steps=[],
                  run_dir=str(config.ROOT), config=env.get("COMPY_PIPELINE_CONFIG"))  # fmt: skip

    def save():
        temporary = manifest.with_suffix(".tmp")
        temporary.write_text(json.dumps(record, indent=2) + "\n")
        temporary.replace(manifest)

    def run(stage, arguments, station=None):
        log = log_dir / f"{manifest.stem}_{stage}{'_' + station if station else ''}.log"
        step = dict(stage=stage, station=station, command=arguments, log=str(log))
        record["steps"].append(step)
        save()
        print(f"{stage} {station or ','.join(args.stations)} → {log}", flush=True)
        with log.open("w") as output:
            completed = subprocess.run([sys.executable, *arguments], cwd=HERE, env=env,
                                       stdout=output, stderr=subprocess.STDOUT)  # fmt: skip
        step["returncode"] = completed.returncode
        if completed.returncode:
            raise RuntimeError(
                f"stage {stage} failed (exit {completed.returncode}); see {log}"
            )
        if stage == "D":
            # A missing calibration is an explicit nominal/alternate-gain fallback in C/E.
            gain_file = config.WORK / station / "calibration" / "compy_gain.json"
            gain = json.loads(gain_file.read_text())
            if gain.get("error"):
                step["calibration_fallback"] = gain["error"]
                print(
                    f"{station}: calibration unavailable; C/E use their documented fallback",
                    flush=True,
                )
        if stage == "C":
            summary = config.RESULTS / station / f"{station}_summary.json"
            values = json.loads(summary.read_text())["results"].values()
            if not any(value.get("n_windows", 0) > 0 for value in values):
                raise RuntimeError(f"{station}: no compliance windows passed selection")
        save()

    save()
    try:
        if args.fetch:
            run("0_inventory", ["fetch_inventory.py"])
            run("0_catalog", ["fetch_catalog.py"])
        selected = STAGES[first : last + 1]
        if "A" in selected:
            run(
                "A", ["a_preprocess.py", *args.stations, "--workers", str(args.workers)]
            )
        for station in args.stations:
            for stage, script in (
                ("B", "b_clean.py"),
                ("D", "d_calibrate.py"),
                ("C", "c_compliance.py"),
                ("E", "e_invert.py"),
            ):
                if stage in selected:
                    command = [script, station]
                    if stage == "E":
                        command += ["--processes", str(args.processes)]
                    run(stage, command, station)
        if "F" in selected:
            run("F", ["f_figures.py", args.stations[0]])
        if "G" in selected:
            run("G", ["g_report.py"])
    except (OSError, ValueError, RuntimeError, KeyError, KeyboardInterrupt) as exc:
        record.update(status="failed", error=str(exc))
        save()
        print(f"Pipeline failed: {exc}\nRun manifest: {manifest}", file=sys.stderr)
        return 1
    record.update(
        status="complete", finished_utc=datetime.now(timezone.utc).isoformat()
    )
    save()
    print(f"Pipeline complete. Run manifest: {manifest}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
