"""``compy-pressure-atlas``: the global pressure-noise atlas from an SDS archive.

Subcommands:

``native``     native-rate spectra of every channel of a plan table, with one
               sequential reader of the archive (slow disks collapse under
               parallel reads) and a pool of workers; resumable.
``distances``  distance to the coastline and to the shelf break of every channel.
``build``      selection, screens, models and figures from the native-rate spectra.

Plan table for ``native`` (CSV): ``seed_id``, ``fs``, ``xml`` (StationXML with
the channel's response) and either ``files`` (day files separated by ``;``) or
nothing, in which case the files are found under ``--sds`` as
``<year>/<NET>/<STA>/<CHA>.D/*``. Channel table for ``build``: see
``compy_pressure.atlas.pipeline``.

Examples::

    compy-pressure-atlas native plan.csv --out native/ --cache /tmp/sds --workers 12
    compy-pressure-atlas distances channels.csv --out distances.csv
    compy-pressure-atlas build channels.csv --out atlas/ --config atlas.json --workers 10
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import sys
import threading
import time
import traceback
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


def _log(path, msg):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    if path:
        with open(path, "a") as fh:
            fh.write(line + "\n")


def _native_one(task):
    """Worker: native-rate spectra of one channel from its cached files, one output per response epoch, merged."""
    import warnings

    import numpy as np
    from obspy import read_inventory

    from . import native as nv

    warnings.filterwarnings("ignore")
    t0 = time.time()
    seed, cdir = task["seed_id"], Path(task["cache_dir"])
    res = dict(seed_id=seed, files=len(task["files"]))
    try:
        net, sta, loc, cha = seed.split(".")
        inv = read_inventory(task["xml"]).select(network=net, station=sta, channel=cha)
        epochs = sorted((c for n in inv for s in n for c in s if c.response is not None and c.response.response_stages),
                        key=lambda c: c.start_date)  # fmt: skip
        if not epochs:
            raise RuntimeError("no channel epoch with a response")
        local = [
            str(cdir / f"{i:06d}_{Path(f).name}") for i, f in enumerate(task["files"])
        ]
        day = {p: nv.file_day(p) for p in local}
        zs, parts = [], []
        for k, c in enumerate(epochs):
            a = c.start_date.timestamp / 86400 - 1
            b = (c.end_date.timestamp / 86400 + 1) if c.end_date else 1e9
            fl = [p for p in local if a <= day[p] <= b]
            if not fl:
                continue
            part = Path(task["out"] + f".epoch{k}.npz")
            nv.process_channel(fl, c.response, seed, task["fs"], part,
                               meta=dict(epoch_start=str(c.start_date), epoch_end=str(c.end_date)))  # fmt: skip
            parts.append(part)
            zs.append(nv.load(part))
        if not zs:
            raise RuntimeError("no files inside any response epoch")
        seen, rows = set(), []
        for k, z in enumerate(
            zs
        ):  # a window inside two epochs (overlap at the edges) is kept once, first epoch
            for i, s in enumerate(z["starts"]):
                if s not in seen:
                    seen.add(s)
                    rows.append((s, k, i))
        rows.sort()
        pick = lambda key: np.array([zs[k][key][i] for _, k, i in rows])  # noqa: E731
        np.savez_compressed(
            task["out"] + ".part.npz", periods=zs[0]["periods"], starts=np.array([r[0] for r in rows], float),
            levels_db=np.array([zs[k]["levels_db"][i] for _, k, i in rows], np.float32).reshape(-1, len(zs[0]["periods"])),
            maxabs_counts=pick("maxabs_counts"), std_counts=pick("std_counts"), long_context_ok=pick("long_context_ok"),
            epoch_index=np.array([r[1] for r in rows], np.int16), response_db=np.vstack([z["response_db"] for z in zs]),
            sensitivity=np.concatenate([z["sensitivity"] for z in zs]), input_units=zs[0]["input_units"], fs=task["fs"],
            seed_id=seed, product=zs[0]["product"], meta=json.dumps(dict(epochs=[z["meta"] for z in zs])))  # fmt: skip
        os.replace(task["out"] + ".part.npz", task["out"])
        for p in parts:
            p.unlink(missing_ok=True)
        errors = sum(z["meta"]["read_errors"] for z in zs)
        if errors:
            Path(task["out"]).unlink(missing_ok=True)
            raise RuntimeError(f"{errors} day files could not be read")
        res.update(status="ok", windows=len(rows), epochs=len(zs), units=zs[0]["input_units"],
                   tried=sum(z["meta"]["windows_tried"] for z in zs), gap=sum(z["meta"]["windows_gap"] for z in zs),
                   read_errors=sum(z["meta"]["read_errors"] for z in zs))  # fmt: skip
    except Exception as e:
        res.update(
            status="failed",
            error=f"{type(e).__name__}: {e}",
            trace=traceback.format_exc()[-1500:],
        )
    finally:
        shutil.rmtree(cdir, ignore_errors=True)
    res["seconds"] = round(time.time() - t0, 1)
    return res


def run_native(
    plan, out, cache, workers=8, max_cached_gb=150.0, max_ahead=16, sds=None, log=None
):
    """Native-rate spectra of every channel of ``plan`` (DataFrame), one reader and ``workers`` processes."""
    out, cache = Path(out), Path(cache)
    out.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    status = out / "status.jsonl"
    tasks = []
    for r in plan.itertuples():
        net, sta, _, cha = r.seed_id.split(".")
        if isinstance(getattr(r, "files", None), str) and r.files:
            files = r.files.split(";")
        else:
            files = (
                sorted(glob.glob(str(Path(sds) / "*" / net / sta / f"{cha}.D" / "*")))
                if sds
                else []
            )
        o = out / net / f"{r.seed_id}.npz"
        if files and not o.exists():
            tasks.append(dict(seed_id=r.seed_id, fs=float(r.fs), xml=r.xml, files=files, out=str(o),
                              bytes=sum(os.path.getsize(f) for f in files)))  # fmt: skip
    tasks.sort(key=lambda t: -t["bytes"])
    half = (len(tasks) + 1) // 2  # alternate large and small channels to bound memory
    tasks = [
        x for pair in zip(tasks[:half], tasks[::-1][: len(tasks) // 2]) for x in pair
    ] + ([tasks[len(tasks) // 2]] if len(tasks) % 2 else [])
    _log(
        log,
        f"native: {len(tasks)} channels to do, {sum(t['bytes'] for t in tasks) / 1e12:.2f} TB",
    )
    pending, lock, futs = {}, threading.Lock(), []
    with ProcessPoolExecutor(workers) as ex:
        for t in tasks:
            while True:
                with lock:
                    n, cached = len(pending), sum(pending.values())
                if n == 0 or (
                    n < max_ahead and (cached + t["bytes"]) / 1e9 < max_cached_gb
                ):
                    break
                time.sleep(5)
            if not all(Path(f).parent.exists() for f in t["files"][:1]):
                _log(
                    log,
                    f"archive not reachable ({t['files'][0]}); stopping, resume with the same command",
                )
                break
            cdir = cache / t["seed_id"]
            cdir.mkdir(parents=True, exist_ok=True)
            tc, copied = time.time(), 0
            for i, f in enumerate(t["files"]):
                dst = cdir / f"{i:06d}_{Path(f).name}"
                try:
                    size = os.path.getsize(f)
                    if not (
                        dst.exists() and dst.stat().st_size == size
                    ):  # a copy left by an interrupted run is reused
                        shutil.copyfile(f, dst)
                    copied += size
                except OSError as e:
                    _log(log, f"copy error {t['seed_id']} {f}: {e}")
            if copied < t["bytes"]:
                _log(log, f"incomplete copy of {t['seed_id']}; not processed")
                shutil.rmtree(cdir, ignore_errors=True)
                continue
            _log(
                log,
                f"copied {t['seed_id']} {copied / 1e9:.2f} GB in {time.time() - tc:.0f} s",
            )
            with lock:
                pending[t["seed_id"]] = copied

            def done(fu, seed=t["seed_id"]):
                try:
                    r = fu.result()
                except Exception as e:
                    r = dict(seed_id=seed, status="failed", error=f"pool: {e}")
                with lock:
                    pending.pop(seed, None)
                    with open(status, "a") as fh:
                        fh.write(json.dumps(r) + "\n")
                _log(
                    log,
                    f"done {seed}: {r.get('status')} {r.get('windows', '')} {r.get('error', '')}",
                )

            fu = ex.submit(_native_one, dict(t, cache_dir=str(cdir)))
            fu.add_done_callback(done)
            futs.append(fu)
        for fu in futs:
            fu.exception()


def main(argv=None):
    import pandas as pd

    from . import config as cfgmod

    ap = argparse.ArgumentParser(
        prog="compy-pressure-atlas",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("native", help="native-rate spectra from an SDS archive")
    a.add_argument("plan")
    a.add_argument("--out", required=True)
    a.add_argument(
        "--cache", required=True, help="local folder for the copied day files"
    )
    a.add_argument("--sds", help="SDS root, for plan rows without a file list")
    a.add_argument("--workers", type=int, default=8)
    a.add_argument("--max-cached-gb", type=float, default=150.0)
    a.add_argument("--log")
    b = sub.add_parser(
        "distances", help="distance to the coastline and to the shelf break"
    )
    b.add_argument("table")
    b.add_argument("--out", required=True)
    c = sub.add_parser("build", help="selection, screens, models and figures")
    c.add_argument("table")
    c.add_argument("--out", required=True)
    c.add_argument("--config")
    c.add_argument("--workers", type=int, default=8)
    c.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)
    for v in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "MKL_NUM_THREADS",
    ):
        os.environ.setdefault(v, "1")
    if args.cmd == "native":
        run_native(
            pd.read_csv(args.plan),
            args.out,
            args.cache,
            args.workers,
            args.max_cached_gb,
            sds=args.sds,
            log=args.log,
        )
    elif args.cmd == "distances":
        from . import geography

        t = pd.read_csv(args.table)
        d = geography.gauge_distances(t.lat, t.lon, t.depth_m)
        for k, v in d.items():
            t[k] = v.round(2)
        t.to_csv(args.out, index=False)
    else:
        from . import pipeline

        pipeline.build(
            args.table,
            args.out,
            cfgmod.load(args.config),
            args.workers,
            not args.no_figures,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
