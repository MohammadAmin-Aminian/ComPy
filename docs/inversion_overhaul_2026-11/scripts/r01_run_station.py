"""Invert one archive station with ``compy.inversion.auto`` and save the result.

Run: PYTHONPATH=<worktree> python r01_run_station.py NET.STA OUTDIR [--processes 8]
     [--response unchecked|corrected] [--quick] [--set key=value ...]

``--set`` overrides a ``SeafloorModel`` field (numbers), for sensitivity runs.
Output: OUTDIR/<NET>/<STA>/<STA>_auto.{npz,json} and _profile.csv.
"""

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import archive  # noqa: E402

import compy  # noqa: E402
from compy.inversion import auto  # noqa: E402


def run(
    sid,
    out,
    processes=8,
    response="auto",
    quick=False,
    overrides=None,
    tag="auto",
    max_fraction_fc=0.8,
    seed=0,
):
    data, column, info = archive.station_inputs(
        sid, response=response, max_fraction_fc=max_fraction_fc
    )
    model = auto.model_from_crust1(column, **(overrides or {}))
    sigma = auto.GAIN_SIGMA[info["gain_source"]]
    sampler = dict(n_iter=20000, n_burn=8000, max_rounds=1) if quick else {}
    t0 = time.time()
    res = auto.invert(
        data,
        model,
        gain_sigma=sigma,
        processes=processes,
        seed=seed,
        info=dict(info, compy=compy.__file__, version=compy.__version__),
        **sampler,
    )
    net, sta = sid.split(".")
    folder = Path(out) / net / sta
    folder.mkdir(parents=True, exist_ok=True)
    res.info["seconds"] = time.time() - t0
    res.save(folder / f"{sta}_{tag}")
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sid")
    parser.add_argument("out")
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--response", default="auto")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--tag", default="auto")
    parser.add_argument("--max-fraction-fc", type=float, default=0.8)
    parser.add_argument("--set", nargs="*", default=[])
    args = parser.parse_args()
    overrides = {}
    for item in args.set:
        key, value = item.split("=")
        overrides[key] = json.loads(value)
    res = run(
        args.sid,
        args.out,
        args.processes,
        args.response,
        args.quick,
        overrides,
        args.tag,
        args.max_fraction_fc,
    )
    s = res.summary()
    q = s["quality"]
    print(
        f"{args.sid} {args.tag}: {q['label']} chi2 {q['chi2']:.1f}/{q['n_data']} (limit {q['chi2_limit']:.1f}) "
        f"rhat {q['max_rhat']:.3f} ess {q['min_ess']:.0f} gain {q['gain']:.3f} rounds {s['rounds']} "
        f"reach {s['reach_m']:.0f} m resolved {s['resolved_from_m']}-{s['resolved_to_m']} m "
        f"{res.info['seconds']:.0f} s; {'; '.join(q['reasons'])}",
        flush=True,
    )
    for k, v in s["parameters"].items():
        p = s["prior"][k]
        print(
            f"  {k:20s} {v['p16']:9.3f} {v['median']:9.3f} {v['p84']:9.3f}   prior {p['p16']:9.3f} "
            f"{p['median']:9.3f} {p['p84']:9.3f}   rhat {v['rhat']:.3f}"
        )
