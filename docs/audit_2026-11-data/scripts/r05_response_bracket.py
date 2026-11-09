"""Leak corner from the tide alone, and the compliance bias it implies.

r04 shows that the tidal ratio R is flat across O1-S2 (median slope -0.11), as a
leak corner scaled by k predicts (R = 1/k far below the corner). So the tide
gives k directly: k_tide solves predicted_ratio(k) = R_M2 (leak-only model).

For each DPG this script compares k_tide with k_phase (s03) and computes:

* the band-bottom residual phase k_tide predicts vs the one measured;
* the compliance bias at the band bottom from k_tide and from k_phase (after
  the 20-60 mHz gain), which brackets the response bias.

Output: ``13_audit/research/response_bracket.csv``.
"""

import csv
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from compy import response_check as rc  # noqa: E402

ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
OUT = ARCHIVE / "13_audit/research"


def one(task):
    sid, cha, r_m2, k_phase, f, res = task
    from obspy import read_inventory

    net, sta = sid.split(".")
    inv = read_inventory(
        str(ARCHIVE / "00_inventory/stationxml" / f"{net}.xml")
    ).select(station=sta)
    chans = [c for n in inv for s in n for c in s if c.code == cha]
    if not chans:
        return {"sid": sid}
    stages = rc.paz_stages(chans[0])
    f = np.asarray(f)

    def g(lk):
        return np.log(rc.predicted_ratio(stages, np.exp(lk))) - np.log(r_m2)

    try:
        k_tide = float(np.exp(brentq(g, np.log(1e-3), np.log(1e3))))
    except ValueError:
        return {"sid": sid, "k_phase": k_phase, "R_M2": r_m2}
    pred = rc.predicted_residual(stages, f, k_tide)
    return {"sid": sid, "R_M2": r_m2, "k_phase": k_phase, "k_tide": k_tide,
            "res_lo_measured": float(res[0]), "res_lo_from_tide": float(pred[0]),
            "bias_lo_tide_pct": float(rc.compliance_bias(stages, f[:1], k_tide)[0]),
            "bias_lo_phase_pct": float(rc.compliance_bias(stages, f[:1], k_phase)[0]),
            "f_lo_mhz": 1e3 * float(f[0])}  # fmt: skip


def main():
    tide = {
        r["sid"]: r
        for r in csv.DictReader(open(ARCHIVE / "13_audit/dpg_tide/dpg_tide.csv"))
    }
    phase = {
        r["sid"]: r
        for r in csv.DictReader(open(ARCHIVE / "13_audit/phase/phase_metadata.csv"))
    }
    curves = json.loads((ARCHIVE / "13_audit/phase/phase_curves.json").read_text())
    tasks = []
    for sid, t in tide.items():
        if (
            t.get("tide_status") != "ok"
            or t["pressure_type"] != "DPG"
            or float(t["snr"]) < 20
        ):
            continue
        tasks.append((sid, phase[sid]["p_channel"], float(t["R_M2"]), float(t["k_lowest"]),
                      curves[sid]["f"], curves[sid]["residual"]))  # fmt: skip
    with ProcessPoolExecutor(int(os.environ.get("WORKERS", 8))) as ex:
        out = list(ex.map(one, tasks))
    for o in out:
        o["network"] = o["sid"].split(".")[0]
        o["family"] = tide[o["sid"]]["leak_corners_mhz"]
    keys = sorted({k for o in out for k in o})
    with open(OUT / "response_bracket.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(out)
    ok = [o for o in out if "k_tide" in o]
    print(f"{len(ok)} DPGs; median k_tide {np.median([o['k_tide'] for o in ok]):.2f}, "
          f"k_phase {np.median([o['k_phase'] for o in ok]):.2f}; bias tide "
          f"{np.median([o['bias_lo_tide_pct'] for o in ok]):+.1f} %, phase {np.median([o['bias_lo_phase_pct'] for o in ok]):+.1f} %")  # fmt: skip


if __name__ == "__main__":
    main()
