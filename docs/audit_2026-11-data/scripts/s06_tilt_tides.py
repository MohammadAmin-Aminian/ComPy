"""Is the tilt noise driven by tidal currents?

Bottom currents tilt an OBS and add low-frequency horizontal noise that leaks
into the vertical; tidal currents modulate it at M2 (or M4 for the speed of a
purely tidal flow; Corela et al. 2023). The whole-archive run measured the tilt
every hour (``03_tilt/<NET>/<ID>_tilt_hourly.csv``: tilt angle, azimuth and the
variance reduction of the vertical by the tilt correction). Here, for each
inverted station, ``compy.anomaly.periodicity`` tests the hourly tilt angle and
variance reduction for tidal lines against non-tidal control periods.

Output: ``13_audit/tilt/tilt_tides.csv``.
"""

import csv
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from compy import anomaly  # noqa: E402

ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
OUT = Path(os.environ.get("COMPY_AUDIT_OUT", ARCHIVE / "13_audit")) / "tilt"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plan = list(
        csv.DictReader(
            open(ARCHIVE / "12_inversion_all_stations/tables/station_plan.csv")
        )
    )
    rows = []
    for p in plan:
        path = ARCHIVE / "03_tilt" / p["network"] / f"{p['sid']}_tilt_hourly.csv"
        if not path.exists():
            continue
        t, ang, vr = [], [], []
        for r in csv.DictReader(open(path)):
            try:
                t.append(
                    datetime.fromisoformat(
                        r["centre_utc"].replace("Z", "+00:00")
                    ).timestamp()
                )
                ang.append(float(r["angle_deg"]))
                vr.append(float(r["variance_reduction"]))
            except ValueError:
                continue
        t, ang, vr = map(np.asarray, (t, ang, vr))
        ok = np.isfinite(ang) & np.isfinite(vr)
        if ok.sum() < 24 * 30:
            continue
        row = {
            "sid": p["sid"],
            "network": p["network"],
            "tier": p["tier"],
            "n_hours": int(ok.sum()),
        }
        for name, x in (("angle", ang), ("vr", vr)):
            # remove slow drift (levelling, settling) with a 3-day running median first
            order = np.argsort(t[ok])
            tt, xx = t[ok][order], x[ok][order]
            from scipy.ndimage import median_filter

            resid = xx - median_filter(xx, size=73, mode="nearest")
            per = anomaly.periodicity(tt, resid)
            for c in ("M2", "M4", "K1", "S1"):
                v = per[c]
                row[f"{name}_{c}_amp"] = v["amplitude"]
                row[f"{name}_{c}_ratio"] = (
                    v["amplitude"] / v["null_median"]
                    if v["null_median"] > 0
                    else np.nan
                )
                row[f"{name}_{c}_p"] = v["p_tail"]
        rows.append(row)
    keys = sorted({k for r in rows for k in r})
    with open(OUT / "tilt_tides.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    sig = sum(1 for r in rows if r["angle_M2_p"] < 0.01)
    print(f"{len(rows)} stations; tilt angle modulated at M2 (tail p < 0.01) at {sig}")


if __name__ == "__main__":
    main()
