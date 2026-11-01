"""Stage H: temporal stability of the Z/P transfer function per station.

For every accepted window (Stage C, 'bruitfm' gain) the band median of
|S_pz|/S_pp and the band-mean Z and P levels are saved. The largest single
step (median before vs after, at least MIN_SIDE windows each side) is
located; a step larger than STEP_LIMIT in log10 (≈ 10 %) is reported with
the per-event Rayleigh gains before and after it (Stage D), which test the
same sensitivity ratio independently.
"""

import json
import sys
from pathlib import Path

import numpy as np
from obspy import UTCDateTime

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "ComPy"))
sys.path.insert(0, str(ROOT / "scripts"))

import c_compliance as C  # noqa: E402
from compy_processing import spectra  # noqa: E402

STATIONS = ("RR28", "RR29", "RR34", "RR36", "RR38", "RR40", "RR50", "RR52")
MIN_SIDE = 100
STEP_LIMIT = 0.04


def main():
    report = {}
    for station in STATIONS:
        d = dict(
            np.load(ROOT / "results" / station / f"{station}_compliance_bruitfm.npz")
        )
        band, accepted = d["band"], set(d["window_starts"].tolist())
        windows, _ = C.load_windows(station)
        rows = []
        for w in windows:
            t = float(w[0].stats.starttime)
            if t not in accepted:
                continue
            f, zz, pp, zp, _ = spectra(w, C.NSEG)
            sel = (f >= band[0]) & (f <= band[1])
            rows.append((t, float(np.log10(np.median(np.abs(zp[sel]) / pp[sel]))),
                         float(10 * np.log10(np.mean(pp[sel]))),
                         float(10 * np.log10(np.mean(zz[sel] * (2 * np.pi * f[sel]) ** 4)))))  # fmt: skip
        rows.sort()
        arr = np.array(rows)
        np.save(ROOT / "results" / station / f"{station}_window_levels.npy", arr)
        r = arr[:, 1]
        best = max(range(MIN_SIDE, len(r) - MIN_SIDE),
                   key=lambda i: abs(np.median(r[:i]) - np.median(r[i:])))  # fmt: skip
        step = float(np.median(r[best:]) - np.median(r[:best]))
        when = UTCDateTime(arr[best, 0])
        cal = json.loads(
            (ROOT / "work" / station / "calibration" / "compy_gain.json").read_text()
        )
        before = [
            g
            for t, g in zip(cal["events_used"], cal["event_gains"])
            if UTCDateTime(t) < when
        ]
        after = [
            g
            for t, g in zip(cal["events_used"], cal["event_gains"])
            if UTCDateTime(t) >= when
        ]
        report[station] = {
            "step_time": str(when), "step_log10": step, "step_factor": 10**step,
            "windows_before": int(best), "windows_after": int(len(r) - best),
            "significant": abs(step) > STEP_LIMIT,
            "rayleigh_gain_before": float(np.median(before)) if before else None,
            "rayleigh_gain_after": float(np.median(after)) if after else None,
            "events_before": len(before), "events_after": len(after),
            "dP_db": float(np.median(arr[best:, 2]) - np.median(arr[:best, 2])),
            "dZ_db": float(np.median(arr[best:, 3]) - np.median(arr[:best, 3])),
        }  # fmt: skip
        print(station, json.dumps(report[station]), flush=True)
    (ROOT / "results" / "tables" / "stability_steps.json").write_text(
        json.dumps(report, indent=1)
    )


if __name__ == "__main__":
    main()
