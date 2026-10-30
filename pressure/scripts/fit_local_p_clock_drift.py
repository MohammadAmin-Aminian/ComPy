"""Fit a linear pressure-logger clock drift from local-P timing lags.

Input: events.csv of a lag-search run. Only events whose lag passed the onset
refinement are used; the fit uses timing alone, never gains. A Theil-Sen line
is fitted, points more than --outlier-s from it are excluded, and the inliers
are refitted by least squares. The model is valid only between the first and
last fitted events (no extrapolation). Weak fits are refused.
Output: a new JSON file (never overwritten) for the `timing_drift` option.
"""

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from obspy import UTCDateTime
from scipy.stats import theilslopes

from compy_pressure.local_p_io import fresh_file

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10


def fit(times, lags, reference, outlier_s):
    days = np.array([(t - reference) / 86400 for t in times])
    lags = np.asarray(lags, dtype=float)
    slope = theilslopes(lags, days)[0]
    # Robust intercept; scipy's default median(y) - slope*median(x) can be
    # moved by a single outlier sitting at the median abscissa.
    intercept = np.median(lags - slope * days)
    inlier = np.abs(lags - (intercept + slope * days)) <= outlier_s
    if inlier.sum() < 2 or np.ptp(days[inlier]) == 0:
        raise ValueError("Need at least two distinct inlier times")
    slope, intercept = np.polyfit(days[inlier], lags[inlier], 1)
    residual = lags[inlier] - (intercept + slope * days[inlier])
    return slope, intercept, days, inlier, residual


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--outlier-s", type=float, default=1.0)
    parser.add_argument("--min-events", type=int, default=8)
    parser.add_argument("--min-span-days", type=float, default=60.0)
    parser.add_argument("--max-rms-s", type=float, default=0.5)
    args = parser.parse_args()
    args.output = fresh_file(args.output)
    with args.events.open(newline="") as stream:
        rows = [
            r
            for r in csv.DictReader(stream)
            if r.get("timing_lag_s") not in (None, "")
            and r.get("reason") != "lag_refinement_failed"
        ]
    times = [UTCDateTime(r["origin_utc"]) for r in rows]
    lags = [float(r["timing_lag_s"]) for r in rows]
    if len(rows) < args.min_events:
        raise SystemExit(f"Only {len(rows)} refined lags; need {args.min_events}")
    reference = UTCDateTime(min(times).date)
    slope, intercept, days, inlier, residual = fit(
        times, lags, reference, args.outlier_s
    )
    span = float(days[inlier].max() - days[inlier].min())
    rms = float(np.sqrt(np.mean(residual**2)))
    if inlier.sum() < args.min_events or span < args.min_span_days:
        raise SystemExit(f"Weak support: {inlier.sum()} inliers over {span:.0f} days")
    if rms > args.max_rms_s:
        raise SystemExit(f"Residual RMS {rms:.2f} s exceeds {args.max_rms_s} s")
    fitted = [t for t, keep in zip(times, inlier, strict=True) if keep]
    model = {
        "reference_utc": str(reference),
        "intercept_s": float(intercept),
        "slope_s_per_day": float(slope),
        "slope_ppm": float(slope / 86400 * 1e6),
        "valid_from_utc": str(min(fitted)),
        "valid_to_utc": str(max(fitted)),
        "inliers": int(inlier.sum()),
        "excluded": [
            {"origin_utc": str(t), "lag_s": lag}
            for t, lag, keep in zip(times, lags, inlier, strict=True)
            if not keep
        ],
        "residual_rms_s": rms,
        "residual_max_s": float(np.max(np.abs(residual))),
        "source": str(args.events.resolve()),
        "source_sha256": hashlib.sha256(args.events.read_bytes()).hexdigest(),
        "created_utc": datetime.now(UTC).isoformat(),
        "method": "Theil-Sen slope, median intercept, outlier cut, LS refit",
        "outlier_s": args.outlier_s,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        stream.write(json.dumps(model, indent=2) + "\n")
    print(json.dumps(model, indent=2))


if __name__ == "__main__":
    main()
