"""Stage C: compliance function of one station for each DPG gain (ComPy 2.2).

Method (Aminian et al. 2025, section 4.2, with the 2.2 stacked estimator):
  * Non-overlapping WINDOW_S windows of cleaned BHZ [m] and BDH [Pa] that do
    not touch earthquake or local-event spans.
  * Gates in the compliance band: median amplitude coherence >= threshold,
    mean pressure PSD 10–50 dB re 1 Pa²/Hz, mean acceleration PSD -170 to
    -150 dB re 1 (m/s²)²/Hz. Threshold 0.9 if it keeps >= MIN_WINDOWS
    windows, otherwise 0.8 (as in the paper).
  * Band: widest run where the mean coherence of accepted windows reaches
    the threshold, searched in (BAND_SEARCH_LOW, f_c) with f_c=sqrt(g/2πH);
    the gates are re-evaluated once in that band.
  * Compliance from stacked spectra (k |S_pz|/S_pp, gravity corrected),
    Bendat & Piersol uncertainty, phase, per-window mean and scatter.
  * Monthly stacks with the same band and threshold test temporal stability.

Gains: 'bruitfm' from config.BRUITFM_GAINS; 'compy' from
work/<STA>/calibration/compy_gain.json (Stage D). Missing gains are skipped.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from obspy import UTCDateTime, read, read_inventory

from config import BRUITFM_GAINS, COMPY, INVENTORY, P, RESULTS, WORK, Z, seed_id

sys.path.insert(0, str(COMPY))

import compy  # noqa: E402
from compy_numerics import GRAVITY  # noqa: E402
from compy_processing import (  # noqa: E402
    SelectionCriteria,
    coherent_band,
    stacked_compliance,
)

WINDOW_S = 2100.0  # 35 min (paper section 4.2)
NSEG = 1024  # 819.2 s Welch segments at 1.25 Hz, df = 1.22 mHz
MIN_WINDOWS = 500
THRESHOLDS = (0.9, 0.8)
BAND_SEARCH_LOW = 0.004
INITIAL_BAND_LOW = 0.007
INITIAL_BAND_HIGH_FRACTION = 0.8
OUTPUT_BAND = (0.002, 0.04)
PAPER_TABLE1 = {}  # published gain and band per station, if any


def gains_for(station):
    gains = {}
    if BRUITFM_GAINS is not None and Path(BRUITFM_GAINS).exists():
        with open(BRUITFM_GAINS) as fh:
            for row in csv.DictReader(fh):
                # Only gains BRUIT-FM actually applied (XE: microseism method);
                # "uncorrected" rows carry g = 1 and are not calibrations.
                if (
                    row["Seed_id"] == seed_id(station)
                    and row.get("gain_status") == "applied"
                ):
                    gains["bruitfm"] = float(row["gain_applied"])
    path = WORK / station / "calibration" / "compy_gain.json"
    if path.exists():
        info = json.loads(path.read_text())
        if info.get("gain"):
            gains["compy"] = float(info["gain"])
    return gains


def load_windows(station):
    clean = WORK / station / "clean"
    spans = json.loads((clean / "spans.json").read_text())
    bad = [
        (UTCDateTime(s), UTCDateTime(e))
        for s, e in spans["earthquakes"] + spans["local"]
    ]
    bad.sort()
    windows, dropped = [], 0
    for path in sorted(clean.glob(f"{station}_seg*.mseed")):
        st = read(str(path))
        st = st.select(channel=Z) + st.select(channel=P)
        for tr in st:
            tr.data = tr.data.astype(np.float64)
        for tr in st:
            tr.data = np.ma.filled(np.ma.masked_invalid(tr.data), np.nan)
        for w in _windows(st):
            t0, t1 = w[0].stats.starttime, w[0].stats.endtime
            if not all(np.all(np.isfinite(tr.data)) for tr in w):
                dropped += 1
                continue
            if any(s < t1 and e > t0 for s, e in bad):
                dropped += 1
                continue
            windows.append(w)
    return windows, dropped


def _windows(st):
    """Non-overlapping windows; NaN samples are allowed here and rejected later."""
    fs = st[0].stats.sampling_rate
    width = int(round(WINDOW_S * fs))
    start = max(tr.stats.starttime for tr in st)
    end = min(tr.stats.endtime for tr in st)
    n = int(round((end - start) * fs)) + 1
    return [
        st.slice(
            start + i / fs, start + (i + width - 1) / fs, nearest_sample=False
        ).copy()
        for i in range(0, n - width + 1, width)
    ]


def run_gain(windows, depth, gain, fc):
    """Threshold choice, coherence band and final stacked estimate."""
    initial = (INITIAL_BAND_LOW, INITIAL_BAND_HIGH_FRACTION * fc)
    attempts = {}
    for threshold in THRESHOLDS:
        criteria = SelectionCriteria(
            threshold, (10.0, 50.0), (-170.0, -150.0), f"gji2025_{threshold}"
        )
        try:
            first = stacked_compliance(windows, depth, gain, initial, criteria, NSEG)
        except ValueError as exc:
            attempts[threshold] = {"error": str(exc)}
            continue
        band = coherent_band(first.all_freqs, first.mean_coherence_all, threshold,
                             (BAND_SEARCH_LOW, fc)) or initial  # fmt: skip
        if band[1] - band[0] < 0.002:
            band = initial
        try:
            final = stacked_compliance(windows, depth, gain, band, criteria, NSEG,
                                       output_band=OUTPUT_BAND)  # fmt: skip
        except ValueError as exc:
            attempts[threshold] = {"error": str(exc), "band": band}
            continue
        attempts[threshold] = {"n_windows": final.n_windows, "band": band}
        if final.n_windows >= MIN_WINDOWS or threshold == THRESHOLDS[-1]:
            return final, threshold, attempts
    return None, None, attempts


def monthly(windows, result, depth, gain):
    months = {}
    accepted = {float(t) for t in result.window_starts}
    for w in windows:
        t = w[0].stats.starttime
        if float(t) in accepted:
            months.setdefault(f"{t.year}-{t.month:02d}", []).append(w)
    rows = {}
    for key, ws in sorted(months.items()):
        if len(ws) < 20:
            continue
        r = stacked_compliance(ws, depth, gain, result.band, result.criteria, NSEG,
                               output_band=result.band)  # fmt: skip
        rows[key] = {"n": r.n_windows, "compliance": r.compliance.tolist(),
                     "uncertainty": r.uncertainty.tolist()}  # fmt: skip
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("station")
    parser.add_argument(
        "--end", help="use only windows before this UTC time (stable epoch)"
    )
    args = parser.parse_args()
    station = args.station
    inv = read_inventory(str(INVENTORY)).select(station=station)
    depth = -inv[0][0].elevation
    fc = float(np.sqrt(GRAVITY / (2 * np.pi * depth)))
    windows, dropped = load_windows(station)
    gains = gains_for(station)
    out = RESULTS / station
    if args.end:
        # Stable epoch: windows before the sensitivity-ratio change (Stage H),
        # and the ComPy gain as the median of per-event gains in that epoch.
        end = UTCDateTime(args.end)
        windows = [w for w in windows if w[0].stats.endtime <= end]
        cal = json.loads(
            (WORK / station / "calibration" / "compy_gain.json").read_text()
        )
        before = [
            g
            for t, g in zip(cal["events_used"], cal["event_gains"])
            if UTCDateTime(t) < end
        ]
        if before:
            gains["compy"] = float(np.median(before))
        out = out / "epoch"
    out.mkdir(parents=True, exist_ok=True)
    summary = {"station": station, "depth_m": depth, "f_c_hz": fc,
               "windows_total": len(windows), "windows_dropped_events": dropped,
               "window_s": WINDOW_S, "nseg": NSEG, "gains": gains,
               "paper_table1": PAPER_TABLE1.get(station), "results": {},
               "epoch_end": args.end}  # fmt: skip
    for label, gain in gains.items():
        result, threshold, attempts = run_gain(windows, depth, gain, fc)
        entry = {"gain": gain, "attempts": {str(k): v for k, v in attempts.items()}}
        if result is None:
            summary["results"][label] = entry
            continue
        in_band = (result.freqs >= result.band[0]) & (result.freqs <= result.band[1])
        entry.update({
            "threshold": threshold, "band_hz": result.band,
            "n_windows": result.n_windows, "n_averages": result.n_averages,
            "rejected": result.rejected,
            "median_coherence_in_band": float(np.median(result.coherence[in_band])),
            "median_phase_deg_in_band": float(np.median(np.abs(result.phase_deg[in_band]))),
            "median_relative_uncertainty": float(np.median(
                result.uncertainty[in_band] / result.compliance[in_band])),
            "monthly": monthly(windows, result, depth, gain),
        })  # fmt: skip
        summary["results"][label] = entry
        np.savez(
            out / f"{station}_compliance_{label}.npz",
            freqs=result.freqs, compliance=result.compliance,
            uncertainty=result.uncertainty, coherence=result.coherence,
            phase_deg=result.phase_deg, window_mean=result.window_mean,
            window_std=result.window_std, band=np.array(result.band),
            all_freqs=result.all_freqs, mean_coherence_all=result.mean_coherence_all,
            S_zz=result.stack["S_zz"], S_pp=result.stack["S_pp"],
            S_pz=result.stack["S_pz"], gain=gain, depth=depth,
            window_starts=np.array([float(t) for t in result.window_starts]),
        )  # fmt: skip
        with open(out / f"{station}_compliance_{label}.csv", "w") as fh:
            fh.write(
                "freq_hz,compliance_1_per_pa,uncertainty_1_per_pa,coherence,phase_deg,in_band\n"
            )
            for i, f in enumerate(result.freqs):
                fh.write(f"{f:.6f},{result.compliance[i]:.6e},{result.uncertainty[i]:.6e},"
                         f"{result.coherence[i]:.4f},{result.phase_deg[i]:.2f},{int(in_band[i])}\n")  # fmt: skip
        print(f"{station} {label}: g={gain:.4f} thr={threshold} band={result.band} "
              f"n={result.n_windows} rej={result.rejected}", flush=True)  # fmt: skip
    (out / f"{station}_summary.json").write_text(
        json.dumps(summary, indent=1, default=str)
    )


if __name__ == "__main__":
    main()
