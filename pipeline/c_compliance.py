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
  * Covariance of the stack from its own windows (compy.inversion.errors:
    bootstrap over windows, slow level and tilt from 30-day blocks), saved with
    the result as cov_rel for the inversion (ComPy 3).

Gains: 'bruitfm' from config.BRUITFM_GAINS; 'compy' from
work/<STA>/calibration/compy_gain.json (Stage D). Missing gains are skipped.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import config
import numpy as np
from config import BRUITFM_GAINS, COMPY, INVENTORY, RESULTS, WORK, P, Z, seed_id
from obspy import UTCDateTime, read, read_inventory

sys.path.insert(0, str(COMPY))

from compy.compliance import (  # noqa: E402
    SelectionCriteria,
    coherent_band,
    spectra,
    stacked_compliance,
)
from compy.inversion.errors import stack_covariance, usable_bins  # noqa: E402
from compy.numerics import GRAVITY  # noqa: E402
from compy.response_check import paz_stages  # noqa: E402

WINDOW_S = 2100.0  # 35 min (paper section 4.2)
NSEG = 1024  # 819.2 s Welch segments at 1.25 Hz, df = 1.22 mHz
MIN_WINDOWS = 500
THRESHOLDS = (0.9, 0.8)
BAND_SEARCH_LOW = 0.004
INITIAL_BAND_LOW = 0.007
INITIAL_BAND_HIGH_FRACTION = 0.8
#: Frequency range saved. The upper limit rises to the band's upper edge where the
#: band goes beyond it (shelf stations): before ComPy 3 the saved curve stopped at
#: 40 mHz while the band and the monthly curves went on, so stage E lost its
#: month-to-month covariance there.
OUTPUT_BAND = (0.002, 0.04)
PAPER_TABLE1 = {}  # published gain and band per station, if any
#: Coherence gate of ``stacked_compliance``: "inclusive" (Aminian et al. 2025,
#: the default) or "exclusive" (per-bin gate without the selection bias; see
#: docs/audit_2026-11/REPORT.md). Set COMPLIANCE_GATE in the run configuration.
GATE = getattr(config, "COMPLIANCE_GATE", "inclusive")
#: Save ``<STA>_window_levels.npy`` (per accepted window: time, log10 median
#: |S_pz|/S_pp, pressure dB, acceleration dB in the band), which
#: ``compy.results`` loads and ``compy.anomaly.scan_levels`` needs for steps,
#: tides and regimes. A new file only; no result changes.
SAVE_WINDOW_LEVELS = getattr(config, "SAVE_WINDOW_LEVELS", True)
#: Gain label whose levels are saved, in this order of preference.
LEVEL_LABELS = ("compy", "bruitfm", "nominal")


def gains_for(station):
    """Available pressure gains of a station."""
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
        if not info.get("error") and info.get("gain"):
            gains["compy"] = float(info["gain"])
    if not gains:
        # No calibration: nominal response (g = 1), flagged; stage E samples
        # the gain over a wide range (e_invert.GAIN_BOUNDS["nominal"]).
        gains["nominal"] = 1.0
    return gains


def load_windows(station):
    """Cleaned 35-min windows outside event spans, and the number dropped."""
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
            first = stacked_compliance(
                windows, depth, gain, initial, criteria, NSEG, gate=GATE
            )
        except ValueError as exc:
            attempts[threshold] = {"error": str(exc)}
            continue
        band = coherent_band(first.all_freqs, first.mean_coherence_all, threshold,
                             (BAND_SEARCH_LOW, fc)) or initial  # fmt: skip
        if band[1] - band[0] < 0.002:
            band = initial
        try:
            final = stacked_compliance(windows, depth, gain, band, criteria, NSEG,
                                       output_band=(OUTPUT_BAND[0], max(OUTPUT_BAND[1], band[1])),
                                       gate=GATE)  # fmt: skip
        except ValueError as exc:
            attempts[threshold] = {"error": str(exc), "band": band}
            continue
        attempts[threshold] = {"n_windows": final.n_windows, "band": band}
        if final.n_windows >= MIN_WINDOWS or threshold == THRESHOLDS[-1]:
            return final, threshold, attempts
    return None, None, attempts


def accepted_spectra(windows, result):
    """Spectra of the accepted windows in ``result.band`` (pressure before gain).

    Returns:
        tuple: frequencies, start times, and S_zz, S_pp, S_pz as (windows, frequencies).
    """
    accepted = {float(t) for t in result.window_starts}
    starts, zz, pp, zp = [], [], [], []
    f_band = None
    for w in windows:
        t = float(w[0].stats.starttime)
        if t not in accepted:
            continue
        f, a, b, x, _ = spectra(w, NSEG)
        sel = (f >= result.band[0] - 1e-12) & (f <= result.band[1] + 1e-12)
        f_band = f[sel]
        starts.append(t)
        zz.append(a[sel])
        pp.append(b[sel])
        zp.append(x[sel])
    return f_band, np.array(starts), np.array(zz), np.array(pp), np.array(zp)


def window_levels(windows, result, gain):
    """Per-window levels of the accepted windows, in ``result.band``.

    Returns:
        ndarray: Rows (POSIX start time, log10 median |S_pz|/S_pp,
        pressure PSD dB re 1 Pa^2/Hz, acceleration PSD dB re 1 (m/s^2)^2/Hz),
        the format of ``compy.results.StationProducts.window_levels``.
    """
    return _levels(accepted_spectra(windows, result), gain)


def _levels(spec, gain):
    f, starts, zz, pp, zp = spec
    rows = []
    for t, a, b, x in zip(starts, zz, pp * gain**2, zp * gain):
        rows.append((t, float(np.log10(np.median(np.abs(x) / b))),
                     float(10 * np.log10(np.mean(b))),
                     float(10 * np.log10(np.mean(a * (2 * np.pi * f) ** 4)))))  # fmt: skip
    return np.array(rows, float).reshape(-1, 4)


def stack_errors(spec, depth, gain):
    """Covariance of the stacked compliance from its windows, or None with too few.

    Returns:
        dict: ``cov_freqs``, ``cov_rel`` (relative covariance over the band),
        ``cov_usable`` (bins that ``usable_bins`` keeps), ``slow_level``,
        ``slow_tilt`` and ``cov_block_days``, as saved with the result.
    """
    f, starts, _, pp, zp = spec
    if len(starts) < 20:
        return None
    stat = stack_covariance(f, pp, zp, starts, depth, gain=gain)
    return dict(cov_freqs=f, cov_rel=stat["cov_rel"], cov_usable=usable_bins(stat["sd_rel"]),
                slow_level=stat["slow_level"], slow_tilt=stat["slow_tilt"],
                cov_block_days=stat["block_days"] or 0.0)  # fmt: skip


def response_check(result, stages, calibration_band=(0.02, 0.06)):
    """Low-frequency check of the pressure response from the Z/P phase. Nothing is applied.

    The Z/P phase of an elastic seafloor is 0 or 180 degrees. A differential
    pressure gauge whose leak corner differs from the metadata leaves a residual
    that grows towards low frequency (Bell et al. 2015, Fig. 14; Doran & Laske
    2021, section 2.2; Webb et al. 2026 measure leak time constants of 100 to
    300 s that change with depth and temperature). The factor ``k`` on the
    lowest metadata pole that fits the residual, and the compliance bias it
    implies after a gain fitted in ``calibration_band``, go to the summary so
    that a run shows how large the effect is at each station.

    Args:
        result (ComplianceResult): Stacked compliance of the station.
        stages (list): ``compy.response_check.paz_stages`` of the pressure channel.

    Returns:
        dict or None: None without poles-and-zeros stages or with fewer than
        three bins in the band.
    """
    from compy.anomaly import residual_phase
    from compy.response_check import _pole_hz, compliance_bias, fit_leak_factor

    in_band = (result.freqs >= result.band[0]) & (result.freqs <= result.band[1])
    poles = [_pole_hz(p, rad) for _, ps, rad in stages for p in ps if abs(p) > 0]
    if not poles or in_band.sum() < 3:
        return None
    f = result.freqs[in_band]
    residual, polarity = residual_phase(result.phase_deg[in_band])
    fit = fit_leak_factor(stages, f, residual)
    bias = compliance_bias(
        stages, f[[0, -1]], fit["k"], calibration_band=calibration_band
    )
    return dict(fit, polarity_deg=polarity, median_residual_deg=float(np.median(residual)),
                metadata_corner_hz=float(min(poles)), fitted_corner_hz=float(min(poles) * fit["k"]),
                bias_percent_band_low=float(bias[0]), bias_percent_band_high=float(bias[1]),
                applied=False)  # fmt: skip


def monthly(windows, result, depth, gain):
    """Monthly compliance stacks of the accepted windows."""
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
                               output_band=result.band, gate=GATE)  # fmt: skip
        rows[key] = {"n": r.n_windows, "compliance": r.compliance.tolist(),
                     "uncertainty": r.uncertainty.tolist()}  # fmt: skip
    return rows


def main():
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("station")
    parser.add_argument(
        "--end", help="use only windows before this UTC time (stable epoch)"
    )
    args = parser.parse_args()
    station = args.station
    inv = read_inventory(str(INVENTORY)).select(station=station)
    depth = -inv[0][0].elevation
    pressure = inv.select(channel=P)
    stages = (
        paz_stages(pressure[0][0][0]) if len(pressure) and len(pressure[0][0]) else []
    )
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
               "window_s": WINDOW_S, "nseg": NSEG, "gains": gains, "gate": GATE,
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
            "response_check": response_check(
                result, stages, getattr(config, "CALIBRATION_BAND", (0.02, 0.06))),
        })  # fmt: skip
        summary["results"][label] = entry
        spec = accepted_spectra(windows, result)
        errors = stack_errors(spec, depth, gain) or {}
        if errors:
            entry.update(slow_level=errors["slow_level"], slow_tilt=errors["slow_tilt"],
                         noisy_bins_cut=int((~errors["cov_usable"]).sum()))  # fmt: skip
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
            **errors,
        )  # fmt: skip
        with open(out / f"{station}_compliance_{label}.csv", "w") as fh:
            fh.write(
                "freq_hz,compliance_1_per_pa,uncertainty_1_per_pa,coherence,phase_deg,in_band\n"
            )
            for i, f in enumerate(result.freqs):
                fh.write(f"{f:.6f},{result.compliance[i]:.6e},{result.uncertainty[i]:.6e},"
                         f"{result.coherence[i]:.4f},{result.phase_deg[i]:.2f},{int(in_band[i])}\n")  # fmt: skip
        level_label = next((x for x in LEVEL_LABELS if x in gains), None)
        if SAVE_WINDOW_LEVELS and label == level_label:
            np.save(out / f"{station}_window_levels.npy", _levels(spec, gain))
        print(f"{station} {label}: g={gain:.4f} thr={threshold} band={result.band} "
              f"n={result.n_windows} rej={result.rejected}", flush=True)  # fmt: skip
    (out / f"{station}_summary.json").write_text(
        json.dumps(summary, indent=1, default=str)
    )


if __name__ == "__main__":
    main()
