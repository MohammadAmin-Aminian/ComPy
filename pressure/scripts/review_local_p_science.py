"""Bounded scientific review of local-P noise windows, timing and the gate.

Read-only on saved campaign products and on the waveforms of the nine accepted
XO.LT02 HDH events. No threshold, rule, gain, master or saved result changes;
outputs go to a fresh results/LOCAL_P/science_review/ folder.

A. Documented noise-window rule (method note §4.3) versus the implementation.
B. Empirical contamination test: near-onset noise window versus an earlier,
   equal-length background window, and gain sensitivity to that choice.
C. Timing: residual delay and constant phase from saved cross-spectra; held-out
   validation of the LT02 clock-drift line.
D. Synthetic leakage test of the pre-multiple gate with real WD66 responses and
   adversarial post-gate arrivals.
E. Field consistency of pre-multiple diagnostics (WD70, WD66) against geometry.
"""

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from obspy import UTCDateTime
from scipy.stats import spearmanr

from compy_pressure import local_p_check as lc
from compy_pressure import local_p_physics as lp
from compy_pressure.local_p_reference import simulate_plane_wave
from scripts.fit_local_p_clock_drift import fit as fit_drift

UTC = timezone.utc  # datetime.UTC needs Python 3.11; ComPy supports 3.10

CAMPAIGNS = Path.home() / "Projects/bruit-fm-localp/results/LOCAL_P/campaigns"
INVENTORY = Path.home() / "Desktop/Data_Bruit_FM/INVENTORIES/XO.xml"
ARCHIVE = Path.home() / "Desktop/Data_Bruit_FM/Data/XO"
CATALOG = (
    Path.home()
    / "Projects/bruit-fm-localp/results/LOCAL_P/catalogs"
    / "usgs_local_XO_20180520_20190101_m2p5.csv"
)
LT02 = CAMPAIGNS / "apg_timing_v2/XO.LT02/XO.LT02..HDH__XO.LT02..HHZ"
LT02_V1 = CAMPAIGNS / "all_v1/XO.LT02/XO.LT02..HDH__XO.LT02..HHZ"
PAIRS = {
    "XO.LT02..HDH": LT02,
    "XO.WD66..HDH": CAMPAIGNS / "all_v1/XO.WD66/XO.WD66..HDH__XO.WD66..HHZ",
    "XO.WD70..HDH": CAMPAIGNS / "all_v1/XO.WD70/XO.WD70..HDH__XO.WD70..HHZ",
    "YV.RR36.00.BDH": CAMPAIGNS / "all_v1/YV.RR36/YV.RR36.00.BDH__YV.RR36.00.BHZ",
    "YV.RR40.00.BDH": CAMPAIGNS / "all_v1/YV.RR40/YV.RR40.00.BDH__YV.RR40.00.BHZ",
}
BANDS = [(0.5, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, 8.0), (8.0, 16.0)]


def events(folder):
    with (folder / "events.csv").open(newline="") as stream:
        return list(csv.DictReader(stream))


def number(value):
    return float(value) if value not in ("", None) else None


# ---------------------------------------------------------------- A. rule audit


def documented_boundary(row):
    """Predicted P minus the pre-registered onset-search half-width (§4.1-4.3)."""
    t_pred = UTCDateTime(row["predicted_p_utc"])
    travel = t_pred - UTCDateTime(row["origin_utc"])
    geometry = SimpleNamespace(p_travel_s=travel)
    event = {
        "horizontal_error_km": number(row["horizontal_error_km"]),
        "depth_error_km": number(row["depth_error_km"]),
    }
    return t_pred - lc.search_half_width(geometry, event, lc.DEFAULTS)


def rule_audit(rows):
    out = []
    for row in rows:
        if row["status"] != "accepted":
            continue
        noise_end = UTCDateTime(row["noise_start_utc"]) + float(row["window_seconds"])
        boundary = documented_boundary(row)
        onset = UTCDateTime(row["onset_vertical_utc"])
        out.append(
            {
                "event_id": row["event_id"],
                "excess_s": round(noise_end - boundary, 3),
                "violates_documented_rule": noise_end > boundary,
                "onset_minus_predicted_s": round(
                    onset - UTCDateTime(row["predicted_p_utc"]), 3
                ),
                "noise_end_before_onset_s": round(onset - noise_end, 3),
            }
        )
    return out


def campaign_rule_counts():
    counts = {}
    for summary in sorted(CAMPAIGNS.glob("all_v1/*/*/events.csv")):
        rows = [r for r in events(summary.parent) if r["status"] == "accepted"]
        if rows:
            audit = rule_audit(rows)
            key = summary.parent.name.split("__")[0]
            counts[key] = [
                sum(a["violates_documented_rule"] for a in audit),
                len(audit),
            ]
    return counts


# ------------------------------------------------- B. contamination/sensitivity


def lt02_context():
    p_epochs = lc.channel_epochs(INVENTORY, "XO.LT02..HDH")
    v_epochs = lc.channel_epochs(INVENTORY, "XO.LT02..HHZ")
    return SimpleNamespace(
        p_epochs=p_epochs,
        v_epochs=v_epochs,
        p_index=lc.DayFileIndex([ARCHIVE / "LT02"], "XO.LT02..HDH"),
        v_index=lc.DayFileIndex([ARCHIVE / "LT02"], "XO.LT02..HHZ"),
        rate=100.0,
        quantum=lc.duration_quantum([100.0, 125.0]),
        catalog=lc.load_local_catalog(CATALOG),
        model=lc.TauPyModel("iasp91"),
    )


def band_median(f, values, band, mask):
    sel = mask & (f >= band[0]) & (f < band[1]) & np.isfinite(values)
    return float(np.median(values[sel])) if sel.sum() >= 3 else None


def contamination(row, ctx, gap_s=30.0, depth=152.0):
    """Compare the near-onset noise window with an earlier background window."""
    o = lc.DEFAULTS
    rate = ctx.rate
    n = int(round(float(row["window_seconds"]) * rate))
    window_start = UTCDateTime(row["window_start_utc"])
    noise_start = UTCDateTime(row["noise_start_utc"])
    early_start = noise_start - gap_s - n / rate
    # Control: a third window as far before `early` as `early` is before `near`,
    # giving the natural background non-stationarity for the ratio test.
    earlier_start = early_start - gap_s - n / rate
    seg_start = earlier_start - 10.0
    seconds = lc.quantize(
        window_start + n / rate + 10.0 - seg_start, ctx.quantum, down=False
    )
    epoch_p = lc.epoch_for(ctx.p_epochs, seg_start, seg_start + seconds)
    epoch_v = lc.epoch_for(ctx.v_epochs, seg_start, seg_start + seconds)
    event = {"id": row["event_id"], "time": UTCDateTime(row["origin_utc"])}
    other = lc.overlapping_event(
        event,
        ctx.catalog,
        ctx.model,
        epoch_v.station.latitude,
        epoch_v.station.longitude,
        earlier_start - o["overlap_before_s"],
        window_start + n / rate,
    )
    p_seg, v_seg, _, _ = lc.load_pair(
        ctx.p_index,
        ctx.v_index,
        epoch_p,
        epoch_v,
        seg_start,
        seconds,
        rate,
        float(row["timing_lag_s"]),
        o,
    )
    index = lambda t: int(round((t - seg_start) * rate))  # noqa: E731
    signal = slice(index(window_start), index(window_start) + n)
    near = slice(index(noise_start), index(noise_start) + n)
    early = slice(index(early_start), index(early_start) + n)
    earlier = slice(index(earlier_start), index(earlier_start) + n)
    responses = lc.responses_for(epoch_p, epoch_v)
    eta = float(row["eta_taup"])
    result = {"event_id": row["event_id"], "other_event_in_span": other}
    gains = {}
    for label, noise in (("near", near), ("early", early), ("earlier", earlier)):
        bands, arrays, _ = lc.analyse_event(
            p_seg, v_seg, signal, noise, rate, responses, "steady_state", depth, eta, o
        )
        gains[label] = {
            f"{b['band_hz'][0]:g}-{b['band_hz'][1]:g}": b.get("gain_noise_corrected")
            for b in bands
        }
        if label == "near":
            f = arrays["frequency_hz"]
            near_pp, near_vv = arrays["n_obs_obs"], arrays["n_vv"]
        elif label == "early":
            early_pp, early_vv = arrays["n_obs_obs"], arrays["n_vv"]
        else:
            earlier_pp, earlier_vv = arrays["n_obs_obs"], arrays["n_vv"]
    in_band = (f >= 0.5) & (f < 2.0)
    result["noise_ratio_pressure_0p5_2hz"] = float(
        np.median(near_pp[in_band] / early_pp[in_band])
    )
    result["noise_ratio_vertical_0p5_2hz"] = float(
        np.median(near_vv[in_band] / early_vv[in_band])
    )
    result["control_ratio_pressure_0p5_2hz"] = float(
        np.median(early_pp[in_band] / earlier_pp[in_band])
    )
    result["control_ratio_vertical_0p5_2hz"] = float(
        np.median(early_vv[in_band] / earlier_vv[in_band])
    )
    result["gains_near_noise"] = gains["near"]
    result["gains_early_noise"] = gains["early"]
    result["saved_gains"] = {
        "0.5-1": number(row["0.5_1hz_gain_noise_corrected"]),
        "1-2": number(row["1_2hz_gain_noise_corrected"]),
    }
    return result


# ------------------------------------------------------------------- C. timing


def residual_timing(folder):
    """Per-event delay tau and constant phase phi0 fitted to obs-minus-pred phase."""
    rows = {r["event_id"]: r for r in events(folder) if r["status"] == "accepted"}
    data = np.load(folder / "spectra.npz")
    taus = np.arange(-0.5, 0.5001, 0.001)
    out = []
    for event_id in rows:
        f = data[f"{event_id}__frequency_hz"]
        valid = data[f"{event_id}__valid"] & (f >= 0.5) & (f <= 16)
        if valid.sum() < 5:
            continue
        phase = np.deg2rad(data[f"{event_id}__phase_obs_minus_pred_deg"][valid])
        fv = f[valid]
        power = np.abs(
            np.exp(
                1j * (phase[None, :] - 2 * np.pi * fv[None, :] * taus[:, None])
            ).mean(axis=1)
        )
        k = int(np.argmax(power))
        residual = phase - 2 * np.pi * fv * taus[k]
        phi0 = np.angle(np.mean(np.exp(1j * residual)))
        out.append(
            {
                "event_id": event_id,
                "bins": int(valid.sum()),
                "f_range_hz": [round(float(fv.min()), 2), round(float(fv.max()), 2)],
                "tau_s": round(float(taus[k]), 3),
                "phi0_deg": round(float(np.rad2deg(phi0)), 1),
                "phase_resultant": round(float(power[k]), 3),
                "raw_mean_phase_deg": round(
                    float(np.rad2deg(np.angle(np.mean(np.exp(1j * phase))))), 1
                ),
            }
        )
    return out


def drift_holdout():
    rows = [
        r
        for r in events(LT02_V1)
        if r.get("timing_lag_s") not in ("", None)
        and r["reason"] != "lag_refinement_failed"
    ]
    rows.sort(key=lambda r: r["origin_utc"])
    times = [UTCDateTime(r["origin_utc"]) for r in rows]
    lags = np.array([float(r["timing_lag_s"]) for r in rows])
    reference = UTCDateTime(times[0].date)
    days = np.array([(t - reference) / 86400 for t in times])
    half = len(rows) // 2
    slope, intercept, _, inlier, _ = fit_drift(
        times[:half], lags[:half], reference, 1.0
    )
    forward = lags[half:] - (intercept + slope * days[half:])
    loo = []
    for i in range(len(rows)):
        keep = [j for j in range(len(rows)) if j != i]
        s, c, _, _, _ = fit_drift([times[j] for j in keep], lags[keep], reference, 1.0)
        loo.append(lags[i] - (c + s * days[i]))
    loo = np.array(loo)
    core = np.abs(loo) < 1.0
    refinement = [
        float(r["timing_lag_s"]) - float(r["timing_lag_coarse_s"])
        for r in events(LT02)
        if r["status"] == "accepted"
    ]
    return {
        "events": len(rows),
        "fit_first_half": {
            "n": half,
            "slope_s_per_day": slope,
            "intercept_s": intercept,
        },
        "forward_residuals_s": [round(x, 3) for x in forward],
        "forward_rms_s": float(np.sqrt(np.mean(forward**2))),
        "leave_one_out_rms_s_without_gross": float(np.sqrt(np.mean(loo[core] ** 2))),
        "leave_one_out_gross_outliers": int((~core).sum()),
        "v2_refinement_shifts_s": [round(x, 3) for x in refinement],
    }


# -------------------------------------------------------- D. synthetic leakage


def synthetic_leakage(gain=0.8):
    """Gate leakage with real WD66 responses and adversarial post-gate arrivals."""
    depth, rate = 4940.0, 100.0
    t, p, v, onset = simulate_plane_wave(
        depth=depth,
        dz=4.0,
        duration=26.0,
        peak_frequency=6.0,
        source_below_seafloor=800.0,
        solid_speed=3000.0,
        solid_density=2200.0,
    )
    grid = np.arange(0, t[-1], 1 / rate)
    p, v = np.interp(grid, t, p), np.interp(grid, t, v)
    pad = int(40 * rate)
    p = np.concatenate([np.zeros(pad), p])
    v = np.concatenate([np.zeros(pad), v])
    i_on = int(round(onset * rate)) + pad
    o = lc.DEFAULTS
    p_epoch = lc.channel_epochs(INVENTORY, "XO.WD66..HDH")[0]
    v_epoch = lc.channel_epochs(INVENTORY, "XO.WD66..HHZ")[0]
    responses = lc.responses_for(p_epoch, v_epoch)
    tau_min = 2 * depth * 0.95 * lp.eta_min_for_p() / o["sound_speed_range"][1]
    gate = o["pre_multiple_safety"] * tau_min + o["lead_s"]
    true_multiple = lp.round_trip_time(depth)
    peak = np.max(np.abs(v[i_on : i_on + int(gate * rate)]))
    rng = np.random.default_rng(11)

    def ricker(center, frequency=3.0):
        x = (np.pi * frequency * (np.arange(len(v)) / rate - center)) ** 2
        return (1 - 2 * x) * np.exp(-x)

    def to_counts(x, which):
        nfft = 1 << int(np.ceil(np.log2(2 * len(x))))
        f = np.fft.rfftfreq(nfft, 1 / rate)
        rp, rv, _, _ = responses(f)
        response = rp if which == "p" else rv
        return np.fft.irfft(np.fft.rfft(x, nfft) * response, nfft)[: len(x)]

    cases = []
    start = i_on / rate - o["lead_s"]
    for injection in (
        "none",
        "S_like_eta0.9",
        "reversed_sign",
        "pressure_free_vertical",
    ):
        vv, pp = v.copy(), p.copy()
        t_s = start + gate + 0.3
        if injection != "none":
            burst = 5 * peak * ricker(t_s)
            vv = vv + burst
            if injection == "S_like_eta0.9":
                pp = pp + 1028 * 1500 / 0.9 * burst
            elif injection == "reversed_sign":
                pp = pp - 1028 * 1500 * burst
        for label, length in (
            ("campaign_gate", gate),
            ("gate_minus_1s", gate - 1.0),
            ("past_multiple", true_multiple + 2.0),
        ):
            n = int(round(length * rate))
            s0 = int(round(start * rate))
            noise = slice(
                s0 - int(o["noise_gap_s"] * rate) - n, s0 - int(o["noise_gap_s"] * rate)
            )
            p_counts = to_counts(pp / gain, "p") + 1e-6 * rng.standard_normal(len(pp))
            v_counts = to_counts(vv, "v") + 1e-6 * rng.standard_normal(len(vv))
            bands, _, _ = lc.analyse_event(
                p_counts,
                v_counts,
                slice(s0, s0 + n),
                noise,
                rate,
                responses,
                "pre_multiple",
                depth,
                1.0,
                o,
            )
            cases.append(
                {
                    "injection": injection,
                    "gate": label,
                    "gate_s": round(length, 3),
                    "bias_percent_by_band": {
                        f"{b['band_hz'][0]:g}-{b['band_hz'][1]:g}": round(
                            100 * (b["gain_noise_corrected"] / gain - 1), 2
                        )
                        for b in bands
                        if "gain_noise_corrected" in b
                    },
                }
            )
    return {
        "true_gain": gain,
        "gate_s": gate,
        "true_first_multiple_s": true_multiple,
        "injection_time_after_window_start_s": round(gate + 0.3, 3),
        "cases": cases,
    }


# ----------------------------------------------------------- E. field geometry


def event_quality(folder, event_id, band):
    """Median pressure SNR and coherence over an event's valid bins in a band."""
    data = np.load(folder / "spectra.npz")
    f = data[f"{event_id}__frequency_hz"]
    sel = data[f"{event_id}__valid"] & (f >= band[0]) & (f < band[1])
    snr = data[f"{event_id}__s_obs_obs"][sel] / data[f"{event_id}__n_obs_obs"][sel]
    return float(np.median(snr)), float(np.median(data[f"{event_id}__coherence"][sel]))


def field_consistency(folder, bands=("1_2hz", "2_4hz")):
    rows = [r for r in events(folder) if r["status"] == "accepted"]
    out = {}
    for band in bands:
        sel = [r for r in rows if r.get(f"{band}_gain_noise_corrected")]
        g = np.array([float(r[f"{band}_gain_noise_corrected"]) for r in sel])
        if len(g) < 6:
            continue
        stats = {
            "n": len(g),
            "median": float(np.median(g)),
            "relative_mad": float(np.median(np.abs(g - np.median(g))) / np.median(g)),
        }
        for key in (
            "distance_km",
            "depth_km",
            "magnitude",
            "eta_taup",
            "window_seconds",
        ):
            x = np.array([float(r[key]) for r in sel])
            rho, pval = spearmanr(x, g)
            stats[f"spearman_{key}"] = [round(float(rho), 3), round(float(pval), 3)]
        limits = tuple(float(x) for x in band.replace("hz", "").split("_"))
        quality = np.array([event_quality(folder, r["event_id"], limits) for r in sel])
        for i, key in enumerate(("pressure_snr", "coherence")):
            rho, pval = spearmanr(quality[:, i], g)
            stats[f"spearman_{key}"] = [round(float(rho), 3), round(float(pval), 3)]
        stats["tests_per_band"] = 7
        out[band] = stats
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--skip-waveforms", action="store_true")
    args = parser.parse_args()
    out = args.output_dir.expanduser().resolve()
    if "/results/LOCAL_P/" not in f"{out}/" or out.exists():
        raise SystemExit("Output must be a new folder under results/LOCAL_P/")
    report = {"created_utc": datetime.now(UTC).isoformat(), "inputs": {}}
    for folder in [*PAIRS.values(), LT02_V1]:
        for name in ("events.csv", "spectra.npz", "summary.json"):
            path = folder / name
            report["inputs"][str(path)] = lc.sha256(path)
    lt02_rows = events(LT02)
    report["A_rule_audit_lt02"] = rule_audit(lt02_rows)
    report["A_rule_counts_all_v1"] = campaign_rule_counts()
    if not args.skip_waveforms:
        ctx = lt02_context()
        report["B_contamination_lt02"] = [
            contamination(r, ctx) for r in lt02_rows if r["status"] == "accepted"
        ]
    report["C_residual_timing"] = {k: residual_timing(v) for k, v in PAIRS.items()}
    report["C_drift_holdout"] = drift_holdout()
    report["D_synthetic_leakage"] = synthetic_leakage()
    report["E_field_consistency"] = {
        k: field_consistency(PAIRS[k]) for k in ("XO.WD70..HDH", "XO.WD66..HDH")
    }
    after = {p: lc.sha256(p) for p in report["inputs"]}
    report["inputs_unchanged"] = after == report["inputs"]
    out.mkdir(parents=True)
    (out / "science_review.json").write_text(
        json.dumps(report, indent=2, default=str) + "\n"
    )
    print(
        json.dumps(
            {k: v for k, v in report.items() if k not in {"inputs"}},
            indent=1,
            default=str,
        )[:20000]
    )


if __name__ == "__main__":
    main()
