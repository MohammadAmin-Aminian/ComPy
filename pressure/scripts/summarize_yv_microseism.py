"""Summarize a completed YV microseism run: phase screen and gain comparison.

Reads only the run folder and tracked metadata; recomputes no spectra and
writes no gain. Theory requires a real P/Z transfer function for any layered
elastic medium and plane-wave mixture: phase 0 deg (sign -1, first branch) or
180 deg. Each station's coherent-bin phase is compared with four hypotheses:
nominal response, the recorded positive-real Laplace pressure pole reflected
to the left half-plane (phase -2 atan(p/omega)), and both with reversed
relative polarity. Labels are screening labels, not scores.
"""

import argparse
import csv
import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from compy_pressure.microseism_check import sha256
from compy_pressure.ppsd import atomic_json

ROOT = Path(__file__).resolve().parents[1]
RESPONSE_FLAGS = ROOT / "docs/yv_rayleigh_independent_validation.json"
COHERENT = 0.5  # phase standard error ~ sqrt((1-g2)/(2 n g2)) ~ 3 deg at n=167
PHASE_BAND_HZ = (0.05, 0.30)
FITS_DEG, UNEXPLAINED_DEG = 10.0, 20.0
INFORMATIVE_CONDITIONING = 1.5
# Amplitude is unaffected by pole reflection or polarity; unexplained phase is
# a response or signal problem, so its amplitude is not interpreted.
EXPLAINED = {
    "nominal",
    "reflected_pole",
    "reversed_nominal",
    "reversed_reflected_pole",
}


def wrap(degrees):
    return (np.asarray(degrees) + 180.0) % 360.0 - 180.0


def circular_mean(degrees):
    radians = np.radians(degrees)
    return float(np.degrees(np.arctan2(np.sin(radians).mean(), np.cos(radians).mean())))


def reflected_pole_phase(f, poles):
    """Phase of D/P_nominal if each recorded +p pole is truly at -p."""
    return -sum(np.degrees(2 * np.arctan(p / (2 * np.pi * f))) for p in poles)


def phase_screen(f, phase, coherence, poles, first_resonance, vertical_poles=()):
    """Vertical poles enter with opposite sign because H = D / P."""
    m = (coherence >= COHERENT) & (f >= PHASE_BAND_HZ[0]) & (f <= PHASE_BAND_HZ[1])
    if m.sum() < 10:
        return {"label": "insufficient_coherent_bins", "coherent_bins": int(m.sum())}
    f, phase = f[m], phase[m]
    pole = reflected_pole_phase(f, poles) - reflected_pole_phase(f, vertical_poles)
    models = {
        "nominal": np.zeros_like(f),
        "reflected_pole": pole,
        "reversed_nominal": np.full_like(f, 180.0),
        "reversed_reflected_pole": 180.0 + pole,
    }
    rms = {k: float(np.sqrt(np.mean(wrap(phase - v) ** 2))) for k, v in models.items()}
    if not poles and not vertical_poles:  # pole hypotheses equal nominal ones
        rms.pop("reflected_pole")
        rms.pop("reversed_reflected_pole")
    best = min(rms, key=rms.get)
    below = f < first_resonance
    result = {
        "coherent_bins": int(m.sum()),
        "frequency_hz": [float(f.min()), float(f.max())],
        "circular_mean_phase_deg": circular_mean(phase),
        "bins_below_first_resonance": int(below.sum()),
        "circular_mean_phase_below_resonance_deg": (
            circular_mean(phase[below]) if below.any() else None
        ),
        "rms_residual_deg": rms,
        "best_hypothesis": best,
    }
    if rms[best] > UNEXPLAINED_DEG:
        result["label"] = "unexplained_phase"
    elif rms[best] > FITS_DEG:
        result["label"] = f"approximate_{best}"
    else:
        result["label"] = best
    return result


def load_bins(path):
    with Path(path).open() as stream:
        rows = list(csv.DictReader(stream))
    return {
        key: np.array([float(r[key]) for r in rows])
        for key in ("frequency_hz", "coherence_squared", "phase_deg")
    }


def summarize(run, docs=None, show=False):
    run = Path(run)
    flags = json.loads(RESPONSE_FLAGS.read_text())["response_flags"]
    poles = {
        r["station"]: [p[0] for p in r["pressure_positive_real_poles"] if p[0] > 0]
        for r in flags
    }
    stations = json.loads((run / "summary.json").read_text())["stations"]
    rows, details = [], {}
    for record in stations:
        station = record["station"]
        row = {"station": station, "status": record["status"]}
        if record["status"] != "checked":
            row["reason"] = record.get("paired_reason")
            rows.append(row)
            continue
        check = json.loads(Path(record["check_summary"]).read_text())
        bins = load_bins(check["files"]["bins"])
        screen = phase_screen(
            bins["frequency_hz"],
            bins["phase_deg"],
            bins["coherence_squared"],
            poles.get(station, []),
            record["first_resonance_hz"],
        )
        details[station] = screen
        conditioning = record["best_sub_band_conditioning_ratio"]
        row.update(
            depth_m=record["depth_m"],
            interval_start=record["start_utc"][:10],
            positive_pressure_poles_rad_s=";".join(
                f"{p:g}" for p in poles.get(station, [])
            ),
            phase_label=screen["label"],
            phase_mean_deg=screen.get("circular_mean_phase_deg"),
            phase_rms_best_deg=(
                min(screen["rms_residual_deg"].values())
                if "rms_residual_deg" in screen
                else None
            ),
            phase_below_resonance_deg=screen.get(
                "circular_mean_phase_below_resonance_deg"
            ),
            long_wave_coherence_squared=record["long_wave_coherence_squared"],
            conditioning_ratio=conditioning,
            gain_informative=bool(
                conditioning <= INFORMATIVE_CONDITIONING
                and screen["label"] in EXPLAINED
            ),
            best_sub_band_hz="-".join(f"{v:.4f}" for v in record["best_sub_band_hz"]),
            microseism_g_h1=record["best_sub_band_mode0_water_h1"],
            microseism_g_h2=record["best_sub_band_mode0_water_h2"],
            upper_bound_h2=record["best_sub_band_upper_bound_h2"],
            baseline_g=record.get("baseline_g"),
            trial_g=record.get("trial_g"),
            published_g=record.get("published_g"),
        )
        for source in ("baseline", "trial", "published"):
            if record.get(f"{source}_g") is not None:
                row[f"{source}_vs_h1_pct"] = record.get(f"{source}_best_h1_pct")
                row[f"{source}_vs_h2_pct"] = record.get(f"{source}_best_h2_pct")
        rows.append(row)

    fields = sorted({k for r in rows for k in r}, key=lambda k: (k != "station", k))
    with (run / "phase_and_gain_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    counts = {}
    for r in rows:
        key = r.get("phase_label", r["status"])
        counts[key] = counts.get(key, 0) + 1
    summary = {
        "run": str(run),
        "run_summary_sha256": sha256(run / "summary.json"),
        "response_flags": str(RESPONSE_FLAGS),
        "response_flags_sha256": sha256(RESPONSE_FLAGS),
        "coherent_bin_threshold": COHERENT,
        "phase_band_hz": list(PHASE_BAND_HZ),
        "label_thresholds_deg": {"fits": FITS_DEG, "unexplained": UNEXPLAINED_DEG},
        "informative_conditioning_ratio": INFORMATIVE_CONDITIONING,
        "phase_label_counts": counts,
        "stations": details,
        "status": "diagnostic_only; screening labels are not scores",
    }
    atomic_json(run / "phase_and_gain_summary.json", summary)
    figure = plot(rows, run / "phase_and_gain_summary.png")
    if show:
        plt.show()
    plt.close(figure)
    if docs:
        docs = Path(docs)
        for name, target in (
            ("phase_and_gain_summary.csv", docs / "yv_microseism_check.csv"),
            ("phase_and_gain_summary.json", docs / "yv_microseism_check.json"),
            ("phase_and_gain_summary.png", docs / "figures/yv_microseism_check.png"),
        ):
            shutil.copyfile(run / name, target)
    return rows, summary


def plot(rows, destination):
    done = [r for r in rows if r["status"] == "checked"]
    x = np.arange(len(done))
    figure, axes = plt.subplots(2, 1, figsize=(15, 10), sharex=True)
    colours = {
        "nominal": "tab:green",
        "reflected_pole": "tab:blue",
        "reversed_reflected_pole": "tab:red",
        "reversed_nominal": "tab:red",
    }
    for i, r in enumerate(done):
        colour = colours.get(r["phase_label"].replace("approximate_", ""), "k")
        if r["phase_mean_deg"] is not None:
            axes[0].plot(i, r["phase_mean_deg"], "o", color=colour)
    for value in (0, 180, -180):
        axes[0].axhline(value, color="gray", lw=0.8)
    axes[0].set(
        ylabel="circular mean P/Z phase (deg)",
        ylim=(-200, 200),
        title="Coherent-bin phase (0.05–0.30 Hz, coh² ≥ 0.5): green nominal, "
        "blue reflected pressure pole, red reversed polarity, black unexplained",
    )
    for i, r in enumerate(done):
        if not r["gain_informative"]:
            continue
        axes[1].vlines(i, r["microseism_g_h1"], r["microseism_g_h2"], lw=5, alpha=0.5)
        for source, style in (
            ("baseline", "ko"),
            ("trial", "C1x"),
            ("published", "g^"),
        ):
            if r.get(f"{source}_g") is not None:
                axes[1].plot(i, r[f"{source}_g"], style)
    axes[1].plot([], [], "ko", label="baseline Rayleigh g")
    axes[1].plot([], [], "C1x", label="trial Rayleigh g")
    axes[1].plot([], [], "g^", label="Aminian (2025)")
    axes[1].plot([], [], lw=5, alpha=0.5, label="microseism H1–H2 (informative)")
    axes[1].set(ylabel="pressure multiplier g", ylim=(0.3, 1.3), xticks=x)
    axes[1].set_xticklabels([r["station"] for r in done], rotation=90)
    axes[1].legend(fontsize=8, ncol=4)
    for axis in axes:
        axis.grid(True, alpha=0.3)
    figure.tight_layout()
    figure.savefig(destination, dpi=130)
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run", type=Path, help="Completed MICROSEISM_YV run folder")
    parser.add_argument("--docs", type=Path, help="Copy curated outputs into docs/")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    rows, summary = summarize(args.run, args.docs, args.show)
    print(json.dumps(summary["phase_label_counts"]))


if __name__ == "__main__":
    main()
