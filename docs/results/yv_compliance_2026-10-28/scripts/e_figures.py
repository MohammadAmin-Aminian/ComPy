"""Stage E: figures and tables for the YV compliance report.

Reads results/<STA>/ (Stage C), work/<STA>/clean/ (Stage B) and
work/<STA>/calibration/ (Stage D); writes figures/ and results/tables/.
Colours: reference categorical slots 1–3 (blue, orange, aqua), at most three
series per panel; the paper value is drawn in neutral grey as a reference.
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from obspy import UTCDateTime, read  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "ComPy"))
sys.path.insert(0, str(ROOT / "scripts"))
from compy_processing import spectra  # noqa: E402
from compy_numerics import FREE_AIR_GRADIENT, wavenumber  # noqa: E402
import c_compliance as C  # noqa: E402

FIG = ROOT / "figures"
TABLES = ROOT / "results" / "tables"
STATIONS = ("RR28", "RR29", "RR34", "RR36", "RR38", "RR40", "RR50", "RR52")
GROUPS = {
    "Somali Plate": ("RR28", "RR29", "RR34"),
    "SWIR": ("RR36", "RR38", "RR40"),
    "RTJ and CIR": ("RR50", "RR52"),
}
SLOTS = ("#2a78d6", "#eb6834", "#1baf7a")
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d3"
GAIN_LABEL = {"bruitfm": "BRUIT-FM gain (gains.csv)", "compy": "ComPy 2.2 calibration"}
GAIN_COLOR = {"bruitfm": SLOTS[0], "compy": SLOTS[1]}

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "lines.linewidth": 2,
    "savefig.bbox": "tight", "figure.dpi": 110,
})  # fmt: skip


def load(station, epoch=True):
    """Stable-epoch result when it exists (recommended), else full period."""
    base = ROOT / "results" / station
    if epoch and (base / "epoch" / f"{station}_summary.json").exists():
        base = base / "epoch"
    summary_path = base / f"{station}_summary.json"
    if not summary_path.exists():
        return None, {}
    summary = json.loads(summary_path.read_text())
    data = {}
    for label in summary["results"]:
        p = base / f"{station}_compliance_{label}.npz"
        if p.exists():
            data[label] = dict(np.load(p))
    return summary, data


def station_figure(station, summary, data, max_windows=300):
    label = "bruitfm" if "bruitfm" in data else next(iter(data))
    d = data[label]
    starts = d["window_starts"]
    rng = np.random.default_rng(0)
    pick = set(
        rng.choice(starts, size=min(max_windows, len(starts)), replace=False).tolist()
    )
    windows, _ = C.load_windows(station)
    if summary.get("epoch_end"):
        windows = [
            w
            for w in windows
            if w[0].stats.endtime <= UTCDateTime(summary["epoch_end"])
        ]
    acc, pre, coh = [], [], []
    for w in windows:
        if float(w[0].stats.starttime) in pick:
            f, zz, pp, _, c = spectra(w, C.NSEG)
            acc.append(zz * (2 * np.pi * f) ** 4)
            pre.append(pp * float(d["gain"]) ** 2)
            coh.append(c)
    band = d["band"]
    thr = summary["results"][label]["threshold"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), constrained_layout=True)
    for ax, values, ylabel in (
        (axes[0, 0], acc, "Vertical acceleration PSD [dB re 1 (m/s²)²/Hz]"),
        (axes[0, 1], pre, "Calibrated pressure PSD [dB re 1 Pa²/Hz]"),
    ):
        values = 10 * np.log10(np.asarray(values)[:, 1:])
        ax.semilogx(f[1:], values.T, color=SLOTS[0], alpha=0.04, lw=0.7)
        ax.semilogx(
            f[1:],
            np.median(values, axis=0),
            color=INK,
            lw=2,
            label="Median of shown windows",
        )
        shown = values[:, (f[1:] >= 0.002) & (f[1:] <= 0.3)]
        low, high = np.percentile(shown, [0.5, 99.5])
        ax.set(xlabel="Frequency [Hz]", ylabel=ylabel, xlim=(0.002, 0.3),
               ylim=(low - 5, high + 5))  # fmt: skip
        ax.legend(loc="upper right")
    ax = axes[1, 0]
    ax.semilogx(f[1:], np.asarray(coh)[:, 1:].T, color=SLOTS[0], alpha=0.04, lw=0.7)
    ax.semilogx(d["all_freqs"][1:], d["mean_coherence_all"][1:], color=INK, lw=2,
                label="Mean over accepted windows")  # fmt: skip
    ax.axhline(thr, color=MUTED, ls="--", lw=1, label=f"Threshold {thr}")
    for x in band:
        ax.axvline(x, color=MUTED, ls=":", lw=1)
    ax.set(
        xlabel="Frequency [Hz]",
        ylabel="Amplitude coherence Z–P",
        ylim=(0, 1),
        xlim=(0.002, 0.3),
    )
    ax.legend(loc="lower left")
    ax = axes[1, 1]
    for gl, dd in data.items():
        sel = (dd["freqs"] >= dd["band"][0]) & (dd["freqs"] <= dd["band"][1])
        ax.errorbar(dd["freqs"][sel], dd["compliance"][sel], yerr=dd["uncertainty"][sel],
                    color=GAIN_COLOR[gl], marker="o", ms=4, lw=2, capsize=2,
                    label=f"{GAIN_LABEL[gl]}, g = {float(dd['gain']):.3f}")  # fmt: skip
        ax.fill_between(dd["freqs"][sel], (dd["window_mean"] - dd["window_std"])[sel],
                        (dd["window_mean"] + dd["window_std"])[sel], color=GAIN_COLOR[gl],
                        alpha=0.12, lw=0)  # fmt: skip
    paper = summary.get("paper_table1")
    if paper:
        for x in paper[1]:
            ax.axvline(x, color="#9a9994", ls="--", lw=1)
    ax.set(xlabel="Frequency [Hz]", ylabel="Normalized compliance η [1/Pa]")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.legend(loc="upper left", fontsize=8.5)
    res = summary["results"][label]
    epoch = (
        f" (stable epoch to {summary['epoch_end']})" if summary.get("epoch_end") else ""
    )
    fig.suptitle(
        f"YV.{station}{epoch} — depth {summary['depth_m']:.0f} m, {res['n_windows']} accepted 35-min "
        f"windows, band {band[0]:.4f}–{band[1]:.4f} Hz (dotted); paper band dashed grey. "
        f"Shading: ±1σ window scatter; bars: Bendat & Piersol error.", fontsize=9.5,
    )  # fmt: skip
    fig.savefig(FIG / f"{station}_compliance_panels.png")
    plt.close(fig)


def network_figure(all_data):
    fig, axes = plt.subplots(
        1, 3, figsize=(15, 4.8), sharey=True, constrained_layout=True
    )
    for ax, (group, stations) in zip(axes, GROUPS.items()):
        for color, station in zip(SLOTS, stations):
            d = all_data.get(station, {}).get("bruitfm")
            if d is None:
                continue
            sel = (d["freqs"] >= d["band"][0]) & (d["freqs"] <= d["band"][1])
            ax.errorbar(d["freqs"][sel], d["compliance"][sel], yerr=d["uncertainty"][sel],
                        color=color, marker="o", ms=4, capsize=2, label=station)  # fmt: skip
            ax.annotate(station, (d["freqs"][sel][-1], d["compliance"][sel][-1]),
                        xytext=(4, 0), textcoords="offset points", color=INK, fontsize=9)  # fmt: skip
        ax.set(xscale="log", yscale="log", title=group, xlabel="Frequency [Hz]")
        ax.legend(loc="lower right")
    axes[0].set_ylabel("Normalized compliance η [1/Pa]")
    fig.suptitle(
        "YV broadband stations, BRUIT-FM gains (cf. Aminian et al. 2025, Fig. 8)"
    )
    fig.savefig(FIG / "network_compliance.png")
    plt.close(fig)


def gain_figure(rows):
    stations = [r["station"] for r in rows]
    x = np.arange(len(stations))
    fig, ax = plt.subplots(figsize=(10, 4.2), constrained_layout=True)
    series = (("bruitfm", "BRUIT-FM (gains.csv)", SLOTS[0]), ("compy", "ComPy 2.2", SLOTS[1]),
              ("paper", "Paper Table 1", SLOTS[2]))  # fmt: skip
    for i, (key, label, color) in enumerate(series):
        y = [r.get(f"gain_{key}") or np.nan for r in rows]
        err = [r.get(f"gain_{key}_std") or 0 for r in rows]
        ax.errorbar(
            x + (i - 1) * 0.2,
            y,
            yerr=err,
            fmt="o",
            color=color,
            ms=7,
            capsize=3,
            label=label,
        )
    ax.set_xticks(x, stations)
    ax.set(ylabel="DPG gain g (pressure × g)", ylim=(0.4, 1.0))
    ax.legend(ncols=3, loc="upper left")
    fig.savefig(FIG / "gain_comparison.png")
    plt.close(fig)


def monthly_figure(station, summary, data):
    d = data.get("bruitfm")
    months = summary["results"].get("bruitfm", {}).get("monthly", {})
    if d is None or len(months) < 2:
        return
    sel = (d["freqs"] >= d["band"][0]) & (d["freqs"] <= d["band"][1])
    f = d["freqs"][sel]
    fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
    for key, m in months.items():
        ax.plot(f, m["compliance"], color=SLOTS[0], alpha=0.35, lw=1)
    ax.errorbar(f, d["compliance"][sel], yerr=d["uncertainty"][sel], color=INK, lw=2,
                capsize=2, label="All accepted windows")  # fmt: skip
    ax.plot(
        [],
        [],
        color=SLOTS[0],
        alpha=0.5,
        lw=1,
        label=f"{len(months)} monthly stacks (≥20 windows)",
    )
    ax.set(xlabel="Frequency [Hz]", ylabel="Normalized compliance η [1/Pa]",
           title=f"YV.{station}: temporal stability")  # fmt: skip
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.legend()
    fig.savefig(FIG / f"{station}_monthly.png")
    plt.close(fig)


def tilt_figure(station):
    path = ROOT / "work" / station / "clean" / "tilt_hourly.csv"
    if not path.exists():
        return
    rows = list(csv.DictReader(open(path)))
    if not rows:
        return
    t = np.array([UTCDateTime(r["centre_utc"]).matplotlib_date for r in rows])
    ang = np.array([float(r["angle_deg"]) for r in rows])
    az = np.array([float(r["azimuth_deg"]) for r in rows])
    fig, axes = plt.subplots(
        2, 1, figsize=(10, 5.5), sharex=True, constrained_layout=True
    )
    axes[0].plot(t, az, ".", ms=1.5, color=SLOTS[0])
    axes[0].set(
        ylabel="Tilt azimuth [°]", title=f"YV.{station}: hourly tilt (CleanRotator)"
    )
    axes[1].plot(t, ang, ".", ms=1.5, color=SLOTS[0])
    axes[1].set(ylabel="Tilt angle [°]")
    axes[1].xaxis_date()
    fig.savefig(FIG / f"{station}_tilt.png")
    plt.close(fig)


def stability_figure():
    steps_path = TABLES / "stability_steps.json"
    steps = json.loads(steps_path.read_text()) if steps_path.exists() else {}
    fig, axes = plt.subplots(
        4, 2, figsize=(13, 11), sharex=True, constrained_layout=True
    )
    for ax, station in zip(axes.flat, STATIONS):
        path = ROOT / "results" / station / f"{station}_window_levels.npy"
        if not path.exists():
            continue
        arr = np.load(path)
        t = np.array([UTCDateTime(x).matplotlib_date for x in arr[:, 0]])
        ratio = 10 ** (arr[:, 1] - np.median(arr[:, 1]))
        ax.plot(t, ratio, ".", ms=1.5, color=SLOTS[0])
        step = steps.get(station, {})
        title = station
        if step.get("significant"):
            ax.axvline(
                UTCDateTime(step["step_time"]).matplotlib_date, color=SLOTS[1], lw=1.5
            )
            title += f": step ×{step['step_factor']:.2f} on {step['step_time'][:10]}"
        ax.set(title=title, ylim=(0.6, 2.2), yscale="log")
        ax.axhline(1, color=MUTED, lw=1)
        ax.xaxis_date()
    for ax in axes[:, 0]:
        ax.set_ylabel("|S_pz|/S_pp ÷ median")
    fig.suptitle(
        "Per-window Z/P transfer function in the compliance band (accepted windows)"
    )
    fig.savefig(FIG / "stability_all_stations.png")
    plt.close(fig)


def map_figure(summaries):
    from bruitfm.maps import PLATE, add_land
    from obspy import read_inventory

    inv = read_inventory(C.INVENTORY)
    fig, ax = plt.subplots(figsize=(8, 6.5), subplot_kw={"projection": PLATE})
    extent = (45, 75, -38, -15)
    add_land(ax, extent, scale="50m")
    for station in STATIONS:
        s = inv.select(station=station)[0][0]
        ax.plot(
            s.longitude,
            s.latitude,
            "^",
            ms=10,
            color=SLOTS[0],
            mec="white",
            transform=PLATE,
        )
        ax.annotate(station, (s.longitude, s.latitude), xytext=(5, 5), textcoords="offset points",
                    transform=PLATE, color=INK, fontsize=9)  # fmt: skip
    ax.set_title("YV broadband OBS used for compliance")
    fig.savefig(FIG / "station_map.png")
    plt.close(fig)


def main():
    FIG.mkdir(exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    all_data, rows, summaries = {}, [], {}
    for station in STATIONS:
        summary, data = load(station)
        cal = ROOT / "work" / station / "calibration" / "compy_gain.json"
        cal = json.loads(cal.read_text()) if cal.exists() else {}
        row = {"station": station,
               "gain_compy": cal.get("gain"), "gain_compy_std": cal.get("gain_std"),
               "gain_compy_events": len(cal.get("event_gains", []))}  # fmt: skip
        paper = C.PAPER_TABLE1[station]
        row.update(gain_paper=paper[0], gain_paper_std=0.05 * paper[0],
                   paper_band=f"{paper[1][0]}–{paper[1][1]}")  # fmt: skip
        tilt_figure(station)
        full_summary, full_data = load(station, epoch=False)
        if full_summary and "bruitfm" in full_data:
            d = full_data["bruitfm"]
            sel = (d["freqs"] >= d["band"][0]) & (d["freqs"] <= d["band"][1])
            row["full_period_bruitfm_median_eta"] = float(
                np.median(d["compliance"][sel])
            )
            row["full_period_bruitfm_windows"] = full_summary["results"]["bruitfm"][
                "n_windows"
            ]
        row["epoch_end"] = summary.get("epoch_end") if summary else None
        if summary is None:
            rows.append(row)
            continue
        summaries[station] = summary
        all_data[station] = data
        row["depth_m"] = summary["depth_m"]
        row["gain_bruitfm"] = summary["gains"].get("bruitfm")
        for label, res in summary["results"].items():
            if "band_hz" not in res:
                row[f"{label}_status"] = "no windows"
                continue
            d = data[label]
            sel = (d["freqs"] >= d["band"][0]) & (d["freqs"] <= d["band"][1])
            row[f"{label}_band"] = f"{res['band_hz'][0]:.4f}–{res['band_hz'][1]:.4f}"
            row[f"{label}_threshold"] = res["threshold"]
            row[f"{label}_windows"] = res["n_windows"]
            row[f"{label}_median_eta"] = float(np.median(d["compliance"][sel]))
            row[f"{label}_median_rel_err"] = res["median_relative_uncertainty"]
            row[f"{label}_median_phase"] = res["median_phase_deg_in_band"]
            eta = d["compliance"][sel]
            row[f"{label}_median_rel_sem"] = float(np.median(
                d["window_std"][sel] / np.sqrt(res["n_windows"]) / eta))  # fmt: skip
            row[f"{label}_median_rel_window_std"] = float(
                np.median(d["window_std"][sel] / eta)
            )
            months = [
                np.array(m["compliance"]) for m in res.get("monthly", {}).values()
            ]
            if len(months) >= 3:
                row[f"{label}_monthly_rel_spread"] = float(
                    np.median(np.std(months, axis=0) / eta)
                )
                row[f"{label}_months"] = len(months)
            # Size of fix D1 at this station: (1 - mean coherence) x gravity term / eta.
            f = d["freqs"][sel]
            k = wavenumber(2 * np.pi * f, float(d["depth"]))
            omega2 = (2 * np.pi * f) ** 2
            allf = d["all_freqs"]
            ratio = np.interp(f, allf, np.abs(d["S_pz"]) / d["S_pp"])
            measured = k * omega2 * ratio / (omega2 + FREE_AIR_GRADIENT)
            coh = np.interp(f, allf, d["mean_coherence_all"])
            row[f"{label}_d1_bias_percent"] = float(
                np.median((1 - coh) * (eta - measured) / eta) * 100
            )
            row[f"{label}_gravity_term_percent"] = float(
                np.median((eta - measured) / eta) * 100
            )
        try:
            station_figure(station, summary, data)
        except Exception as exc:  # keep the other figures
            print(f"{station} panel figure failed: {exc!r}")
        monthly_figure(station, summary, data)
        rows.append(row)
    network_figure(all_data)
    gain_figure(rows)
    stability_figure()
    try:
        map_figure(summaries)
    except Exception as exc:
        print(f"map failed: {exc!r}")
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k != "station", k))
    with open(TABLES / "compliance_summary.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, keys)
        writer.writeheader()
        writer.writerows(rows)
    (TABLES / "compliance_summary.json").write_text(
        json.dumps(rows, indent=1, default=str)
    )
    print(json.dumps(rows, indent=1, default=str))


if __name__ == "__main__":
    main()
