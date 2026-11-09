"""Figures and numbers of the whole-archive audit (docs/audit_2026-11-data).

Reads the outputs of s01-s06 in ``13_audit`` and writes PNGs to
``docs/audit_2026-11-data/figures`` and ``numbers.json`` next to the report.
Style and palette: ``compy.plotting.diagnostics`` (categorical hues in fixed
order, marker shapes as secondary encoding, one-hue ramps).
"""

import csv
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
from compy import anomaly  # noqa: E402
from compy.plotting import diagnostics as dg  # noqa: E402

ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
AUDIT = Path(os.environ.get("COMPY_AUDIT_OUT", ARCHIVE / "13_audit"))
DOC = HERE.parent
FIG = DOC / "figures"
C = dg.CATEGORICAL
MARK = ("o", "s", "^", "D")
NUM = {}


def rows(path):
    return list(csv.DictReader(open(path))) if Path(path).exists() else []


def fnum(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan


def style(ax):
    dg._style(ax)


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / name, dpi=140, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
def fig_phase():
    ph = rows(AUDIT / "phase/phase_metadata.csv")
    curves = json.loads((AUDIT / "phase/phase_curves.json").read_text())
    usable = [
        r for r in ph if abs(fnum(r["res_lo"])) < 45 and abs(fnum(r["res_hi"])) < 45
    ]
    groups = {"APG (absolute)": "APG", "DPG": "DPG", "Hydrophone": "Hydrophone"}
    fig, (a, b) = plt.subplots(
        1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.1, 1]}
    )
    rng = np.random.default_rng(0)
    out = {}
    for i, (key, name) in enumerate(groups.items()):
        sel = [r for r in usable if r["pressure_type"] == key]
        deep = [r for r in sel if fnum(r["depth_m"]) >= 1500]
        for subset, filled in ((sel, False), (deep, True)):
            y = np.array([fnum(r["res_lo"]) for r in subset])
            x = i + rng.uniform(-0.28, 0.28, len(y))
            a.scatter(x, y, s=22 if filled else 14, marker=MARK[i], facecolor=C[i] if filled else "none",
                      edgecolor=C[i], linewidth=0.9, alpha=0.9, zorder=3 if filled else 2)  # fmt: skip
        y_all = np.array([fnum(r["res_lo"]) for r in sel])
        y_deep = np.array([fnum(r["res_lo"]) for r in deep])
        out[name] = {"n": len(sel), "median_all": float(np.median(y_all)), "n_deep": len(deep),
                     "median_deep": float(np.median(y_deep)) if len(y_deep) else None,
                     "positive_fraction": float(np.mean(y_all > 0))}  # fmt: skip
        a.text(
            i,
            33,
            f"n = {len(sel)}\nmedian {np.median(y_all):+.1f}°",
            ha="center",
            va="top",
            fontsize=9,
            color=dg.INK,
        )
    a.axhline(0, color=dg.MUTED, lw=0.8)
    a.set_xticks(range(3), list(groups.values()))
    a.set_ylim(-45, 35)
    a.set_ylabel("Residual Z/P phase at the band bottom [°]")
    a.set_title("Every inverted station (filled: ≥ 1500 m)", fontsize=10, loc="left")
    style(a)
    # deep 7D Trillium Compact: APG vs DPG curves on a common axis
    fgrid = np.linspace(0.0075, 0.02, 30)
    for i, (key, name) in enumerate((("APG (absolute)", "APG"), ("DPG", "DPG"))):
        sel = [r for r in usable if r["network"] == "7D" and r["sensor_family"] == "Trillium Compact"
               and r["pressure_type"] == key and fnum(r["depth_m"]) >= 1500 and r["quality"] == "good"]  # fmt: skip
        m = []
        for r in sel:
            c = curves[r["sid"]]
            f, y = np.array(c["f"]), np.array(c["residual"])
            v = np.interp(fgrid, f, y, left=np.nan, right=np.nan)
            m.append(v)
            b.plot(1e3 * f, y, color=C[i], lw=0.6, alpha=0.35)
        m = np.array(m)
        med = np.nanmedian(m, axis=0)
        med[np.sum(np.isfinite(m), axis=0) < 4] = np.nan  # at least four stations
        b.plot(1e3 * fgrid, med, color=C[i], lw=2.2, marker=MARK[i], ms=4, markevery=5,
               label=f"7D Trillium Compact + {name} (n = {len(sel)})")  # fmt: skip
        out[f"7D_TC_deep_{name}"] = {
            "n": len(sel),
            "median_at_8mhz": float(np.nanmedian(m[:, 1])),
        }
    x = [fnum(r["res_lo"]) for r in usable if r["network"] == "7D" and r["sensor_family"] == "Trillium Compact"
         and fnum(r["depth_m"]) >= 1500 and r["quality"] == "good" and r["pressure_type"] == "APG (absolute)"]  # fmt: skip
    y = [fnum(r["res_lo"]) for r in usable if r["network"] == "7D" and r["sensor_family"] == "Trillium Compact"
         and fnum(r["depth_m"]) >= 1500 and r["quality"] == "good" and r["pressure_type"] == "DPG"]  # fmt: skip
    out["7D_TC_deep_mannwhitney_p"] = float(stats.mannwhitneyu(x, y).pvalue)
    b.axhline(0, color=dg.MUTED, lw=0.8)
    b.set_xscale("log")
    b.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    b.xaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    b.set_xlabel("Frequency [mHz]")
    b.set_ylabel("Residual phase [°]")
    b.set_title(f"Same seismometer, same margin, ≥ 1500 m: gauge decides (Mann-Whitney p = {out['7D_TC_deep_mannwhitney_p']:.0e})",
                fontsize=10, loc="left")  # fmt: skip
    b.legend(frameon=False, fontsize=9)
    style(b)
    save(fig, "phase_by_gauge.png")
    # by network for DPGs
    net = {}
    for r in usable:
        if r["pressure_type"] == "DPG":
            net.setdefault(r["network"], []).append(fnum(r["res_lo"]))
    out["dpg_by_network"] = {k: {"n": len(v), "median": float(np.median(v)), "positive": int(np.sum(np.array(v) > 0))}
                             for k, v in sorted(net.items())}  # fmt: skip
    seis = {}
    for r in usable:
        if r["pressure_type"] == "DPG":
            seis.setdefault(r["sensor_family"], []).append(fnum(r["res_lo"]))
    out["dpg_by_seismometer"] = {
        k: {"n": len(v), "median": float(np.median(v))} for k, v in seis.items()
    }
    out["n_stations"], out["n_usable"] = len(ph), len(usable)
    NUM["phase"] = out


def family(corners):
    return {"2": "SIO DPG (2 mHz)", "1.68;33.9": "LDEO/WHOI DPG (1.7, 34 mHz)",
            "2.9;7.7": "INSU DPG (2.9, 7.7 mHz)"}.get(corners)  # fmt: skip


def fig_tide():
    td = [
        r for r in rows(AUDIT / "dpg_tide/dpg_tide.csv") if r.get("tide_status") == "ok"
    ]
    dpg = [r for r in td if r["pressure_type"] == "DPG" and fnum(r["snr"]) > 10]
    fig, (a, b) = plt.subplots(1, 2, figsize=(13, 5))
    fams = ["SIO DPG (2 mHz)", "LDEO/WHOI DPG (1.7, 34 mHz)", "INSU DPG (2.9, 7.7 mHz)"]
    for i, fam in enumerate(fams):
        sel = [r for r in dpg if family(r["leak_corners_mhz"]) == fam]
        x = np.array([20 * np.log10(fnum(r["R_pred_lowest"])) for r in sel])
        y = np.array([20 * np.log10(fnum(r["R_M2"])) for r in sel])
        a.scatter(
            x,
            y,
            s=24,
            marker=MARK[i],
            color=C[i],
            edgecolor="white",
            linewidth=0.5,
            label=f"{fam}, n = {len(sel)}",
            zorder=3,
        )
    x_all = np.array([20 * np.log10(fnum(r["R_pred_lowest"])) for r in dpg])
    y_all = np.array([20 * np.log10(fnum(r["R_M2"])) for r in dpg])
    rho, p = stats.spearmanr(x_all, y_all)
    off = float(np.median(y_all - x_all))
    lim = (-15, 40)
    a.plot(
        lim,
        lim,
        color=dg.MUTED,
        lw=1,
        ls="--",
        label="metadata scaled exactly as the phase says",
    )
    a.plot(
        lim,
        np.array(lim) + off,
        color=dg.MUTED,
        lw=1,
        ls=":",
        label=f"same + {off:.1f} dB (median offset)",
    )
    for r in dpg:
        if r["sid"] in ("YV.RR40", "YV.RR52", "7B.NN04"):
            a.annotate(r["sid"], (20 * np.log10(fnum(r["R_pred_lowest"])), 20 * np.log10(fnum(r["R_M2"]))),
                       xytext=(6, -10), textcoords="offset points", fontsize=8, color=dg.INK)  # fmt: skip
    a.set_xlim(-15, 40)
    a.set_ylim(-15, 40)
    a.set_xlabel("Predicted from the compliance-band phase [dB]")
    a.set_ylabel("Measured: M2 after deconvolution / EOT20 [dB]")
    a.set_title(
        f"Two independent observables agree (Spearman {rho:.2f}, p = {p:.0e}, n = {len(dpg)})",
        fontsize=10,
        loc="left",
    )
    a.legend(frameon=False, fontsize=8, loc="upper left")
    style(a)
    # within-network: k vs R_M2
    within = {}
    for i, net in enumerate(("7B", "XE", "XO", "7D")):
        sel = [r for r in dpg if r["network"] == net]
        k = np.array([fnum(r["k_lowest"]) for r in sel])
        y = np.array([20 * np.log10(fnum(r["R_M2"])) for r in sel])
        rho_n, p_n = stats.spearmanr(k, y)
        within[net] = {"n": len(sel), "spearman": float(rho_n), "p": float(p_n)}
        b.scatter(k, y, s=22, marker=MARK[i], color=C[i], edgecolor="white", linewidth=0.5,
                  label=f"{net}: Spearman {rho_n:+.2f} (n = {len(sel)})")  # fmt: skip
    b.set_xscale("log")
    b.set_xlabel("Leak-corner factor k fitted to the phase (true / metadata)")
    b.set_ylabel("M2 after deconvolution / EOT20 [dB]")
    b.set_title(
        "Within one network, station-to-station: lower corner, larger tide",
        fontsize=10,
        loc="left",
    )
    b.axhline(0, color=dg.MUTED, lw=0.8)
    b.axvline(1, color=dg.MUTED, lw=0.8)
    b.legend(frameon=False, fontsize=8)
    style(b)
    save(fig, "dpg_tide_test.png")
    apg_like = [r for r in td if r["pressure_type"] == "Hydrophone"]
    NUM["tide"] = {"n_ok": len(td), "n_dpg_snr10": len(dpg), "spearman": float(rho), "p": float(p),
                   "median_offset_db": off, "within_network": within,
                   "median_R_M2_db_dpg": float(np.median(y_all)),
                   "positive_R_M2_fraction": float(np.mean(y_all > 0)),
                   "hydrophone_median_R_M2_db": float(np.median([20 * np.log10(fnum(r["R_M2"])) for r in apg_like])) if apg_like else None,
                   "YV": {r["sid"]: {"k": fnum(r["k_lowest"]), "pred_db": 20 * np.log10(fnum(r["R_pred_lowest"])),
                                     "meas_db": 20 * np.log10(fnum(r["R_M2"])), "snr": fnum(r["snr"])}
                          for r in td if r["network"] == "YV"}}  # fmt: skip


def fig_gate():
    st = [
        r
        for r in rows(AUDIT / "anomaly/stations.csv")
        if r.get("gate_bias_median_pct") not in (None, "")
    ]
    if not st:
        return
    fig, a = plt.subplots(figsize=(7.5, 4.6))
    nets = {"YV": 0, "XE": 1, "7B": 2}
    other = [r for r in st if r["network"] not in nets]
    a.scatter([100 * fnum(r["acceptance"]) for r in other], [-fnum(r["gate_bias_median_pct"]) for r in other],
              s=12, color=dg.MUTED, alpha=0.5, label=f"other networks (n = {len(other)})")  # fmt: skip
    for net, i in nets.items():
        sel = [r for r in st if r["network"] == net]
        a.scatter([100 * fnum(r["acceptance"]) for r in sel], [-fnum(r["gate_bias_median_pct"]) for r in sel],
                  s=26, marker=MARK[i], color=C[i], edgecolor="white", linewidth=0.5, label=f"{net} (n = {len(sel)})")  # fmt: skip
    a.axhline(0, color=dg.MUTED, lw=0.8)
    a.axhspan(1.4, 4.1, color=C[0], alpha=0.08, lw=0)
    a.text(
        98,
        2.75,
        "PR #21 Monte Carlo\nestimate for YV",
        ha="right",
        va="center",
        fontsize=8,
        color=dg.MUTED,
    )
    a.set_xlabel("Accepted windows [%]")
    a.set_ylabel("Inclusive (published) over exclusive gate − 1 [%]")
    a.set_ylim(-3, 6)
    a.set_title(
        "Coherence-gate selection bias measured on real data", fontsize=10, loc="left"
    )
    a.legend(frameon=False, fontsize=8, loc="lower right", bbox_to_anchor=(1, 0.58))
    style(a)
    save(fig, "gate_bias_measured.png")
    v = np.array([-fnum(r["gate_bias_median_pct"]) for r in st])
    NUM["gate"] = {"n": len(st), "median_pct": float(np.median(v)), "p10": float(np.percentile(v, 10)),
                   "p90": float(np.percentile(v, 90)), "abs_ge_1pct": int(np.sum(np.abs(v) >= 1)),
                   "by_network": {n: float(np.median([-fnum(r["gate_bias_median_pct"]) for r in st if r["network"] == n]))
                                  for n in sorted({r["network"] for r in st})},
                   "spearman_vs_acceptance": float(stats.spearmanr([fnum(r["acceptance"]) for r in st], v)[0])}  # fmt: skip


def fig_which_sensor():
    d = np.load(AUDIT / "level_steps/YV_relative_levels.npz")
    names = list(d["stations"])
    t = np.array([datetime.fromtimestamp(x, tz=timezone.utc) for x in d["times"]])
    panels = (
        ("RR50", "2013-07-30"),
        ("RR36", "2013-10-21"),
        ("RR52", "2013-09-14"),
        ("RR28", "2013-09-14"),
    )
    fig, axes = plt.subplots(len(panels), 1, figsize=(11, 9), sharex=True)
    import pandas as pd

    for ax, (sta, day) in zip(axes, panels):
        i = names.index(sta)
        series = (("Z_ms", "Z, 0.1–0.3 Hz"), ("P_ms", "P, 0.1–0.3 Hz"),
                  ("Z_ig", "Z, 7–15 mHz"), ("P_ig", "P, 7–15 mHz"))  # fmt: skip
        for j, (key, lab) in enumerate(series):
            x = d[key][i]
            ok = np.isfinite(x)
            if ok.sum() < 20:
                continue
            # each series relative to its own median, 7-day running median of valid bins only
            sm = pd.Series(x[ok] - np.median(x[ok]), index=pd.DatetimeIndex(t[ok]))
            sm = sm.rolling("7D", center=True, min_periods=5).median()
            ax.plot(sm.index, sm.values, color=C[j], lw=1.3, label=lab)
        ax.axvline(datetime.fromisoformat(day).replace(tzinfo=timezone.utc), color=dg.INK, lw=0.8, ls="--")
        ax.set_ylabel(f"{sta}\nre network [dB]")
        ax.set_ylim(-8, 10)
        style(ax)
    axes[0].legend(frameon=False, fontsize=8, ncol=4, loc="upper left")
    axes[0].set_title(
        "Which sensor changed? Station level minus the YV median (7-day running median)",
        fontsize=10,
        loc="left",
    )
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    save(fig, "which_sensor_yv.png")
    NUM["which_sensor_yv"] = {r["station"] + " " + r["date"]: {k: r[k] for k in r if k.endswith("_db") or k.endswith("_z")}
                              for r in rows(AUDIT / "level_steps/dated_steps.csv") if r["listed"] == "True"}  # fmt: skip


def fig_epochs():
    st = rows(AUDIT / "anomaly/stations.csv")
    rec = []
    for r in st:
        for s in json.loads(r.get("steps") or "[]"):
            rec.append(
                (r["sid"], r["tier"], r["quality"], s, fnum(r["epoch_shift_pct"]))
            )
    if not rec:
        return
    shifts = {r["sid"]: fnum(r["epoch_shift_pct"]) for r in st
              if np.isfinite(fnum(r["epoch_shift_pct"])) and int(r.get("n_steps") or 0) > 0}  # fmt: skip
    sensors = {}
    for sid, tier, q, s, sh in rec:
        sensors.setdefault(sid, set()).add(s["sensor"])
    top = sorted(shifts.items(), key=lambda kv: -abs(kv[1]))[:30]
    fig, a = plt.subplots(figsize=(8, 7.5))
    order = {"Z": 0, "P": 1, "both": 2, "neither": 3}
    y = np.arange(len(top))[::-1]
    for yy, (sid, v) in zip(y, top):
        s = sorted(sensors.get(sid, {"neither"}), key=lambda k: order[k])[0]
        a.barh(yy, v, color=C[order[s]], height=0.7)
    a.set_yticks(y, [k for k, _ in top], fontsize=8)
    a.axvline(0, color=dg.MUTED, lw=0.8)
    a.set_xlabel("Stack used by stage E minus stack before the first |Z/P| step [%]")
    a.set_xscale("symlog", linthresh=2)
    handles = [plt.Rectangle((0, 0), 1, 1, color=C[i]) for i in range(4)]
    a.legend(handles, ["vertical level stepped", "pressure level stepped", "both", "neither (ratio only)"],
             frameon=False, fontsize=8, loc="lower right")  # fmt: skip
    a.set_title("Epoch mixing: the 30 largest shifts", fontsize=10, loc="left")
    style(a)
    save(fig, "epoch_mixing.png")
    tier1 = {r["sid"] for r in st if r["tier"] == "1"}
    good = {r["sid"] for r in st if r["quality"] == "good"}
    sh = np.array(list(shifts.values()))
    NUM["epochs"] = {"stations_scanned": len(st), "stations_with_steps": len(shifts),
                     "abs_shift_ge_1pct": int(np.sum(np.abs(sh) >= 1)), "abs_shift_ge_5pct": int(np.sum(np.abs(sh) >= 5)),
                     "tier1_abs_ge_1pct": sorted(k for k, v in shifts.items() if abs(v) >= 1 and k in tier1),
                     "tier1_good_abs_ge_5pct": sorted(k for k, v in shifts.items() if abs(v) >= 5 and k in tier1 and k in good),
                     "sensor_counts": {k: sum(1 for s in sensors.values() if k in s) for k in order}}  # fmt: skip


def fig_tilt():
    tt = rows(AUDIT / "tilt/tilt_tides.csv")
    nets = sorted(
        {r["network"] for r in tt},
        key=lambda n: -sum(1 for r in tt if r["network"] == n),
    )
    nets = [n for n in nets if sum(1 for r in tt if r["network"] == n) >= 6]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    out = {}
    for ax, (col, title) in zip(axes, (("vr_M2_p", "Tilt-correction variance reduction at M2"),
                                       ("angle_M2_p", "Tilt angle at M2"), ("angle_S1_p", "Tilt angle at S1 (24 h)"))):  # fmt: skip
        frac = [
            np.mean([fnum(r[col]) < 0.01 for r in tt if r["network"] == n])
            for n in nets
        ]
        ax.bar(range(len(nets)), 100 * np.array(frac), color=C[0], width=0.65)
        ax.set_xticks(range(len(nets)), nets)
        ax.set_title(title, fontsize=10, loc="left")
        style(ax)
        out[col] = {n: float(f) for n, f in zip(nets, frac)}
        out[col + "_all"] = int(sum(fnum(r[col]) < 0.01 for r in tt))
    axes[0].set_ylabel("Stations with a line (tail p < 0.01) [%]")
    save(fig, "tilt_tides.png")
    out["n"] = len(tt)
    NUM["tilt"] = out


def fig_ranking():
    an = rows(AUDIT / "anomaly/anomalies.csv")
    if not an:
        return
    recs = [anomaly.Anomaly(r["station"], r["category"], r["kind"], r["summary"], fnum(r["significance"]),
                            fnum(r["effect_pct"]), r.get("time", "")) for r in an]  # fmt: skip
    fig = dg.plot_anomaly_ranking(recs, top=30, title="Whole archive: top 30 anomalies")
    save(fig, "anomaly_ranking.png")
    kinds = {}
    for r in an:
        kinds[r["kind"]] = kinds.get(r["kind"], 0) + 1
    NUM["anomalies"] = {
        "n": len(an),
        "stations": len({r["station"] for r in an}),
        "by_kind": kinds,
    }


def fig_multistation():
    from compy.results import load_station

    for net, ncols in (("YV", 4), ("XE", 6)):
        plan = [
            r
            for r in rows(ARCHIVE / "12_inversion_all_stations/tables/station_plan.csv")
            if r["network"] == net
        ]
        st = {}
        for p in plan:
            sub = AUDIT / "stage_c/wide_inclusive" / net / "results" / p["station"]
            if p["epoch_end"]:
                sub = sub / "epoch"
            try:
                st[p["station"]] = load_station(
                    sub, p["station"], p["wide_gain_used"] or None
                )
            except FileNotFoundError:
                continue
        if len(st) < 3:
            continue
        fig = dg.plot_multistation(
            st,
            ncols=ncols,
            title=f"{net}: compliance ± total σ (wide gates, as inverted)",
        )
        save(fig, f"{net}_multistation.png")
        fig = dg.plot_noise_models(
            st,
            kind="vertical",
            title=f"{net}: stacked vertical PSD of accepted windows",
        )
        save(fig, f"{net}_noise_vertical.png")


def main():
    for f in (
        fig_phase,
        fig_tide,
        fig_which_sensor,
        fig_tilt,
        fig_gate,
        fig_epochs,
        fig_ranking,
        fig_multistation,
    ):
        try:
            f()
            print("ok", f.__name__, flush=True)
        except Exception as exc:  # noqa: BLE001
            import traceback

            traceback.print_exc()
            print("FAILED", f.__name__, exc, flush=True)
    old = (
        json.loads((DOC / "numbers.json").read_text())
        if (DOC / "numbers.json").exists()
        else {}
    )
    old.update(NUM)
    (DOC / "numbers.json").write_text(json.dumps(old, indent=1, default=float))


if __name__ == "__main__":
    main()
