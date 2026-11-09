"""Figures and numbers of the research part (RESEARCH.md), from r01-r06 outputs."""

import csv
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from compy.plotting import diagnostics as dg  # noqa: E402

ARCHIVE = Path(os.environ.get("COMPY_ARCHIVE", Path.home() / "Desktop/ComPy"))
R = ARCHIVE / "13_audit/research"
DOC = HERE.parent
FIG = DOC / "figures"
C, MARK = dg.CATEGORICAL, ("o", "s", "^", "D", "v")
NUM = {}


def save(fig, name):
    fig.savefig(FIG / name, dpi=140, bbox_inches="tight")
    plt.close(fig)


def dates(t):
    return [datetime.fromtimestamp(x, tz=timezone.utc) for x in t]


def fig_tide_reference():
    m = pd.read_csv(R / "monthly_tide.csv")
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    a = axes[0, 0]
    for i, (sid, day) in enumerate(
        (
            ("YV.RR50", "2013-07-30"),
            ("YV.RR52", "2013-09-19"),
            ("YV.RR28", "2013-09-16"),
            ("YV.RR40", "2013-09-14"),
        )
    ):
        g = m[m.sid == sid].sort_values("t_mid")
        a.plot(
            dates(g.t_mid),
            100 * (g.M2 / g.M2.median() - 1),
            color=C[i],
            marker=MARK[i],
            ms=4,
            lw=1.3,
            label=sid,
        )
    a.axvline(
        datetime(2013, 7, 30, tzinfo=timezone.utc), color=dg.MUTED, ls="--", lw=0.8
    )
    a.axvline(
        datetime(2013, 9, 14, tzinfo=timezone.utc), color=dg.MUTED, ls=":", lw=0.8
    )
    a.text(
        datetime(2013, 7, 31, tzinfo=timezone.utc),
        52,
        "RR50 Z step",
        fontsize=8,
        color=dg.MUTED,
    )
    a.text(
        datetime(2013, 9, 15, tzinfo=timezone.utc),
        44,
        "Sept. |Z/P|\nsteps",
        fontsize=8,
        color=dg.MUTED,
    )
    a.set_ylabel("M2 in pressure counts, 30-day fits [% of median]")
    a.set_ylim(-35, 60)
    a.set_title("YV: the tide as a pressure-only reference", fontsize=10, loc="left")
    a.legend(frameon=False, fontsize=8, ncol=2)
    a.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    dg._style(a)
    a = axes[0, 1]
    for i, sid in enumerate(("8J.UP38", "8J.UP45", "8J.UP44", "8J.UP04", "8J.UP16")):
        g = m[m.sid == sid].sort_values("t_mid")
        a.plot(
            dates(g.t_mid),
            20 * np.log10(g.M2 / g.M2.median()),
            color=C[i],
            marker=MARK[i],
            ms=4,
            lw=1.3,
            label=sid,
        )
    a.set_ylabel("M2 in hydrophone counts [dB re median]")
    a.set_title(
        "8J: low-frequency hydrophone response changes during deployment",
        fontsize=10,
        loc="left",
    )
    a.legend(frameon=False, fontsize=8, ncol=2)
    a.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    dg._style(a)
    # 8J: monthly M2 vs monthly |Z/P| and phase, pooled (station-demeaned)
    tp = pd.read_csv(R / "tide_phase_monthly.csv")
    rows = []
    for sid in tp[tp.sid.str.startswith("8J")].sid:
        net, sta = sid.split(".")
        path = list(
            (ARCHIVE / f"13_audit/stage_c/wide_inclusive/{net}/results/{sta}").glob(
                f"{sta}_window_table_*.npz"
            )
        )
        if not path:
            continue
        w = np.load(path[0])
        ok = w["coherence"] > 0.8
        t, tf, ph = w["time"][ok], w["log10_tf"][ok], w["phase_deg"][ok]
        ph = (ph - np.median(ph) + 180) % 360 - 180
        g = m[m.sid == sid]
        for tm, m2 in zip(g.t_mid, g.M2):
            s = np.abs(t - tm) < 15 * 86400
            if s.sum() >= 20:
                rows.append(
                    (
                        sid,
                        20 * np.log10(m2),
                        100 * (10 ** np.median(tf[s]) - 1),
                        np.median(ph[s]),
                    )
                )
    d = pd.DataFrame(rows, columns=["sid", "m2_db", "tf", "ph"])
    for col in ("m2_db", "ph"):
        d[col] = d[col] - d.groupby("sid")[col].transform("median")
    d["tf_pct"] = 100 * (
        np.log(1 + d["tf"] / 100)
        - d.assign(x=np.log(1 + d["tf"] / 100)).groupby("sid")["x"].transform("median")
    )
    for a, col, lab in ((axes[1, 0], "tf_pct", "|Z/P| in the compliance band [% re station median]"),
                        (axes[1, 1], "ph", "Z/P phase in the compliance band [° re station median]")):  # fmt: skip
        a.scatter(d.m2_db, d[col], s=10, color=C[0], alpha=0.45, edgecolor="none")
        rho, p = stats.spearmanr(d.m2_db, d[col])
        a.set_xlabel("M2 in hydrophone counts [dB re station median]")
        a.set_ylabel(lab)
        a.set_title(
            f"8J, {d.sid.nunique()} stations, monthly: Spearman {rho:+.2f} (n = {len(d)})",
            fontsize=10,
            loc="left",
        )
        a.axhline(0, color=dg.MUTED, lw=0.8)
        a.axvline(0, color=dg.MUTED, lw=0.8)
        dg._style(a)
        NUM[f"8J_monthly_M2_vs_{col}"] = {
            "spearman": float(rho),
            "p": float(p),
            "n": len(d),
        }
    axes[1, 0].set_ylim(-80, 80)
    axes[1, 1].set_ylim(-40, 40)
    save(fig, "research_tide_reference.png")
    ts = pd.read_csv(R / "tide_at_steps.csv")
    NUM["tide_at_steps"] = {
        "n": len(ts),
        "verdicts": ts.verdict.value_counts().to_dict(),
    }
    NUM["8J_phase_tracks_M2"] = int(
        ((tp.sid.str.startswith("8J")) & (tp.p_phase < 0.01)).sum()
    )
    NUM["8J_stations_tested"] = int(tp.sid.str.startswith("8J").sum())


def fig_bracket():
    b = pd.read_csv(R / "response_bracket.csv").dropna(subset=["k_tide"])
    cons = pd.read_csv(R / "tide_constituents.csv")
    fam = {
        "2": ("SIO DPG", 0),
        "1.68;33.9": ("LDEO/WHOI DPG", 1),
        "2.9;7.7": ("INSU DPG", 2),
    }
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    a = axes[0]
    n = cons.n_slope.dropna()
    a.hist(n, bins=np.linspace(-1.5, 1.5, 31), color=C[0], edgecolor="white")
    a.axvline(0, color=dg.INK, lw=1, ls="--")
    a.axvline(1, color=dg.MUTED, lw=1, ls=":")
    a.text(0.03, a.get_ylim()[1] * 0.92, "scaled leak corner\n(n = 0)", fontsize=8)
    a.text(
        0.97,
        a.get_ylim()[1] * 0.6,
        "one missing\nhigh-pass\n(n = 1)",
        fontsize=8,
        color=dg.MUTED,
        ha="right",
    )
    a.set_xlabel("Slope n of R(f) across O1, K1, M2, S2")
    a.set_ylabel("DPGs")
    a.set_title(
        f"Tidal ratio is flat in frequency (median n = {np.median(n):+.2f})",
        fontsize=10,
        loc="left",
    )
    dg._style(a)
    a = axes[1]
    for key, (lab, i) in fam.items():
        s = b[b.family == key]
        a.scatter(
            s.res_lo_from_tide,
            s.res_lo_measured,
            s=22,
            marker=MARK[i],
            color=C[i],
            edgecolor="white",
            lw=0.5,
            label=f"{lab} (n = {len(s)})",
        )
    lim = (-20, 30)
    a.plot(lim, lim, color=dg.MUTED, ls="--", lw=1)
    rho, p = stats.spearmanr(b.res_lo_from_tide, b.res_lo_measured)
    a.set_xlim(lim)
    a.set_ylim(lim)
    a.set_xlabel("Band-bottom phase predicted from the M2 tide alone [°]")
    a.set_ylabel("Band-bottom phase measured in the compliance stack [°]")
    a.set_title(
        f"Tide predicts the compliance phase (Spearman {rho:.2f})",
        fontsize=10,
        loc="left",
    )
    a.legend(frameon=False, fontsize=8)
    dg._style(a)
    a = axes[2]
    for key, (lab, i) in fam.items():
        s = b[b.family == key]
        a.scatter(
            s.bias_lo_phase_pct,
            s.bias_lo_tide_pct,
            s=22,
            marker=MARK[i],
            color=C[i],
            edgecolor="white",
            lw=0.5,
            label=lab,
        )
    lim = (-15, 12)
    a.plot(lim, lim, color=dg.MUTED, ls="--", lw=1)
    a.set_xlim(lim)
    a.set_ylim(lim)
    a.set_xlabel("Band-bottom bias from the phase [%]")
    a.set_ylabel("Band-bottom bias from the tide [%]")
    a.set_title("Two estimates of the response bias", fontsize=10, loc="left")
    dg._style(a)
    save(fig, "research_response_bracket.png")
    NUM["constituent_slope"] = {
        "n": int(len(n)),
        "median": float(np.median(n)),
        "p25": float(np.percentile(n, 25)),
        "p75": float(np.percentile(n, 75)),
    }
    NUM["bracket"] = {"n": len(b), "spearman_phase": float(rho), "p": float(p),
                      "median_phase_offset_deg": float(np.median(b.res_lo_measured - b.res_lo_from_tide)),
                      "bias_tide_median": float(b.bias_lo_tide_pct.median()), "bias_tide_p10": float(b.bias_lo_tide_pct.quantile(0.1)),
                      "bias_tide_p90": float(b.bias_lo_tide_pct.quantile(0.9)), "bias_phase_median": float(b.bias_lo_phase_pct.median()),
                      "k_tide_median": float(b.k_tide.median()), "k_phase_median": float(b.k_phase.median())}  # fmt: skip


def fig_events():
    e = pd.read_csv(R / "earthquake_epochs.csv")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for a, (col, lab) in zip(
        axes[:2],
        (("tilt_angle_z", "Tilt angle change"), ("zp_level_z", "|Z/P| change")),
    ):
        s = e.dropna(subset=[col])
        a.scatter(
            s.shaking,
            s[col].clip(-15, 15),
            s=10,
            color=C[0],
            alpha=0.5,
            edgecolor="none",
        )
        a.axhspan(-3, 3, color=dg.GRID, alpha=0.6, lw=0)
        a.set_xlabel("Shaking proxy M − 1.66 log10(distance km)")
        a.set_ylabel(f"{lab}: robust z vs random times (clipped ±15)")
        p = e[col.replace("_z", "_p")]
        frac = float(np.mean(p.dropna() < 0.01))
        a.set_title(
            f"{lab}: {100 * frac:.1f} % of {p.notna().sum()} pairs at p < 0.01",
            fontsize=10,
            loc="left",
        )
        dg._style(a)
        NUM[f"quake_{col}"] = {"pairs": int(p.notna().sum()), "frac_p01": frac,
                               "spearman_abs_z_vs_shaking": float(stats.spearmanr(s.shaking, s[col].abs())[0])}  # fmt: skip
    # 7K bottom-current event
    a = axes[2]
    for i, sta in enumerate(("AN1", "AN3", "AN6")):
        d = np.load(ARCHIVE / f"04_ppsd/7K/{sta}/7K.{sta}_cleaned.npz")
        f = d["freqs"]
        b = (f >= 0.005) & (f <= 0.035)
        t = pd.to_datetime(d["1_times"], unit="s")
        h = np.median(d["1_psd_db"][:, b], 1)
        sel = (t >= "2019-05-20") & (t <= "2019-06-01")
        a.plot(t[sel], h[sel] - np.median(h), color=C[i], lw=1.2, label=f"7K.{sta} H1")
        tz = pd.to_datetime(d["Z_times"], unit="s")
        z = np.median(d["Z_psd_db"][:, b], 1)
        sz = (tz >= "2019-05-20") & (tz <= "2019-06-01")
        a.plot(tz[sz], z[sz] - np.median(z), color=C[i], lw=1, ls=":")
    a.set_ylabel("Level re station median, 5–35 mHz [dB] (solid H1, dotted Z)")
    a.set_title(
        "7K, May 2019: array-wide bottom-current event", fontsize=10, loc="left"
    )
    a.legend(frameon=False, fontsize=8)
    a.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    dg._style(a)
    save(fig, "research_events.png")
    ep = list(csv.DictReader(open(R / "coherence/episodes.csv")))
    hole = pd.read_csv(R / "coherence/hole.csv")
    NUM["coherence"] = {"stations": len(pd.read_csv(R / "coherence/stations.csv")), "episode_bins": len(ep),
                        "hole_stations": len(hole), "hole_min_median": float(hole.hole_min.median()),
                        "hole_min_p99": float(hole.hole_min.quantile(0.99)), "hole_filled": hole[hole.hole_min > 0.5].sid.tolist()}  # fmt: skip


def corrected_inversion():
    path = R / "inversion_corrected/comparison.csv"
    if path.exists():
        NUM["corrected_inversion_yv"] = list(csv.DictReader(open(path)))


def main():
    for f in (fig_tide_reference, fig_bracket, fig_events, corrected_inversion):
        f()
        print("ok", f.__name__, flush=True)
    old = json.loads((DOC / "numbers.json").read_text())
    old["research"] = NUM
    (DOC / "numbers.json").write_text(json.dumps(old, indent=1, default=float))


if __name__ == "__main__":
    main()
