"""Stage F: figures and tables in the style of Aminian et al. (2025, GJI).

Usage: python f_figures.py [STATION_FOR_PROCESSING_FIGURES]

Writes to figures/:
  fig01_map.png          station map on bathymetry (paper Fig. 1)
  fig03_psd_stages.png   vertical PSDs at the processing stages (Fig. 3)
  fig04_tilt.png         hourly tilt angle and azimuth (Fig. 4)
  fig05_calibration.png  Rayleigh-wave DPG calibration (Fig. 5)
  fig06_coherogram.png   one month of pressure, vertical and coherence (Fig. 6)
  fig07_windows.png      selected windows: PSDs, coherence, compliance (Fig. 7)
  fig08_compliance.png   normalized compliance of all stations (Fig. 8)
  fig09_heatmaps.png     Vs posterior heat maps (Fig. 9)
  fig10_profiles.png     median Vs by sub-array and serpentinization (Fig. 10)
  figS_misfit.png        misfit traces of all chains
  <STA>_inversion.png    data fit, misfit and heat map of one station
and results/tables/table1.{csv,md} (paper Table 1) and inversion.{csv,md}.

Depths are true metres (docs/INVERSION.md, F1); serpentinization uses the
Fig. 10 caption relation through compy_inversion.serpentinization (F2) and
is computed per posterior draw, so its interval is shown.
"""

import csv
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from obspy import UTCDateTime, read, read_inventory  # noqa: E402

from config import FIGURES, INVENTORY, NETWORK, RESULTS, ROOT, STATIONS, WORK  # noqa: E402
from e_invert import build_model  # noqa: E402

import compy_inversion as ci  # noqa: E402
from compy_processing import spectra  # noqa: E402

VIRIDIS = plt.get_cmap("viridis")
CMAP = LinearSegmentedColormap.from_list(
    "aminian2025", list(zip([0, 0.4, 1], [VIRIDIS(0.0), VIRIDIS(0.5), VIRIDIS(1.0)]))
)  # fmt: skip
Z = np.arange(0.0, 8000.0 + 1, 50.0)
VS_BINS = np.linspace(0, 5000, 101)
NSEG = 1024
plt.rcParams.update(
    {"font.family": "serif", "font.size": 12, "axes.grid": True, "grid.alpha": 0.35}
)
TABLES = RESULTS / "tables"


#: Station groups for Figs. 8 and 10 (config.GROUPS: prefix -> label); default: first two letters.
GROUPS = getattr(__import__("config"), "GROUPS", {
    "WW": "West (WW)", "WC": "West-central (WC)", "CC": "Central (CC)", "EC": "East-central (EC)", "EE": "East (EE)",
})  # fmt: skip


def sub_array(station):
    return GROUPS.get(station[:2], station[:2])


def coords():
    inv = read_inventory(str(INVENTORY))
    return {s.code: (s.latitude, s.longitude, -s.elevation) for s in inv[0]}


def title(station, c):
    lat, lon, _ = c[station]
    return f"{station} ({abs(lat):.1f}°{'S' if lat < 0 else 'N'} - {abs(lon):.1f}°{'W' if lon < 0 else 'E'})"


# --- loading -------------------------------------------------------------------


def compliance_summary(station):
    path = RESULTS / station / f"{station}_summary.json"
    return json.loads(path.read_text()) if path.exists() else None


def load_posterior(station, variant="final"):
    path = RESULTS / station / "inversion" / f"{station}_{variant}.npz"
    if not path.exists():
        return None, None
    d = np.load(path)
    # Same model as e_invert.VARIANTS: free crustal Vp/Vs except for "main".
    kwargs = {} if variant == "main" else dict(vp_vs_range=(1.7, 2.4))
    model, column = build_model(station, **kwargs)
    data = ci.ComplianceData(d["freqs"], d["eta"], float(d["depth"]), d["covariance"], station=station,
                             groups=d["groups"], group_names=tuple(d["group_names"]))  # fmt: skip
    names = list(d["names"])
    post = ci.Posterior(data=data, model=model, prior_only=False, theta=d["theta"], alpha=d["alpha"],
                        loglik=d["loglik"], roughness=d["roughness"], pred=d["pred"],
                        acceptance=np.zeros(1), acceptance_dr=np.zeros(1), settings={},
                        trace_chi2=d["trace_chi2"], names_extra=names[len(model.names):])  # fmt: skip
    return post, column


MIN_WINDOWS = (
    500  # Aminian et al. (2025, section 4.2): at least 500 windows per station
)


def n_windows(station):
    summ = compliance_summary(station)
    if not summ:
        return 0
    return max((r.get("n_windows", 0) for r in summ["results"].values()), default=0)


def usable(station):
    """Compliance from >= 500 windows and a converged final run that fits (chi2 < 2N, R-hat < 1.1)."""
    path = RESULTS / station / "inversion" / f"{station}_final.json"
    if not path.exists() or n_windows(station) < MIN_WINDOWS:
        return False
    j = json.loads(path.read_text())
    return (
        j["max_rhat"] is not None
        and j["max_rhat"] < 1.1
        and j["chi2"]["median"] < 2 * j["n_data"]
    )


# --- CRUST1.0 helpers -------------------------------------------------------------


def crust1_steps(column, zmax=8000.0):
    m = ci.crust1_model(column)
    tops = np.concatenate([[0.0], np.cumsum(m[:-1, 0])])
    z = np.repeat(np.append(tops, zmax), 2)[1:-1]
    return np.repeat(m[:, 3], 2), z


def crust1_vs_at(column, z):
    m = ci.crust1_model(column)
    tops = np.concatenate([[0.0], np.cumsum(m[:-1, 0])])
    return m[np.searchsorted(tops, z, side="right") - 1, 3]


# --- Fig. 1: map ----------------------------------------------------------------------


def fig01_map(c, status):
    from scipy.io import netcdf_file

    lats = [v[0] for v in c.values()]
    lons = [v[1] for v in c.values()]
    extent = (min(lons) - 0.6, max(lons) + 0.6, min(lats) - 0.6, max(lats) + 0.6)
    fig, ax = plt.subplots(figsize=(10, 8.5))
    bathy = (
        ROOT / "metadata" / "bathy" / f"etopo_{NETWORK.lower()}.nc"
    )  # optional ETOPO subset
    if bathy.exists():
        nc = netcdf_file(str(bathy), "r", mmap=False)
        lat, lon, z = (
            nc.variables["latitude"][:],
            nc.variables["longitude"][:],
            nc.variables["altitude"][:],
        )
        im = ax.pcolormesh(lon, lat, z / 1000, cmap="Blues_r", shading="auto", vmin=np.percentile(z, 2) / 1000,
                           vmax=np.percentile(z, 98) / 1000)  # fmt: skip
        plt.colorbar(im, ax=ax, label="Elevation [km]", fraction=0.04)
    for s, (la, lo, _) in c.items():
        st = status.get(s, "not processed")
        style = {"inverted": ("^", "red"), "compliance only": ("^", "orange"),
                 "not processed": ("v", "0.5")}[st]  # fmt: skip
        ax.plot(lo, la, style[0], color=style[1], ms=11, mec="black", mew=0.8)
        ax.text(lo + 0.08, la + 0.08, s, fontsize=9)
    for st, (mk, col) in {"inverted": ("^", "red"), "compliance only": ("^", "orange"),
                          "not processed": ("v", "0.5")}.items():  # fmt: skip
        ax.plot([], [], mk, color=col, mec="black", ms=10, label=st)
    ax.legend(loc="lower left", fontsize=10)
    ax.set(
        xlabel="Longitude [°]", ylabel="Latitude [°]", xlim=extent[:2], ylim=extent[2:]
    )
    ax.set_aspect("equal")
    ax.set_title(f"{NETWORK} broadband OBS: compliance status")
    regional = ROOT / "metadata" / "bathy" / "etopo_region.nc"  # optional inset
    if regional.exists():
        nc2 = netcdf_file(str(regional), "r", mmap=False)
        ins = ax.inset_axes([0.66, 0.66, 0.33, 0.33])
        ins.contourf(nc2.variables["longitude"][:], nc2.variables["latitude"][:], nc2.variables["altitude"][:],
                     levels=[-12000, 0, 9000], colors=["#cfe3f3", "#8c8c8c"])  # fmt: skip
        ins.add_patch(plt.Rectangle((extent[0], extent[2]), extent[1] - extent[0], extent[3] - extent[2],
                                    fill=False, color="red", lw=1.5))  # fmt: skip
        ins.set_xticks([])
        ins.set_yticks([])
        ins.grid(False)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig01_map.png", dpi=150)
    plt.close(fig)


# --- Fig. 3: PSD stages ------------------------------------------------------------------------


def _psd_windows(paths, channel, starts, width=3600):
    out = []
    for path in paths:
        st = read(str(path)).select(channel=channel)
        for tr in st:
            for t0 in starts:
                if tr.stats.starttime <= t0 and t0 + width <= tr.stats.endtime:
                    seg = tr.slice(t0, t0 + width - tr.stats.delta).data.astype(float)
                    if len(seg) >= 4000 and np.all(np.isfinite(seg)):
                        f, p = _welch(seg)
                        out.append(p)
    return f if out else None, np.array(out)


def _welch(x, fs=1.25):
    from scipy.signal import welch

    return welch(x - np.mean(x), fs=fs, nperseg=NSEG, noverlap=NSEG // 2)


def fig03_psd(station):
    clean_dir = WORK / station / "clean"
    rng = np.random.default_rng(0)
    phys = sorted((WORK / station / "phys").glob(f"{station}.*.mseed"))
    phys = [
        phys[i]
        for i in sorted(rng.choice(len(phys), min(60, len(phys)), replace=False))
    ]
    starts = []
    for p in phys:
        d = read(str(p), headonly=True)[0].stats.starttime
        day = UTCDateTime(d.year, d.month, d.day)
        starts += [day + h * 3600 for h in range(0, 24, 3)]
    f, raw = _psd_windows(phys, "BHZ", starts)
    f2, cleaned = _psd_windows(
        sorted(clean_dir.glob(f"{station}_seg*.mseed")), "BHZ", starts
    )
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    acc = (2 * np.pi * f) ** 4
    for ax, curves, label in (
        (axes[0], raw, "Response removed (Z)"),
        (axes[0], cleaned, "Cleaned (Z′)"),
    ):
        db = 10 * np.log10(curves * acc)
        ax.semilogx(f, np.median(db, axis=0), lw=2.5, label=label)
        ax.fill_between(f, *np.percentile(db, [10, 90], axis=0), alpha=0.2)
    for k, nlnm in enumerate(_nlnm(f)):
        axes[0].semilogx(
            *nlnm, color="0.3", ls="--", lw=1, label="NLNM / NHNM" if k == 0 else None
        )
    axes[0].axvspan(
        0.005, 0.015, color="palegreen", alpha=0.35, label="compliance band"
    )
    axes[0].set(xlim=(0.002, 0.6), ylim=(-200, -90), xlabel="Frequency [Hz]",
                ylabel="PSD [dB re 1 (m/s²)²/Hz]", title=f"(a) {station} vertical, median and 10–90 %")  # fmt: skip
    axes[0].legend(fontsize=10, loc="upper right")
    _, praw = _psd_windows(phys, "BDH", starts)
    gain = json.loads(
        (WORK / station / "calibration" / "compy_gain.json").read_text()
    ).get("gain", 1.0)
    db = 10 * np.log10(praw * gain**2)
    axes[1].semilogx(
        f,
        np.median(db, axis=0),
        lw=2.5,
        color="tab:purple",
        label=f"Pressure × gain ({gain:.3f})",
    )
    axes[1].fill_between(
        f, *np.percentile(db, [10, 90], axis=0), color="tab:purple", alpha=0.2
    )
    axes[1].axvspan(0.005, 0.015, color="palegreen", alpha=0.35)
    axes[1].set(xlim=(0.002, 0.6), xlabel="Frequency [Hz]", ylabel="PSD [dB re 1 Pa²/Hz]",
                title="(b) Calibrated pressure")  # fmt: skip
    axes[1].legend(fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig03_psd_stages.png", dpi=140)
    plt.close(fig)


def _nlnm(f):
    from obspy.signal.spectral_estimation import get_nhnm, get_nlnm

    out = []
    for fn in (get_nlnm, get_nhnm):
        p, db = fn()
        out.append((1 / p, db))
    return out


# --- Fig. 4: tilt ---------------------------------------------------------------------------------


def fig04_tilt(station):
    rows = list(csv.DictReader(open(WORK / station / "clean" / "tilt_hourly.csv")))
    t = np.array([UTCDateTime(r["centre_utc"]).matplotlib_date for r in rows])
    ang = np.array([float(r["angle_deg"]) for r in rows])
    az = np.array([float(r["azimuth_deg"]) for r in rows])
    vr = np.array([float(r["variance_reduction"]) for r in rows])
    fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
    for ax, y, label in (
        (axes[0], az, "Azimuth [°]"),
        (axes[1], ang, "Incident angle [°]"),
    ):
        sc = ax.scatter(
            t, y, c=np.log10(np.clip(vr, 1e-3, None)), s=4, cmap="viridis_r"
        )
        ax.set_ylabel(label)
    axes[1].set_yscale("log")
    plt.colorbar(sc, ax=axes, label="log10 variance reduction")
    axes[1].xaxis_date()
    axes[0].set_title(f"(a, b) {station}: hourly tilt minimizing vertical variance")
    fig.savefig(FIGURES / "fig04_tilt.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


# --- Fig. 5: calibration -----------------------------------------------------------------------


def fig05_calibration(station):
    path = WORK / station / "calibration" / f"{station}_compy_calibration.npz"
    if not path.exists():
        return
    d = np.load(path)
    f = d["frequencies"]
    band = (f >= d["band"][0] - 1e-9) & (f <= d["band"][1] + 1e-9)
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
    positive = f > 0  # skip a DC bin; stage D's grid starts in the band
    for coh in d["coherences"]:
        axes[0].semilogx(f[positive], coh[positive], color="0.6", lw=0.8)
    axes[0].semilogx(
        f[positive],
        np.median(d["coherences"], axis=0)[positive],
        color="blue",
        lw=2.5,
        label="Median",
    )
    for x in d["band"]:
        axes[0].axvline(x, color="red", ls="--")
    axes[0].set(xlim=(0.005, 0.1), ylim=(0, 1), xlabel="Frequency [Hz]",
                ylabel="Magnitude-squared coherence",
                title=f"(a) {station}: P–Z coherence, {len(d['coherences'])} events")  # fmt: skip
    axes[0].legend()
    axes[1].semilogx(
        f[band],
        d["measured_ratio"][band],
        color="green",
        lw=2.5,
        label="Measured P/a ÷ ρh",
    )
    axes[1].semilogx(
        f[band], d["theoretical_ratio"], color="purple", lw=2.5, label="Theoretical"
    )
    axes[1].semilogx(f[band], d["measured_ratio"][band] * d["gain"], "k--", lw=2.5,
                     label=f"Corrected (g = {float(d['gain']):.3f} ± {np.std(d['event_gains'], ddof=1):.3f})")  # fmt: skip
    axes[1].set(
        xlabel="Frequency [Hz]",
        ylabel="Normalized spectral ratio",
        title="(b) Spectral ratio",
    )
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "fig05_calibration.png", dpi=140)
    plt.close(fig)


# --- Figs. 6 and 7: windows ---------------------------------------------------------------------


def _clean_windows(station, t_start, t_end, width=2100.0):
    out = []
    for path in sorted((WORK / station / "clean").glob(f"{station}_seg*.mseed")):
        st = read(str(path))
        st = st.select(channel="BHZ") + st.select(channel="BDH")
        if st[0].stats.endtime < t_start or st[0].stats.starttime > t_end:
            continue
        t0 = max(st[0].stats.starttime, t_start)
        while t0 + width <= min(st[0].stats.endtime, t_end):
            w = st.slice(t0, t0 + width - 0.8).copy()
            if len(w) == 2 and all(np.all(np.isfinite(tr.data)) for tr in w):
                for tr in w:
                    tr.data = tr.data.astype(float)
                out.append(w)
            t0 += width
    return out


def fig06_coherogram(station, summ):
    r = summ["results"]["compy"]
    band = r["band_hz"]
    months = sorted(r["monthly"], key=lambda k: -r["monthly"][k]["n"])
    year, month = map(int, months[0].split("-"))
    t0 = UTCDateTime(year, month, 1)
    t1 = UTCDateTime(year + (month == 12), month % 12 + 1, 1)
    windows = _clean_windows(station, t0, t1)
    gain = r["gain"]
    rows = []
    for w in windows:
        f, zz, pp, zp, coh = spectra(w, NSEG)
        rows.append((w[0].stats.starttime.matplotlib_date, 10 * np.log10(pp * gain**2),
                     10 * np.log10(zz * (2 * np.pi * f) ** 4), coh))  # fmt: skip
    t = np.array([r_[0] for r_ in rows])
    P, A, C = (np.array([r_[i] for r_ in rows]).T for i in (1, 2, 3))
    sel = (f >= band[0]) & (f <= band[1])
    fig, axes = plt.subplots(
        3, 2, figsize=(16, 11), sharex="col", gridspec_kw=dict(width_ratios=[4, 1])
    )
    for k, (img, lab, lim) in enumerate(((P, "Pressure [dB re Pa²/Hz]", (0, 60)),
                                         (A, "Vertical acc. [dB re (m/s²)²/Hz]", (-190, -130)),
                                         (C, "Coherence", (0, 1)))):  # fmt: skip
        ax = axes[k, 0]
        m = ax.pcolormesh(
            t, f, img, shading="auto", vmin=lim[0], vmax=lim[1], cmap="turbo"
        )
        ax.set_yscale("log")
        ax.set_ylim(0.002, 0.1)
        for x in band:
            ax.axhline(x, color="white", ls="--")
        plt.colorbar(m, ax=ax, label=lab, fraction=0.03)
        ax.set_ylabel("Frequency [Hz]")
        val = np.mean(img[sel], axis=0) if k < 2 else np.median(img[sel], axis=0)
        axes[k, 1].plot(val, t, lw=1)
        if k == 2:
            axes[k, 1].axvline(r["threshold"], color="red", ls="--")
            axes[k, 1].fill_betweenx(
                t, 0, 1, where=val >= r["threshold"], color="green", alpha=0.25
            )
            axes[k, 1].set_xlim(0, 1)
        axes[k, 1].set_xlabel("band " + ("mean" if k < 2 else "median"))
    axes[2, 0].xaxis_date()
    axes[0, 0].set_title(
        f"{station}, {months[0]}: spectrograms of 35-min windows (a, c) and coherogram (e)"
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig06_coherogram.png", dpi=130)
    plt.close(fig)


def fig07_windows(station, summ):
    d = np.load(RESULTS / station / f"{station}_compliance_compy.npz")
    starts = d["window_starts"]
    rng = np.random.default_rng(1)
    pick = sorted(rng.choice(starts, min(150, len(starts)), replace=False))
    gain, depth, band = float(d["gain"]), float(d["depth"]), d["band"]
    from compy_processing import gravity_corrected_compliance

    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    segs = [
        read(str(p))
        for p in sorted((WORK / station / "clean").glob(f"{station}_seg*.mseed"))
    ]
    for t0 in pick:
        t0 = UTCDateTime(t0)
        for st in segs:
            if st[0].stats.starttime <= t0 and t0 + 2100 <= st[0].stats.endtime:
                w = (
                    (st.select(channel="BHZ") + st.select(channel="BDH"))
                    .slice(t0, t0 + 2100 - 0.8)
                    .copy()
                )
                for tr in w:
                    tr.data = tr.data.astype(float)
                f, zz, pp, zp, coh = spectra(w, NSEG)
                pp, zp = pp * gain**2, zp * gain
                axes[0, 0].semilogx(
                    f, 10 * np.log10(zz * (2 * np.pi * f) ** 4), color="0.5", lw=0.3
                )
                axes[0, 1].semilogx(f, 10 * np.log10(pp), color="0.5", lw=0.3)
                axes[1, 0].semilogx(f, coh, color="0.5", lw=0.3)
                ok = pp > 0
                eta = gravity_corrected_compliance(
                    f[ok], depth, np.abs(zp[ok]) / pp[ok], np.sqrt(pp[ok])
                )
                axes[1, 1].semilogx(f[ok], eta, color="0.5", lw=0.3)
                break
    inb = (d["freqs"] >= band[0]) & (d["freqs"] <= band[1])
    axes[1, 1].errorbar(d["freqs"][inb], d["compliance"][inb], yerr=d["window_std"][inb], fmt="o",
                        color="red", capsize=4, label="Stacked (error bars: window scatter)")  # fmt: skip
    axes[1, 0].axhline(summ["results"]["compy"]["threshold"], color="black", ls="--")
    for ax in axes.ravel():
        for x in band:
            ax.axvline(x, color="black", ls="--")
        ax.set_xlim(0.002, 0.1)
        ax.set_xlabel("Frequency [Hz]")
    axes[0, 0].set(
        ylabel="Acceleration PSD [dB]",
        title="(a) Corrected vertical acceleration",
        ylim=(-200, -120),
    )
    axes[0, 1].set(
        ylabel="Pressure PSD [dB]", title="(b) Calibrated pressure", ylim=(-10, 70)
    )
    axes[1, 0].set(ylabel="Coherence", title="(c) Coherence", ylim=(0, 1))
    axes[1, 1].set(
        ylabel="Normalized compliance [1/Pa]", title="(d) Compliance", ylim=(0, 1e-10)
    )
    axes[1, 1].legend(fontsize=10)
    fig.suptitle(f"{station}: {len(pick)} of {len(starts)} selected 35-min windows")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig07_windows.png", dpi=130)
    plt.close(fig)


# --- Figs. 8-10: compliance and inversion ----------------------------------------------------------


def heat(post):
    vs = post.vs_draws(Z, max_draws=8000)
    h = np.array([np.histogram(v, VS_BINS)[0] for v in vs.T], float)
    h /= np.maximum(h.max(axis=1, keepdims=True), 1)
    return h, np.median(vs, axis=0)


def heat_axis(ax, post, column, cbar=True, legend=False):
    h, med = heat(post)
    im = ax.imshow(h, extent=(0, 5000, -Z[-1], 0), aspect="auto", cmap=CMAP, vmin=0, vmax=1,
                   origin="upper", interpolation="nearest")  # fmt: skip
    cv, cz = crust1_steps(column)
    ax.plot(cv, -cz, color="white", ls="--", lw=1.8, label="CRUST1.0")
    ax.plot(med, -Z, color="red", lw=2, label="Median")
    ax.set_xlim(0, 5000)
    ax.set_ylim(-8000, 0)
    ax.grid(color="white", alpha=0.25)
    if cbar:
        plt.colorbar(im, ax=ax, fraction=0.06, pad=0.02, label="Probability")
    if legend:
        ax.legend(
            loc="lower left",
            facecolor="0.15",
            labelcolor="white",
            framealpha=0.6,
            fontsize=9,
        )


def fig08_compliance(store, c):
    groups = {}
    for s in store:
        groups.setdefault(sub_array(s), []).append(s)
    order = [g for g in dict.fromkeys(GROUPS.values()) if g in groups] + sorted(
        set(groups) - set(GROUPS.values())
    )
    fig, axes = plt.subplots(1, len(order), figsize=(5 * len(order), 6), sharey=True)
    for ax, g in zip(np.atleast_1d(axes), order):
        for s in groups[g]:
            d = store[s][0].data
            line = ax.plot(d.freqs, d.eta * 1e9, lw=2, label=s)[0]
            ax.fill_between(d.freqs, (d.eta - d.sigma) * 1e9, (d.eta + d.sigma) * 1e9,
                            color=line.get_color(), alpha=0.25, lw=0)  # fmt: skip
        ax.set_yscale("log")
        ax.set_xlim(0.004, 0.017)
        ax.set_title(g)
        ax.set_xlabel("Frequency [Hz]")
        ax.legend(loc="lower right", fontsize=9)
    np.atleast_1d(axes)[0].set_ylabel("Normalized compliance [1/GPa]")
    fig.suptitle(
        "Normalized compliance (±1σ of the full covariance, gain error excluded: gain is sampled)"
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig08_compliance.png", dpi=130)
    plt.close(fig)


def fig09_heatmaps(store, c):
    prefixes = list(GROUPS)
    names = sorted(
        store,
        key=lambda s: (
            prefixes.index(s[:2]) if s[:2] in prefixes else len(prefixes),
            s,
        ),
    )
    ncol = 5
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(
        nrow, ncol, figsize=(4.6 * ncol, 4.6 * nrow), sharey=True, squeeze=False
    )
    for k, ax in enumerate(axes.ravel()):
        if k >= len(names):
            ax.axis("off")
            continue
        s = names[k]
        heat_axis(ax, *store[s], cbar=False, legend=k == 0)
        ax.set_title(title(s, c), fontsize=11)
        if k % ncol == 0:
            ax.set_ylabel("Depth [m]")
        if k >= len(names) - ncol:
            ax.set_xlabel("Shear velocity [m/s]")
    fig.suptitle(
        "Vs posterior (free crustal Vp/Vs, pressure gain sampled within ±5 %); rows normalized"
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig09_heatmaps.png", dpi=120)
    plt.close(fig)


def fig10_profiles(store, c):
    groups = {}
    for s in store:
        groups.setdefault(sub_array(s), []).append(s)
    fig, axes = plt.subplots(1, 3, figsize=(20, 9), sharey=True)
    colors = plt.get_cmap("tab10")
    for k, (g, stations) in enumerate(sorted(groups.items())):
        meds = np.array(
            [np.median(store[s][0].vs_draws(Z, 4000), axis=0) for s in stations]
        )
        for m in meds:
            axes[0].plot(m, -Z, color=colors(k), lw=0.8, alpha=0.6)
        axes[0].plot(
            np.median(meds, axis=0),
            -Z,
            color=colors(k),
            lw=2.8,
            label=f"{g} (n={len(stations)})",
        )
    avg = np.mean([crust1_vs_at(store[s][1], Z) for s in store], axis=0)
    axes[0].plot(avg, -Z, "k--", lw=2, label="CRUST1.0 (mean)")
    axes[0].set(xlim=(0, 5000), xlabel="Shear velocity [m/s]", ylabel="Depth [m]",
                title="(a) Median Vs by sub-array (thin: stations)")  # fmt: skip
    axes[0].legend(loc="lower left", fontsize=10)
    # (b) network: pooled posterior
    pooled = np.concatenate([store[s][0].vs_draws(Z, 1500) for s in store])
    q = np.quantile(pooled, [0.05, 0.25, 0.5, 0.75, 0.95], axis=0)
    axes[1].fill_betweenx(-Z, q[0], q[4], color="tab:blue", alpha=0.2, label="5–95 %")
    axes[1].fill_betweenx(-Z, q[1], q[3], color="tab:blue", alpha=0.35, label="25–75 %")
    axes[1].plot(q[2], -Z, color="tab:blue", lw=2.5, label="Median, all stations")
    mains = [load_posterior(s, "main")[0] for s in store]
    mains = [m for m in mains if m is not None]
    if mains:
        pooled_main = np.concatenate([m.vs_draws(Z, 1500) for m in mains])
        axes[1].plot(np.median(pooled_main, axis=0), -Z, color="tab:orange", lw=2.2,
                     label=f"Median, Vp/Vs = 1.85 ({len(mains)} stations)")  # fmt: skip
    axes[1].plot(avg, -Z, "k--", lw=2, label="CRUST1.0 (mean)")
    axes[1].axvline(4200, color="grey", ls=":", lw=2)
    axes[1].text(4250, -6000, "fresh peridotite", rotation=90, color="grey")
    axes[1].set(
        xlim=(0, 5000),
        xlabel="Shear velocity [m/s]",
        title="(b) Pooled posterior (free Vp/Vs) and Vp/Vs = 1.85",
    )
    axes[1].legend(loc="lower left", fontsize=10)
    beta = ci.serpentinization(pooled)
    qb = np.quantile(beta, [0.05, 0.5, 0.95], axis=0)
    axes[2].fill_betweenx(
        -Z, qb[0], qb[2], color="tab:red", alpha=0.25, label="5–95 % (pooled draws)"
    )
    axes[2].plot(qb[1], -Z, color="tab:red", lw=2.5, label="Median")
    axes[2].set(xlim=(0, 100), xlabel="Serpentinization [%]",
                title="(c) Upper bound if low Vs were only serpentinite")  # fmt: skip
    axes[2].legend(loc="lower right", fontsize=10)
    axes[0].set_ylim(-8000, 0)
    fig.suptitle("Serpentinization from Vs: 4.2 km/s = 0 %, 2.2 km/s = 100 % (Miller & Christensen 1997; "
                 "Aminian et al. 2025, Fig. 10c)")  # fmt: skip
    fig.tight_layout()
    fig.savefig(FIGURES / "fig10_profiles.png", dpi=130)
    plt.close(fig)


def misfit_axis(ax, post, title_=None):
    n = len(post.data.eta)
    it = np.arange(1, post.trace_chi2.shape[1] + 1)
    keep = np.unique(np.geomspace(1, len(it), 2000).astype(int)) - 1
    for trace in post.trace_chi2:
        ax.plot(it[keep], trace[keep] / n, lw=0.6, alpha=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axvline(30000, color="red", ls="--", lw=1.2)
    ax.axhline(1.0, color="red", ls=":", lw=1.2)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("χ² / N")
    if title_:
        ax.set_title(title_, fontsize=10)


def figS_misfit(store, c):
    names = sorted(store)
    ncol = 6
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(
        nrow, ncol, figsize=(3.6 * ncol, 3.0 * nrow), squeeze=False, sharex=True
    )
    for k, ax in enumerate(axes.ravel()):
        if k < len(names):
            misfit_axis(ax, store[names[k]][0], names[k])
        else:
            ax.axis("off")
    fig.suptitle("Misfit of every chain (red dashed: end of burn-in; dotted: χ²/N = 1)")
    fig.tight_layout()
    fig.savefig(FIGURES / "figS_misfit.png", dpi=110)
    plt.close(fig)


def station_figure(s, post, column, c):
    fig, axes = plt.subplots(
        1, 3, figsize=(19, 6), gridspec_kw=dict(width_ratios=[1, 1, 0.9])
    )
    d = post.data
    pred = post.pred.reshape(-1, len(d.eta))
    for i in np.random.default_rng(0).choice(len(pred), 300, replace=False):
        axes[0].plot(d.freqs, pred[i] * 1e9, color="green", lw=0.25, alpha=0.4)
    axes[0].errorbar(
        d.freqs,
        d.eta * 1e9,
        yerr=d.sigma * 1e9,
        fmt="o",
        color="black",
        capsize=4,
        label="Measured",
    )
    axes[0].plot([], [], color="green", label="Posterior models")
    axes[0].set(
        xlabel="Frequency [Hz]",
        ylabel="Normalized compliance [1/GPa]",
        title="(a) Data fit",
    )
    axes[0].legend()
    misfit_axis(axes[1], post, "(b) Misfit function")
    heat_axis(axes[2], post, column, legend=True)
    axes[2].set(
        title="(c) " + title(s, c), xlabel="Shear velocity [m/s]", ylabel="Depth [m]"
    )
    fig.tight_layout()
    fig.savefig(FIGURES / f"{s}_inversion.png", dpi=110)
    plt.close(fig)


# --- tables ------------------------------------------------------------------------------------------------


def q3(x):
    a = np.quantile(x, [0.05, 0.5, 0.95])
    return a


def tables(c, store):
    TABLES.mkdir(parents=True, exist_ok=True)
    rows = []
    for s in STATIONS:
        summ = compliance_summary(s)
        cal_path = WORK / s / "calibration" / "compy_gain.json"
        cal = json.loads(cal_path.read_text()) if cal_path.exists() else {}
        lat, lon, depth = c[s]
        row = dict(station=s, lat=lat, lon=lon, depth_m=depth, gain_compy=cal.get("gain"),
                   gain_compy_std=cal.get("gain_std"), gain_events=len(cal.get("events_used", [])),
                   gain_bruitfm=(summ or {}).get("gains", {}).get("bruitfm"))  # fmt: skip
        label = None
        for lab in ("compy", "bruitfm", "nominal"):
            if summ and "n_windows" in summ["results"].get(lab, {}):
                label = lab
                break
        if label:
            r = summ["results"][label]
            row.update(gain_used=label, band_lo=r["band_hz"][0], band_hi=r["band_hz"][1], threshold=r["threshold"],
                       n_windows=r["n_windows"], coherence=r["median_coherence_in_band"],
                       months=len(r["monthly"]))  # fmt: skip
        inv = RESULTS / s / "inversion" / f"{s}_final.json"
        if inv.exists():
            j = json.loads(inv.read_text())
            row.update(
                chi2=j["chi2"]["median"],
                n_data=j["n_data"],
                rhat=j["max_rhat"],
                ess=j["min_ess"],
            )
            for v in ("main", "final_h2"):
                p = RESULTS / s / "inversion" / f"{s}_{v}.json"
                if p.exists():
                    row[f"chi2_{v}"] = json.loads(p.read_text())["chi2"]["median"]
        if s in store:
            post = store[s][0]
            flat = post.flat()
            names = post.names
            row["vp_vs"] = q3(flat[:, names.index("vp_vs")]).tolist()
            row["sediment_m"] = q3(flat[:, names.index("sediment_thickness")]).tolist()
            if "gain" in names:
                row["gain_sampled"] = (
                    q3(flat[:, names.index("gain")]) * (row.get("gain_compy") or 1)
                ).tolist()
            vs = post.vs_draws(np.array([250.0, 500, 1000, 2000, 3000, 5000]), 6000)
            for k, zz in enumerate((250, 500, 1000, 2000, 3000, 5000)):
                row[f"vs_{zz}"] = q3(vs[:, k]).tolist()
            beta = ci.serpentinization(
                post.vs_draws(np.array([500.0, 1000, 2000, 3000]), 6000)
            )
            for k, zz in enumerate((500, 1000, 2000, 3000)):
                row[f"serp_{zz}"] = q3(beta[:, k]).tolist()
            row["crust1_sediment_m"] = store[s][1].get("sediment_thickness")
            row["crust1_moho_m"] = store[s][1].get("moho_depth")
        rows.append(row)
    (TABLES / "stations.json").write_text(json.dumps(rows, indent=1, default=float))
    # Table 1 (paper style)
    lines = ["| Station | Latitude [°] | Longitude [°] | Depth [m] | Gain factor (ComPy, Rayleigh) | Events | "
             "Gain (BRUIT-FM microseism) | Frequency* [Hz] | Threshold | Windows |",
             "|---|---|---|---|---|---|---|---|---|---|"]  # fmt: skip
    for r in rows:
        g = f"{r['gain_compy']:.3f} ± {r['gain_compy_std']:.3f}" if r.get("gain_compy_std") else (
            f"{r['gain_compy']:.3f}" if r.get("gain_compy") else "–")  # fmt: skip
        band = f"{r['band_lo']:.4f}–{r['band_hi']:.4f}" if "band_lo" in r else "–"
        gb = f"{r['gain_bruitfm']:.3f}" if r.get("gain_bruitfm") else "–"
        lines.append(f"| {r['station']} | {r['lat']:.4f} | {r['lon']:.4f} | {r['depth_m']:.0f} | {g} | "
                     f"{r['gain_events'] or '–'} | {gb} | {band} | {r.get('threshold', '–')} | "
                     f"{r.get('n_windows', '–')} |")  # fmt: skip
    (TABLES / "table1.md").write_text("\n".join(lines) + "\n")
    with open(TABLES / "table1.csv", "w") as fh:
        keys = sorted({k for r in rows for k in r})
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    return rows


def main():
    FIGURES.mkdir(exist_ok=True)
    c = coords()
    store = {}
    for s in STATIONS:
        if usable(s):
            store[s] = load_posterior(s)
    status = {s: ("inverted" if s in store else "compliance only" if compliance_summary(s) else "not processed")
              for s in c}  # fmt: skip
    rows = tables(c, store)
    fig01_map(c, status)
    ref = sys.argv[1] if len(sys.argv) > 1 else "CC11"
    summ = compliance_summary(ref)
    fig03_psd(ref)
    fig04_tilt(ref)
    fig05_calibration(ref)
    fig06_coherogram(ref, summ)
    fig07_windows(ref, summ)
    if store:
        fig08_compliance(store, c)
        fig09_heatmaps(store, c)
        fig10_profiles(store, c)
        figS_misfit(store, c)
        for s in store:
            station_figure(s, *store[s], c)
    print(f"figures done: {len(store)} inverted stations, {len(rows)} table rows")


if __name__ == "__main__":
    main()
