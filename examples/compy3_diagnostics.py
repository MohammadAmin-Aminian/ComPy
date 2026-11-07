"""ComPy 3 audit: diagnostics, anomaly ranking and checks on archived results.

Runs on the stage C/E products committed in ``docs/results`` (no waveforms,
no ObsPy) and writes figures, the ranked anomaly table and the numbers
quoted in ``docs/audit_2026-11/REPORT.md``::

    MPLBACKEND=Agg python examples/compy3_diagnostics.py --output docs/audit_2026-11

Sections (each writes into ``<output>/figures`` and ``<output>/numbers.json``):

1. validation: recompute every published compliance from its stored stack;
2. anomaly scan of the YV and XE stations, ranked table;
3. residual Z/P phase and instrument-response models;
4. epoch mixing (steps inside the stacked period);
5. coherence-gate selection bias (Monte Carlo);
6. water-column assumptions of the Rayleigh-wave calibration;
7. time-frequency figures on a synthetic record (the waveforms are not in
   the repository; ``examples/compy3_waveform_diagnostics.py`` makes the
   same figures from cleaned MiniSEED).
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from compy import anomaly as an  # noqa: E402
from compy import results  # noqa: E402
from compy.compliance import _compliance  # noqa: E402
from compy.plotting import diagnostics as dg  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
YV = ROOT / "docs/results/yv_compliance_2026-10-28/results"
YV_INV = ROOT / "docs/results/yv_inversion_2026-10-29/results"
XE = ROOT / "docs/results/xe_young_orca_2026-10-29/results"
DPI = 110


def save(fig, path):
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    print("wrote", path.relative_to(ROOT) if path.is_relative_to(ROOT) else path)


def iso(t):
    return datetime.fromtimestamp(float(t), tz=timezone.utc).strftime("%Y-%m-%d")


# --------------------------------------------------------------------------
def validation():
    """Recompute compliance and coherence of every archived result from its stack."""
    rows = []
    for path in sorted(ROOT.glob("docs/results/*/results/**/*_compliance_*.npz")):
        with np.load(path) as d:
            f, spp, spz, szz = d["all_freqs"], d["S_pp"], d["S_pz"], d["S_zz"]
            ratio = np.zeros_like(spp)
            np.divide(np.abs(spz), spp, out=ratio, where=spp > 0)
            eta = _compliance(f, float(d["depth"]), ratio, spp, True)
            keep = (f >= d["freqs"][0] - 1e-12) & (f <= d["freqs"][-1] + 1e-12)
            coh = np.abs(spz[keep]) / np.sqrt(szz[keep] * spp[keep])
            rows.append(dict(file=str(path.relative_to(ROOT)),
                             max_rel_compliance=float(np.max(np.abs(eta[keep] / d["compliance"] - 1))),
                             max_abs_coherence=float(np.max(np.abs(coh - d["coherence"])))))  # fmt: skip
    worst = max(r["max_rel_compliance"] for r in rows)
    print(f"validation: {len(rows)} files, worst relative difference {worst:.2e}")
    return dict(files=len(rows), worst_relative=worst,
                worst_coherence=max(r["max_abs_coherence"] for r in rows))  # fmt: skip


# --------------------------------------------------------------------------
def scan(figures):
    yv = results.discover(YV, label="compy", inversion_root=YV_INV)
    xe = results.discover(XE)
    found, infos, nets = [], {}, {}
    for net, group in (("YV", yv), ("XE", xe)):
        infos[net] = {}
        for code, s in group.items():
            a, info = an.scan_station(s)
            found += a
            infos[net][code] = info
        a, ninfo = an.scan_network(group, infos[net], name=net)
        found += a
        nets[net] = ninfo
    ranked = an.rank(found, 2.0)
    return yv, xe, ranked, infos, nets


def figures_from_scan(yv, xe, ranked, infos, nets, figures):
    save(dg.plot_anomaly_ranking(ranked, top=30), figures / "anomaly_ranking.png")
    models = {c: i["phase_models"] for net in infos.values() for c, i in net.items()
              if i.get("phase_models") is not None}  # fmt: skip
    groups = {
        c: ("YV (INSU-IPGP)" if c.startswith("RR") else "XE (Scripps)") for c in models
    }
    save(
        dg.plot_phase_diagnostic(models, groups=groups),
        figures / "phase_diagnostic.png",
    )
    for code in ("RR52", "RR50", "RR28"):
        s = yv[code]
        lv = infos["YV"][code]["levels"]
        save(dg.plot_window_levels(s, lv), figures / f"{code}_window_levels.png")
    for code in ("RR52", "RR50"):
        steps = [
            c.time
            for c in infos["YV"][code]["levels"]["changepoints"]["transfer function"]
        ]
        save(
            dg.plot_compliance_stability(yv[code], steps=steps),
            figures / f"{code}_monthly_stability.png",
        )
        save(
            dg.plot_compliance_confidence(yv[code]),
            figures / f"{code}_error_budget.png",
        )
    save(
        dg.plot_compliance_stability(xe["CC07"]), figures / "CC07_monthly_stability.png"
    )
    ray = {c: i["acceptance"]["rayleigh"] for net in infos.values() for c, i in net.items()
           if "rayleigh" in i["acceptance"]}  # fmt: skip
    save(dg.plot_periodicity(ray, title="Accepted-window clustering by tidal/solar period "
                             "(● empirical and tail p ≤ 0.01)"), figures / "tidal_acceptance_map.png")  # fmt: skip
    save(dg.plot_acceptance_phase(yv), figures / "YV_acceptance_M2.png")
    save(dg.plot_acceptance_phase(yv, period_h=24.0, title="YV accepted windows folded on 24 h (UTC)"),
         figures / "YV_acceptance_S1.png")  # fmt: skip
    per = {
        c: i["levels"]["periodicity"]
        for c, i in infos["YV"].items()
        if "periodicity" in i["levels"]
    }
    save(dg.plot_periodicity(per, value="amplitude",
                             title="|Z/P| modulation by tidal/solar period, YV (● p ≤ 0.01)"),
         figures / "YV_tidal_transfer_function.png")  # fmt: skip
    for net, group in (("YV", yv), ("XE", xe)):
        save(dg.plot_multistation(group, title=f"{net}: compliance (±1σ total) and Z–P coherence"),
             figures / f"{net}_multistation.png")  # fmt: skip
        save(dg.plot_common_mode(nets[net]["common_mode"],
                                 title=f"{net}: monthly deviations (left) and their correlation (right)"),
             figures / f"{net}_common_mode.png")  # fmt: skip
    save(dg.plot_noise_models({**yv, **xe}, kind="vertical", highlight=("RR52", "RR50", "CC11")),
         figures / "noise_vertical.png")  # fmt: skip
    save(dg.plot_noise_models({**yv, **xe}, kind="pressure", highlight=("RR52", "RR50", "CC11")),
         figures / "noise_pressure.png")  # fmt: skip
    return models


# --------------------------------------------------------------------------
def epoch_mixing(yv, infos):
    """Effect of steps inside the stacked periods actually used (full or epoch)."""
    out = {}
    for code, s in yv.items():
        t, r, p_db, _ = np.asarray(s.window_levels, float).T
        cps = infos["YV"][code]["levels"]["changepoints"]["transfer function"]
        epoch = YV / code / "epoch" / f"{code}_summary.json"
        end = json.loads(epoch.read_text())["epoch_end"] if epoch.exists() else None
        end_t = (
            None
            if end is None
            else datetime.fromisoformat(end).replace(tzinfo=timezone.utc).timestamp()
        )
        out[code] = dict(
            steps=[(iso(c.time), round(an._pct_from_log10(c.shift), 2), round(c.significance, 1)) for c in cps],
            published_epoch_end=end,
            full_period_shift_pct=round(an.stack_shift(t, r, p_db, cps), 2),
            used_period_shift_pct=round(an.stack_shift(t, r, p_db, cps, end_t), 2),
            windows_after_first_step_in_used_period=int(np.sum((t >= min([c.time for c in cps], default=np.inf))
                                                               & (t < (end_t or np.inf)))),
        )  # fmt: skip
    return out


def selection_bias():
    """Monte Carlo table of the coherence-gate bias (inclusive vs exclusive)."""
    rows = []
    for coherence in (0.8, 0.85, 0.9, 0.95):
        for sigma in (0.0, 0.5, 1.0):
            inc = [
                an.simulate_gate_bias(coherence, seastate_sigma=sigma, seed=k)
                for k in range(3)
            ]
            exc = [an.simulate_gate_bias(coherence, seastate_sigma=sigma, seed=k, exclusive=True)
                   for k in range(3)]  # fmt: skip
            rows.append(dict(coherence=coherence, seastate_sigma=sigma,
                             acceptance=round(np.mean([x["acceptance"] for x in inc]), 3),
                             inclusive_bias_pct=round(np.mean([x["bias_pct"] for x in inc]), 2),
                             inclusive_sd=round(np.std([x["bias_pct"] for x in inc]), 2),
                             exclusive_bias_pct=round(np.mean([x["bias_pct"] for x in exc]), 2),
                             exclusive_sd=round(np.std([x["bias_pct"] for x in exc]), 2)))  # fmt: skip
    return rows


def station_selection_bias(stations):
    """Bias expected at each station's acceptance, for two sea-state spreads."""
    out = {}
    for code, s in stations.items():
        acc = s.acceptance
        thr = s.summary["results"][s.label].get("threshold", 0.9)
        if not acc:
            continue
        sigma = None
        if (
            s.window_levels is not None
        ):  # spread of accepted pressure levels [ln amplitude]
            sigma = float(np.std(s.window_levels[:, 2]) / (20 / np.log(10)))
        row = dict(acceptance=round(acc, 3), gate=thr)
        for label, sg in (("sigma0.5", 0.5), ("sigma_obs", sigma)):
            if sg is None:
                continue
            e = an.expected_gate_bias(acc, gate=thr, seastate_sigma=sg)
            row[label] = dict(sigma=round(sg, 2), coherence=round(e["coherence"], 3),
                              bias_pct=round(e["bias_pct"], 2))  # fmt: skip
        out[code] = row
    return out


# --------------------------------------------------------------------------
def water_ratio(f, c, depth, density, sound_speed):
    """Zha & Webb (2016, A5) P/a of a Rayleigh wave under a water layer (as compy_pressure)."""
    q2 = sound_speed**-2 - c**-2
    x = 2 * np.pi * f * depth * np.sqrt(np.abs(q2))
    return density * depth * np.where(q2 >= 0, np.tan(x) / x, np.tanh(x) / x)


def water_column(stations):
    """Gain bias from a homogeneous 1028 kg/m³, 1500 m/s water column."""
    disp = np.loadtxt(ROOT / "compy_pressure/reference_data/rr52_historical_rayleigh_dispersion.csv",
                      delimiter=",", skiprows=1)  # fmt: skip
    f = np.linspace(0.02, 0.06, 41)
    c = np.interp(f, disp[:, 0], disp[:, 1])
    beta, rho_deep = (
        4.4e-10,
        1027.8,
    )  # seawater compressibility [1/Pa], deep potential density
    out = {}
    for code, s in stations.items():
        h = s.depth
        rho_mean = rho_deep * (1 + 0.5 * beta * rho_deep * 9.81 * h)
        base = water_ratio(f, c, h, 1028.0, 1500.0)
        row = dict(depth_m=h, column_mean_density=round(rho_mean, 1))
        for cs in (1490.0, 1500.0, 1510.0):
            # g = expected/observed: proportional to the water-layer prediction
            ratio = water_ratio(f, c, h, rho_mean, cs) / base
            row[f"gain_change_pct_c{int(cs)}"] = round(100 * (np.median(ratio) - 1), 2)
        out[code] = row
    return out


# --------------------------------------------------------------------------
def synthetic_record(days=40, fs=0.2, seed=3):
    """Synthetic Z, P, H1, H2 with storms, M2 tidal-current tilt and a P step.

    Vertical = compliance response to infragravity pressure + tilt leakage
    of horizontal noise (azimuth 70°), whose amplitude follows a residual
    current plus an M2 tidal current; the pressure sensitivity drops 7 % on
    day 30. Units are arbitrary; only the diagnostics matter.
    """
    from scipy import signal

    rng = np.random.default_rng(seed)
    n = int(days * 86400 * fs)
    t = np.arange(n) / fs
    b, a = signal.butter(4, [0.004, 0.03], btype="band", fs=fs)
    storm = (
        1
        + 3 * np.exp(-0.5 * ((t / 86400 - 12) / 2.0) ** 2)
        + 2 * np.exp(-0.5 * ((t / 86400 - 26) / 1.5) ** 2)
    )
    p = signal.lfilter(b, a, rng.standard_normal(n)) * storm
    current = np.abs(0.6 + np.cos(2 * np.pi * t / (12.4206 * 3600)))
    bh, ah = signal.butter(2, 0.03, fs=fs)
    tilt = signal.lfilter(bh, ah, rng.standard_normal(n)) * current
    az = np.radians(70.0)
    other = signal.lfilter(bh, ah, rng.standard_normal((2, n)), axis=1) * current
    h1 = np.cos(az) * tilt - np.sin(az) * other[0] + 0.05 * rng.standard_normal(n)
    h2 = np.sin(az) * tilt + np.cos(az) * other[0] + 0.05 * rng.standard_normal(n)
    z = 0.8 * p + 0.6 * tilt + 0.05 * rng.standard_normal(n)
    p_rec = p * np.where(t > 30 * 86400, 0.93, 1.0)
    t0 = datetime(2013, 1, 1, tzinfo=timezone.utc).timestamp()
    return dict(Z=z, P=p_rec, **{"1": h1, "2": h2}), fs, t0


def synthetic_figures(figures):
    from compy import timefreq as tf

    channels, fs, t0 = synthetic_record()
    spec = tf.window_spectra(channels, fs, 7200, 256, t0=t0)
    band = (0.006, 0.02)
    save(dg.plot_spectrogram(spec["times"], spec["f"], {k: spec[k + k] for k in ("Z", "P", "1", "2")},
                             band=band, title="Synthetic record: spectrograms (storms days 12 and 26)"),
         figures / "synthetic_spectrogram.png")  # fmt: skip
    coh = {"Z–P": tf.coherence(spec, "P", "Z"), "Z–H1": tf.coherence(spec, "1", "Z"),
           "Z–H2": tf.coherence(spec, "2", "Z"), "P–H1": tf.coherence(spec, "P", "1")}  # fmt: skip
    save(dg.plot_coherogram(spec["times"], spec["f"], coh, band=band, threshold=0.9,
                            title="Synthetic record: coherograms (M2 tilt noise; P step on day 30)"),
         figures / "synthetic_coherogram.png")  # fmt: skip
    tilt = tf.tilt_direction(spec, (0.005, 0.035))
    save(dg.plot_tilt(spec["times"], tilt["azimuth"], tilt["admittance"], strength=tilt["coherence"],
                      angle_label="Admittance |S_HZ|/S_HH", strength_label="Z–H coherence",
                      title="Synthetic record: tilt azimuth (true 70°) and admittance (true 0.6)"),
         figures / "synthetic_tilt.png")  # fmt: skip
    adm = tf.band_average(spec["f"], tf.admittance(spec, "P", "Z"), band)
    cps = an.changepoints(spec["times"], np.log10(adm), min_size=20)
    return dict(true_azimuth=70.0, median_azimuth=float(np.median(tilt["azimuth"])),
                step_found=[(iso(c.time), round(an._pct_from_log10(c.shift), 1)) for c in cps])  # fmt: skip


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "docs/audit_2026-11")
    parser.add_argument(
        "--skip-mc", action="store_true", help="skip the Monte Carlo sections"
    )
    args = parser.parse_args()
    figures = args.output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    numbers = {"validation": validation()}
    yv, xe, ranked, infos, nets = scan(figures)
    models = figures_from_scan(yv, xe, ranked, infos, nets, figures)
    (args.output / "anomalies.md").write_text(an.to_markdown(ranked) + "\n")
    an.to_csv(ranked, args.output / "anomalies.csv")
    numbers["phase_models"] = {
        c: dict(residual_deg=np.round(m["residual_deg"], 2).tolist(), freqs_mhz=np.round(1e3 * m["freqs"], 2).tolist(),
                rms={k: round(m[k]["rms_deg"], 2) for k in ("constant", "timing", "dpg_corner", "seismometer_corner")},
                dt_s=round(m["timing"]["params"]["dt_s"], 2),
                dpg_bias_band_bottom_pct=[round(m["dpg_corner"]["bias_low_pct"][0], 1),
                                          round(m["dpg_corner"]["bias_high_pct"][0], 1)],
                seis_bias_band_bottom_pct=[round(m["seismometer_corner"]["bias_low_pct"][0], 1),
                                           round(m["seismometer_corner"]["bias_high_pct"][0], 1)],
                rayleigh_slope_dpg=[round(v, 3) for v in m["dpg_corner"]["rayleigh_slope"]],
                rayleigh_slope_seis=[round(v, 3) for v in m["seismometer_corner"]["rayleigh_slope"]])
        for c, m in models.items()
    }  # fmt: skip
    numbers["network"] = {
        net: dict(common_mode_mean_corr=round(n["common_mode"]["network_mean_corr"], 3),
                  common_mode_p=n["common_mode"]["p_value"], phase_sign=n.get("phase_sign"),
                  phase_vs_sediment=n.get("phase_vs_sediment"),
                  synchronous_steps=[[(iso(t), c, round(an._pct_from_log10(x), 1)) for t, c, x in g]
                                     for g in n.get("synchronous_steps", []) if len(g) >= 2])
        for net, n in nets.items()
    }  # fmt: skip
    numbers["epoch_mixing"] = epoch_mixing(yv, infos)
    numbers["water_column"] = water_column({**yv, **xe})
    if not args.skip_mc:
        numbers["selection_bias_mc"] = selection_bias()
        numbers["station_selection_bias"] = station_selection_bias({**yv, **xe})
    numbers["synthetic"] = synthetic_figures(figures)
    numbers["n_anomalies"] = len(ranked)
    (args.output / "numbers.json").write_text(
        json.dumps(numbers, indent=1, default=str) + "\n"
    )
    print("wrote", args.output / "numbers.json")


if __name__ == "__main__":
    main()
