"""Socrates scientific validation checks for docs/reviews/socrates_initial_audit.md.

Read-only. Uses synthetic signals, analytic/ODE solutions and files committed on
review/sentinel; never reads raw waveforms, masters or saved results, and writes
nothing. Run from the repository root with the package importable, e.g.

    PYTHONPATH=.. python docs/reviews/socrates_validation.py

Each section prints the numbers quoted in the audit (V1-V11).
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from obspy.core.inventory.response import Response
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.signal.windows import hann

from compy_pressure.rayleigh_calibration import water_ratio

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
DISPERSION = np.loadtxt(
    ROOT / "../compy_pressure/reference_data/rr52_historical_rayleigh_dispersion.csv",
    delimiter=",",
    skiprows=1,
)
INSU = {"RR28", "RR29", "RR31", "RR34", "RR36", "RR38", "RR40", "RR50", "RR52"}


def lagrangian_ratio(
    f, phase_velocity, depth, *, rho=1028.0, alpha=1500.0, g=9.81, stratified=True
):
    """Seafloor Lagrangian pressure / upward acceleration, acoustic-gravity ODE.

    z up from the floor; dp = p' - rho0*g*xi_z; dp = -rho0*alpha^2*div(xi).
    Free surface dp(H)=0; the gauge moves with the floor and records dp(0).
    """
    w = 2 * np.pi * f
    k = w / phase_velocity

    def rho0(z):
        return rho * np.exp(g * (depth - z) / alpha**2) if stratified else rho

    def rhs(z, y):
        xi, dp = y
        r = rho0(z)
        return [
            dp * (k**2 / (w**2 * r) - 1 / (r * alpha**2)) + g * k**2 / w**2 * xi,
            xi * r * (w**2 - g**2 * k**2 / w**2) - g * k**2 / w**2 * dp,
        ]

    sol = solve_ivp(rhs, (depth, 0.0), [1.0, 0.0], rtol=1e-11, atol=1e-14)
    xi0, dp0 = sol.y[:, -1]
    return dp0 / (-(w**2) * xi0)


def v1_acoustic_gravity():
    print("\nV1  water_ratio against the full acoustic-gravity ODE (c=3950 m/s)")
    for depth in (900, 2000, 2800, 4000, 5000):
        top = 0.07 if depth < 5000 else 0.06
        rows = []
        for f in np.arange(0.03, top + 1e-9, 0.01):
            base = float(water_ratio([f], [3950.0], depth, 1028, 1500)[0])
            no_g = lagrangian_ratio(f, 3950.0, depth, g=1e-12, stratified=False)
            grav = lagrangian_ratio(f, 3950.0, depth, stratified=False)
            full = lagrangian_ratio(f, 3950.0, depth)
            rows.append((no_g / base - 1, grav / base - 1, full / base - 1))
        a = 100 * np.array(rows)
        print(
            f"  H={depth} m: g->0 max|err| {abs(a[:, 0]).max():.1e}%; "
            f"gravity {a[:, 1].min():+.2f}..{a[:, 1].max():+.2f}%; "
            f"gravity+in-situ density {a[:, 2].min():+.2f}..{a[:, 2].max():+.2f}%"
        )


def v2_zha_offset_models():
    print("\nV2  Model variants: change of reported coefficient (1/g), YL peers")
    table = pd.read_csv(DOCS / "yl_all_zha_comparison.csv")
    peers = table[table.relative_difference_percent.notna() & (table.station != "C07W")]
    print(
        "  depth used == Zha Table 1 depth at every station:",
        bool((peers.depth_m == peers.published_depth_m).all()),
    )
    for band in ((0.03, 0.07), (0.04, 0.10)):
        f = np.arange(band[0], band[1] + 1e-9, 0.005)
        c = np.interp(f, DISPERSION[:, 0], DISPERSION[:, 1])
        change = {}
        for _, s in peers.iterrows():
            ref = water_ratio(f, c, s.depth_m, 1028, 1500)
            k = 2 * np.pi * f / c
            variants = {
                "incompressible": 1028 * np.tanh(k * s.depth_m) / k,
                "rho_1000": water_ratio(f, c, s.depth_m, 1000, 1500),
                "in_situ_density": water_ratio(
                    f, c, s.depth_m, 1028.8 + 2.24e-3 * s.depth_m, 1500
                ),
            }
            for name, model in variants.items():
                g_alt = np.dot(ref, model) / np.dot(ref, ref)
                change.setdefault(name, []).append(100 * (1 / g_alt - 1))
        summary = ", ".join(
            f"{k} median {np.median(v):+.1f}% (range {min(v):+.1f}..{max(v):+.1f})"
            for k, v in change.items()
        )
        print(f"  band {band}: {summary}")
    offsets = peers.relative_difference_percent
    print(
        f"  observed offset median {offsets.median():+.2f}%, range "
        f"{offsets.min():+.2f}..{offsets.max():+.2f}%; correlation with depth "
        f"{np.corrcoef(offsets, peers.depth_m)[0, 1]:+.2f}"
    )


def v3_b01w_highpass():
    print("\nV3  Single-pole high-pass, joint amplitude+phase fit (gamma^2 >= 0.9)")
    data = pd.read_csv(DOCS / "sentinel/b01w_boundary_events_spectra.csv")
    for station in ("B01W", "B04W", "A02W"):
        rows = data[(data.station == station) & (data.coherence_squared >= 0.9)]
        for event, e in rows.groupby("event_id"):
            f = e.frequency_hz.to_numpy()
            ratio = e.ratio_over_theory.to_numpy()
            phase = np.radians(e.phase_deg.to_numpy())

            def residual(p, f=f, ratio=ratio, phase=phase):
                gain, corner = p
                h = gain * (1j * f / corner) / (1 + 1j * f / corner)
                return np.r_[np.log(abs(h) / ratio), np.angle(h) - phase]

            fit = least_squares(residual, [0.3, 0.07], bounds=([1e-3, 1e-3], [10, 1]))
            gain, corner = fit.x
            n = len(f)
            print(
                f"  {station} {event[:12]:12s} bins {n:2d}: gain {gain:.3f}, "
                f"corner {corner:.3f} Hz (tau {1 / (2 * np.pi * corner):6.1f} s), "
                f"rms amplitude {100 * np.sqrt(np.mean(fit.fun[:n] ** 2)):.1f}%, "
                f"rms phase {np.degrees(np.sqrt(np.mean(fit.fun[n:] ** 2))):.1f} deg"
            )
            if station == "B01W":
                # Zha & Webb (2016) fit theory*C to the observed ratio, 0.04-0.10 Hz.
                fz = np.arange(0.04, 0.1001, 0.005)
                cz = np.interp(fz, DISPERSION[:, 0], DISPERSION[:, 1])
                tz = water_ratio(fz, cz, 1723, 1028, 1500)
                hz = gain / np.sqrt(1 + (corner / fz) ** 2)
                coefficient = np.dot(tz * hz, tz) / np.dot(tz, tz)
                print(
                    f"      model coefficient over 0.04-0.10 Hz {coefficient:.3f} "
                    "(Zha & Webb 2016 Table 1: 0.19 +/- 0.02)"
                )


def v4_c07w_states():
    print("\nV4  C07W two-state vertical amplitude ratio")
    daily = pd.read_csv(DOCS / "sentinel/c07w_daily_vertical_ratio.csv")
    rng = np.random.default_rng(0)
    for column in ("z_c07w_over_c08w", "z_c07w_over_c09w"):
        normal = daily.loc[daily.state == "normal", column].to_numpy()
        reduced = daily.loc[daily.state == "reduced", column].to_numpy()
        boots = [
            np.median(rng.choice(normal, normal.size))
            / np.median(rng.choice(reduced, reduced.size))
            for _ in range(5000)
        ]
        low, high = np.percentile(boots, [2.5, 97.5])
        ratio = np.median(normal) / np.median(reduced)
        print(
            f"  {column}: normal/reduced {ratio:.3f} "
            f"(95% IID day-bootstrap {low:.3f}-{high:.3f}; serial dependence omitted)"
        )
    gains = pd.read_csv(DOCS / "sentinel/c07w_event_gain_by_state.csv")
    normal = gains.loc[gains.vertical_state == "normal", "pressure_multiplier"]
    reduced = gains.loc[gains.vertical_state == "reduced", "pressure_multiplier"]
    print(
        "  event gains (0.03-0.07 Hz) normal/reduced medians: "
        f"{normal[normal > 0.4].median() / reduced.median():.3f}"
    )


def v5_taper_then_differentiate(trials=300):
    print("\nV5  Bias of i*w*FFT(w*v) vs true acceleration, 9 x 200 s Hann segments")
    rng = np.random.default_rng(1)
    segment, count = 200, 9
    n = segment * count
    full = np.fft.rfftfreq(4 * n, 1.0)
    window = hann(segment, sym=False)
    f = np.fft.rfftfreq(segment, 1.0)
    band = (f >= 0.03) & (f <= 0.07)
    for slope in (0, -2, -4, 2):
        ratios = []
        for _ in range(trials):
            spec = rng.standard_normal(full.size) + 1j * rng.standard_normal(full.size)
            shape = np.where(full > 0, np.maximum(full, 0.005) ** (slope / 2), 0)
            shape[full > 0.2] *= np.exp(-(full[full > 0.2] - 0.2) / 0.02)
            v = np.fft.irfft(spec * shape, 4 * n)[n : 2 * n]
            a = np.fft.irfft(2j * np.pi * full * spec * shape, 4 * n)[n : 2 * n]
            spp = saa = 0
            for i in range(count):
                cut = slice(i * segment, (i + 1) * segment)
                p_fft = np.fft.rfft((a[cut] - a[cut].mean()) * window)[band]
                v_fft = np.fft.rfft((v[cut] - v[cut].mean()) * window)[band]
                spp = spp + abs(p_fft) ** 2
                saa = saa + abs(2j * np.pi * f[band] * v_fft) ** 2
            ratios.append(np.sqrt(spp / saa))
        r = np.mean(ratios, axis=0)
        # Production fit with T=1: g = sum(R*T)/sum(R^2); a true gain is g=1.
        bias = 100 * (r.sum() / np.dot(r, r) - 1)
        print(
            f"  velocity power ~ f^{slope:+d}: R/T at 0.03 Hz {r[0]:.3f}, "
            f"at 0.07 Hz {r[-1]:.3f}; LSQ gain bias {bias:+.2f}%"
        )


def v6_estimator_bias(trials=20000):
    print("\nV6  Gain bias of R, |H1|, |H2| under production QC (K=9, gamma^2>=0.95)")
    rng = np.random.default_rng(2)

    def noise(*shape):
        return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / 2**0.5

    shape = (trials, 9, 9)
    for label, in_p, in_a in (("P noise", 1, 0), ("a noise", 0, 1), ("both", 1, 1)):
        for coherence in (0.95, 0.98):
            level = (
                np.sqrt(1 / np.sqrt(coherence) - 1)
                if in_p and in_a
                else np.sqrt(1 / coherence - 1)
            )
            s = noise(*shape)
            p = s + in_p * level * noise(*shape)
            a = s + in_a * level * noise(*shape)
            spp, saa = (abs(p) ** 2).mean(1), (abs(a) ** 2).mean(1)
            spa = (p * a.conj()).mean(1)
            g2 = abs(spa) ** 2 / (spp * saa)
            valid = g2 >= 0.95
            ok = (np.median(g2, 1) >= 0.95) & (g2.mean(1) >= 0.95) & (valid.sum(1) >= 3)
            out = []
            for name, x in (
                ("R", np.sqrt(spp / saa)),
                ("|H1|", abs(spa) / saa),
                ("|H2|", spp / abs(spa)),
            ):
                xs = np.where(valid, x, 0.0)[ok]
                gain = xs.sum(1) / (xs**2).sum(1)
                out.append(f"{name} {100 * (np.median(gain) - 1):+.2f}%")
            print(
                f"  {label:7s} true gamma^2={coherence:.2f}: accepted "
                f"{100 * ok.mean():5.1f}%; " + ", ".join(out)
            )


def v7_conventions():
    print("\nV7  evalresp conventions for a Pa-input high-pass (pole -p, zero 0)")
    p, f = 0.05, np.array([0.01, 0.03, 0.07])
    w = 2 * np.pi * f
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        kwargs = dict(
            zeros=[0j],
            stage_gain=1000.0,
            input_units="PA",
            output_units="COUNTS",
            normalization_frequency=1.0,
        )
        stable = Response.from_paz(poles=[-p + 0j], **kwargs)
        unstable = Response.from_paz(poles=[p + 0j], **kwargs)
        vel = stable.get_evalresp_response_for_frequencies(f, output="VEL")
        default = stable.get_evalresp_response_for_frequencies(f, output="DEF")
        flipped = unstable.get_evalresp_response_for_frequencies(f, output="DEF")
    print("  Pa response, output VEL == DEF:", bool(np.allclose(vel, default)))
    print(
        "  phase (deg) evalresp",
        np.round(np.angle(default, deg=True), 3),
        "vs numpy-convention iw/(iw+p)",
        np.round(np.degrees(np.arctan(p / w)), 3),
    )
    print(
        "  +p pole: amplitude ratio",
        np.round(abs(flipped) / abs(default), 6),
        "phase shift",
        np.round(np.angle(flipped / default, deg=True), 2),
        "= -2atan(p/w)",
        np.round(-np.degrees(2 * np.arctan(p / w)), 2),
    )


def v8_yv_pole_status():
    print("\nV8  YV station outcome against the metadata positive-real pressure pole")
    summary = json.loads((DOCS / "yv_rayleigh_completed_summary.json").read_text())
    review = json.loads((DOCS / "yv_rayleigh_independent_validation.json").read_text())
    poles = {
        r["station"]: (r["pressure_positive_real_poles"] or [[0.0]])[0][0]
        for r in review["response_flags"]
    }
    rows = []
    for s in summary["stations"]:
        station = s["seed_id"].split(".")[1]
        pole = poles.get(station, 0.0)
        error = np.degrees(2 * np.arctan(pole / (2 * np.pi * np.array([0.03, 0.07]))))
        rows.append(
            dict(
                station=station,
                pool="INSU DPG" if station in INSU else "DEPAS/Geomar hydrophone",
                pole_class="none"
                if pole == 0
                else "p<0.06"
                if pole < 0.06
                else "p>=0.06",
                status=s["status"].split("_")[0],
                gain=s["candidate_pressure_multiplier"],
                error_03=error[0],
                error_07=error[1],
            )
        )
    table = pd.DataFrame(rows)
    print(pd.crosstab([table.pool, table.pole_class], table.status).to_string())
    big = table[table.pole_class == "p>=0.06"]
    print(
        "  p>=0.06: minimum predicted phase error at 0.07 Hz "
        f"{big.error_07.min():.1f} deg"
    )
    gains = table.dropna(subset=["gain"]).groupby("pool").gain.median()
    print("  median candidate gain by pool:", gains.round(3).to_dict())
    level = 20 * np.log10(gains["DEPAS/Geomar hydrophone"] / gains["INSU DPG"])
    print(f"  implied nominal pressure-level excess, INSU over DEPAS: {level:.2f} dB")


def v9_brown_digitization():
    print("\nV9  Repository Brown digitization vs TiSKitPy Pressure_noise_model")
    import tiskitpy

    brown = pd.read_csv(
        ROOT / "../compy_pressure/reference_data/brown_hydroacoustic.csv"
    )
    low, high = tiskitpy.Pressure_noise_model(
        brown.frequency_hz.to_numpy(), as_freqs=True
    )
    dl = brown.low_db_pa2_hz.to_numpy() - low
    dh = brown.high_db_pa2_hz.to_numpy() - high
    beyond = int((abs(dl) > 2).sum() + (abs(dh) > 2).sum())
    print(
        f"  max |difference| low {abs(dl).max():.1f} dB, high {abs(dh).max():.1f} dB; "
        f"points beyond the 2 dB tolerance: {beyond} of {2 * len(dl)}"
    )


def v10_sensitivities():
    print("\nV10 Reference sensitivity: d ln T / d ln H, and steep apparent velocity")
    for depth in (2000, 3000, 4000, 5000):
        out = []
        for f in (0.03, 0.05, 0.07):
            try:
                hi = water_ratio([f], [3950.0], depth * 1.01, 1028, 1500)[0]
                lo = water_ratio([f], [3950.0], depth * 0.99, 1028, 1500)[0]
                out.append(f"{np.log(hi / lo) / np.log(1.01 / 0.99):.2f}")
            except ValueError:
                out.append("pole")
        steep = [
            water_ratio([f], [20000.0], depth, 1028, 1500)[0]
            / water_ratio([f], [3950.0], depth, 1028, 1500)[0]
            for f in (0.03, 0.05, 0.06)
        ]
        print(
            f"  H={depth}: dlnT/dlnH at 0.03/0.05/0.07 Hz {out}; "
            f"T(c=20 km/s)/T(c=3.95 km/s) at 0.03/0.05/0.06 Hz {np.round(steep, 3)}"
        )


def v11_published_insu_gains():
    print("\nV11 BRUIT-FM INSU candidate gains vs Aminian et al. (2025) Table 1")
    # Published DPG gain factors (multiply recorded pressure), stated +/- ~5%.
    published = {
        "RR28": 0.78,
        "RR29": 0.70,
        "RR34": 0.58,
        "RR36": 0.64,
        "RR38": 0.60,
        "RR40": 0.74,
        "RR50": 0.68,
        "RR52": 0.67,
    }
    summary = json.loads((DOCS / "yv_rayleigh_completed_summary.json").read_text())
    ours = {
        s["seed_id"].split(".")[1]: s["candidate_pressure_multiplier"]
        for s in summary["stations"]
    }
    differences = []
    for station, value in published.items():
        differences.append(100 * (ours[station] / value - 1))
        print(
            f"  {station}: published {value:.2f}, BRUIT-FM {ours[station]:.3f} "
            f"({differences[-1]:+.1f}%)"
        )
    print(
        f"  median {np.median(differences):+.1f}%, range "
        f"{min(differences):+.1f}..{max(differences):+.1f}%"
    )


def main():
    v1_acoustic_gravity()
    v2_zha_offset_models()
    v3_b01w_highpass()
    v4_c07w_states()
    v5_taper_then_differentiate()
    v6_estimator_bias()
    v7_conventions()
    v8_yv_pole_status()
    v9_brown_digitization()
    v10_sensitivities()
    v11_published_insu_gains()


if __name__ == "__main__":
    main()
