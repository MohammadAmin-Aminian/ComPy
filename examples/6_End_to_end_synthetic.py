"""From two time series to a shear-velocity profile, without any data file.

Usage (any folder; ComPy 3.2 or later installed):
    python 6_End_to_end_synthetic.py --output worked_example

The script makes fourteen days of seafloor pressure and vertical displacement
at 1.25 Hz for a known crust under 2918 m of water, as stage B of the pipeline
would deliver them, and then does what stages C and E do with real records:

1. cut the record into 35-min windows and select them (``stacked_compliance``
   with the gates of Aminian et al. 2025);
2. choose the band where pressure and displacement are coherent
   (``coherent_band``) and stack the compliance with its error;
3. measure the covariance of the stack from its own windows
   (``stack_covariance``);
4. invert for shear velocity (``auto.invert``) and report what the data resolve.

The synthetic record. Pressure is Gaussian noise with the spectrum of
infragravity waves seen at the seafloor (a plateau cut by 1/cosh^2(kH)), its
level drifting with the sea state, plus sensor noise. The seismometer records
the deformation of the known crust under that pressure, the gravitational
attraction of the waves and the free-air term (the terms the estimator
removes), sensor noise, and, in one window out of twelve, a transient that the
gates must reject. Nothing is taken from the estimator to build the record
except the physical constants.

Outputs in --output: the two spectra, the compliance curve, the Vs profile and
the parameter table as CSV, the posterior draws (``.npz``), a summary
(``.json``) and one figure. About one minute on eight cores.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
from obspy import Stream, Trace, UTCDateTime

import compy
from compy.compliance import GJI2025, coherent_band, spectra, stacked_compliance
from compy.inversion import auto, dram
from compy.inversion.errors import attraction_ratio, stack_covariance, usable_bins
from compy.numerics import FREE_AIR_GRADIENT, GRAVITY, wavenumber

DEPTH = 2918.0  # water depth [m]
RATE = 1.25  # sampling rate after stage A [Hz]
DAYS = 14
WINDOW_S = 2100.0  # 35 min
NSEG = 1024
#: A normal oceanic crust: rows [thickness m, density kg/m3, Vp m/s, Vs m/s].
REFERENCE = np.array([
    [1170.0, 2400.0, 4700.0, 2540.0],
    [1170.0, 2760.0, 6300.0, 3590.0],
    [4560.0, 2960.0, 6900.0, 3940.0],
    [0.0, 3190.0, 7730.0, 4310.0],
])  # fmt: skip


def true_model():
    """The crust to recover: slower than the reference in its top, 120 m of sediment."""
    model = auto.SeafloorModel(REFERENCE)
    truth = model.reference_theta(sediment=120.0)
    first = model.names.index("vs_crust_1")
    truth[first : first + 4] *= [0.75, 0.8, 0.85, 0.9]
    return model, truth


def shaped_noise(rng, n, psd):
    """Gaussian noise of ``n`` samples with the one-sided spectrum ``psd(f)`` [unit^2/Hz]."""
    f = np.fft.rfftfreq(n, 1 / RATE)
    white = np.fft.rfft(rng.standard_normal(n))
    amplitude = np.sqrt(np.maximum(psd(f), 0.0) * RATE / 2)
    amplitude[0] = 0.0
    return np.fft.irfft(white * amplitude, n), f


def synthetic_record(seed=1):
    """Pressure [Pa] and vertical displacement [m] of the known crust, as a Stream."""
    rng = np.random.default_rng(seed)
    model, truth = true_model()
    n = int(DAYS * 86400 * RATE)
    n -= n % 2

    def wave_number(f):
        k = np.zeros_like(f)
        k[1:] = wavenumber(2 * np.pi * f[1:], DEPTH)
        return k

    def pressure_psd(f):  # infragravity waves at the seafloor [Pa^2/Hz]
        return 2000.0 / np.cosh(np.minimum(wave_number(f) * DEPTH, 50.0)) ** 2

    waves, f = shaped_noise(rng, n, pressure_psd)
    # Sea state: the wave pressure drifts by a factor of about two over days.
    drift, _ = shaped_noise(rng, n, lambda f: np.where((f > 0) & (f < 2e-5), 1.0, 0.0))
    waves *= np.exp(0.35 * drift / drift.std())
    # Displacement per unit pressure of the crust, with the terms a seismometer adds.
    # Above 50 mHz no wave pressure reaches this depth (1/cosh^2(kH) < 1e-12),
    # so the transfer function is not needed there.
    grid = np.geomspace(5e-4, 0.05, 200)
    eta = dram.forward(DEPTH, grid, model.layers(truth))
    k = wavenumber(2 * np.pi * grid, DEPTH)
    omega2 = (2 * np.pi * grid) ** 2
    ratio = (
        eta * (omega2 + FREE_AIR_GRADIENT) / k - attraction_ratio(grid, DEPTH)
    ) / omega2
    transfer = np.interp(f, grid, ratio, left=0.0, right=0.0)
    displacement = np.fft.irfft(np.fft.rfft(waves) * transfer, n)
    # Sensor noise: pressure 0.05 Pa^2/Hz; acceleration -178 dB re (m/s^2)^2/Hz.
    pressure = waves + shaped_noise(rng, n, lambda f: np.full_like(f, 0.05))[0]

    def displacement_noise(f):
        out = np.zeros_like(f)
        out[1:] = 10 ** (-178 / 10) / (2 * np.pi * f[1:]) ** 4
        return np.where(f > 5e-4, out, 0.0)

    displacement += shaped_noise(rng, n, displacement_noise)[0]
    # One window in twelve carries a transient on the vertical (an earthquake
    # or a current gust), at a random place of the record.
    width = int(WINDOW_S * RATE)
    hit = rng.choice(n // width, size=(n // width) // 12, replace=False)
    t = np.arange(width) / RATE
    for i in hit:
        burst = 3e-7 * np.exp(-t / 300.0) * np.sin(2 * np.pi * 0.012 * t)
        displacement[i * width : (i + 1) * width] += burst * rng.uniform(0.5, 2.0)
    start = UTCDateTime(2020, 1, 1)
    stats = dict(network="XX", station="SYN", location="00", sampling_rate=RATE,
                 starttime=start)  # fmt: skip
    stream = Stream([Trace(displacement, dict(stats, channel="BHZ")),
                     Trace(pressure, dict(stats, channel="BDH"))])  # fmt: skip
    return stream, model, truth, len(hit)


def measure_compliance(stream):
    """Stages C: windows, gates, band, stacked compliance and its covariance."""
    f_c = float(np.sqrt(GRAVITY / (2 * np.pi * DEPTH)))  # infragravity cut-off
    windows = compy.split_stream(stream, WINDOW_S)
    first = stacked_compliance(windows, DEPTH, 1.0, (0.007, 0.8 * f_c), GJI2025, NSEG)
    band = coherent_band(first.all_freqs, first.mean_coherence_all, 0.9,
                         (0.004, 0.8 * f_c))  # fmt: skip
    result = stacked_compliance(windows, DEPTH, 1.0, band, GJI2025, NSEG,
                                output_band=(0.002, 0.04))  # fmt: skip
    accepted = {float(t) for t in result.window_starts}
    starts, pp, zp, zz_all, pp_all, coh_all = [], [], [], [], [], []
    for window in windows:
        f, zz, p, x, coh = spectra(window, NSEG)
        zz_all.append(zz * (2 * np.pi * np.maximum(f, 1e-9)) ** 4)
        pp_all.append(p)
        coh_all.append(coh)
        if float(window[0].stats.starttime) in accepted:
            in_band = (f >= band[0] - 1e-12) & (f <= band[1] + 1e-12)
            starts.append(float(window[0].stats.starttime))
            pp.append(p[in_band])
            zp.append(x[in_band])
    f_band = f[in_band]
    errors = stack_covariance(
        f_band, np.array(pp), np.array(zp), np.array(starts), DEPTH
    )
    keep = usable_bins(errors["sd_rel"])
    raw = dict(freq=f, acceleration=np.array(zz_all), pressure=np.array(pp_all),
               coherence=np.array(coh_all))  # fmt: skip
    return result, f_band, errors, keep, raw, len(windows)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", default="worked_example", help="output folder")
    parser.add_argument(
        "--seed", type=int, default=1, help="seed of the synthetic record"
    )
    parser.add_argument("--processes", type=int, default=None)
    parser.add_argument("--no-figure", action="store_true")
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    clock = time.time()

    stream, model, truth, n_transients = synthetic_record(args.seed)
    result, f_band, errors, keep, raw, n_windows = measure_compliance(stream)
    t_compliance = time.time() - clock
    eta_true = dram.forward(DEPTH, result.freqs, model.layers(truth))
    in_band = (result.freqs >= result.band[0] - 1e-12) & (
        result.freqs <= result.band[1] + 1e-12
    )
    misfit = result.compliance[in_band] / eta_true[in_band] - 1
    print(f"windows: {n_windows} of 35 min, {result.n_windows} accepted, rejected {result.rejected} "
          f"({n_transients} carried a transient)")  # fmt: skip
    print(f"band: {1e3 * result.band[0]:.1f} to {1e3 * result.band[1]:.1f} mHz; stacked coherence "
          f"{np.median(result.coherence[in_band]):.3f}")  # fmt: skip
    print(f"compliance against the known curve: median {100 * np.median(misfit):+.2f} %, "
          f"largest {100 * np.max(np.abs(misfit)):.2f} %; median error bar "
          f"{100 * np.median(result.uncertainty[in_band] / result.compliance[in_band]):.2f} %")  # fmt: skip

    eta = np.interp(f_band[keep], result.freqs, result.compliance)
    data = auto.compliance_data(f_band[keep], eta, DEPTH,
                                cov_rel=errors["cov_rel"][np.ix_(keep, keep)],
                                station="synthetic")  # fmt: skip
    clock = time.time()
    inversion = auto.invert(data, model, gain_sigma=0.05, processes=args.processes)
    t_inversion = time.time() - clock
    quality = inversion.quality()
    print(f"inversion: {quality['label']}; chi-square {quality['chi2']:.1f} for "
          f"{quality['n_data']} data; max R-hat {quality['max_rhat']:.3f}")  # fmt: skip
    z = np.arange(0.0, 9001.0, 50.0)
    profile = inversion.profile(z)
    vs_true = model.vs_at(truth, z)
    resolved = profile["resolved"]
    inside = (vs_true >= profile["vs"][0]) & (vs_true <= profile["vs"][4])
    print("depth below seafloor   Vs median (5-95 %)        true    resolved")
    for depth in (300.0, 700.0, 1500.0, 3000.0, 5000.0):
        i = int(np.argmin(np.abs(z - depth)))
        print(f"  {depth:6.0f} m             {profile['vs'][2, i]:5.0f} ({profile['vs'][0, i]:.0f} to "
              f"{profile['vs'][4, i]:.0f})".ljust(50) + f"{vs_true[i]:5.0f}    "
              + ("yes" if resolved[i] else "no"))  # fmt: skip
    if resolved.any():
        print(f"resolved from {z[resolved].min():.0f} to {z[resolved].max():.0f} m; the known Vs is "
              f"inside the 5-95 % interval at {100 * inside[resolved].mean():.0f} % of these depths")  # fmt: skip
    print(
        f"time: {t_compliance:.0f} s for the compliance, {t_inversion:.0f} s for the inversion"
    )

    # ---- files ----
    def save(name, header, columns):
        np.savetxt(out / name, np.column_stack(columns), delimiter=",", header=header,
                   comments="", fmt="%.8g")  # fmt: skip

    sel = (raw["freq"] >= 0.002) & (raw["freq"] <= 0.1)
    db = {key: 10 * np.log10(np.quantile(raw[key][:, sel], (0.1, 0.5, 0.9), axis=0))
          for key in ("pressure", "acceleration")}  # fmt: skip
    save("spectra.csv",
         "freq_hz,pressure_p10_db_re_1Pa2_per_Hz,pressure_p50_db,pressure_p90_db,"
         "acceleration_p10_db_re_1m2s-4_per_Hz,acceleration_p50_db,acceleration_p90_db,"
         "mean_window_coherence",
         [raw["freq"][sel], *db["pressure"], *db["acceleration"],
          raw["coherence"][:, sel].mean(0)])  # fmt: skip
    save("compliance.csv",
         "freq_hz,compliance_1_per_pa,uncertainty_1_per_pa,true_compliance_1_per_pa,"
         "coherence,in_band",
         [result.freqs, result.compliance, result.uncertainty, eta_true, result.coherence,
          in_band.astype(int)])  # fmt: skip
    save("vs_profile.csv",
         "depth_below_seafloor_m,vs_p05_m_s,vs_p16_m_s,vs_p50_m_s,vs_p84_m_s,vs_p95_m_s,"
         "prior_p05_m_s,prior_p50_m_s,prior_p95_m_s,true_vs_m_s,resolved",
         [z, *profile["vs"], profile["vs_prior"][0], profile["vs_prior"][2],
          profile["vs_prior"][4], vs_true, resolved.astype(int)])  # fmt: skip
    predicted = np.quantile(inversion.posterior.pred.reshape(-1, len(data.freqs)),
                            (0.05, 0.5, 0.95), axis=0)  # fmt: skip
    save("fit.csv", "freq_hz,observed_1_per_pa,predicted_p05,predicted_p50,predicted_p95",
         [data.freqs, data.eta, *predicted])  # fmt: skip
    summary = inversion.summary()
    with open(out / "parameters.csv", "w") as fh:
        fh.write("parameter,true,median,p16,p84,posterior_width_over_prior_width\n")
        for i, name in enumerate(model.names):
            p = summary["parameters"][name]
            fh.write(f"{name},{truth[i]:.6g},{p['median']:.6g},{p['p16']:.6g},{p['p84']:.6g},"
                     f"{summary['resolution'][name]:.3f}\n")  # fmt: skip
    inversion.save(out / "inversion")
    record = dict(
        compy_version=compy.__version__, water_depth_m=DEPTH, days=DAYS, seed=args.seed,
        windows=n_windows, windows_accepted=int(result.n_windows),
        windows_rejected={k: int(v) for k, v in result.rejected.items()},
        windows_with_transient=int(n_transients),
        band_hz=[float(v) for v in result.band],
        stacked_coherence_median=float(np.median(result.coherence[in_band])),
        compliance_median_difference_percent=float(100 * np.median(misfit)),
        compliance_largest_difference_percent=float(100 * np.max(np.abs(misfit))),
        compliance_median_error_percent=float(
            100 * np.median(result.uncertainty[in_band] / result.compliance[in_band])),
        data_in_inversion=int(len(data.freqs)), quality=quality["label"],
        chi2=float(quality["chi2"]), max_rhat=float(quality["max_rhat"]),
        resolved_from_m=float(z[resolved].min()) if resolved.any() else None,
        resolved_to_m=float(z[resolved].max()) if resolved.any() else None,
        true_vs_inside_90_percent_interval_where_resolved=(
            float(inside[resolved].mean()) if resolved.any() else None),
        seconds_compliance=round(t_compliance, 1), seconds_inversion=round(t_inversion, 1),
    )  # fmt: skip
    (out / "summary.json").write_text(json.dumps(record, indent=1) + "\n")
    if not args.no_figure:
        try:
            import matplotlib

            matplotlib.use("Agg")
            from compy.plotting.auto import station_figure

            station_figure(inversion, out / "inversion.png")
        except ImportError:
            pass
    print(f"written to {out}/: spectra.csv, compliance.csv, fit.csv, vs_profile.csv, "
          "parameters.csv, inversion.npz, inversion.json, summary.json")  # fmt: skip


if __name__ == "__main__":
    main()
