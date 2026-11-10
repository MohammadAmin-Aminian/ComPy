"""Vs inversion of a compliance curve with compy.inversion.auto (ComPy 3).

Usage (repository root):
    python examples/4_Vs_inversion.py                      # offline, synthetic
    python examples/4_Vs_inversion.py RR52_compliance_compy.npz \
        --lat -20.4723 --lon 68.1094 --crust1 path/to/crust1.0

Without arguments a compliance curve is computed from a known model (an oceanic
crust slower than its reference in the top kilometre, under 120 m of sediment),
noise with a realistic covariance is added, and the curve is inverted. The
output shows what the data recover and what they leave to the prior.

With a stage C result (``.npz``) the script does what ``compy-invert`` does:
the covariance saved by stage C, the CRUST1.0 model of the station as the
reference, and the default priors (docs/INVERSION.md).

About 2 minutes on 8 cores for the synthetic case (shorter with --quick).
"""

import argparse

import numpy as np

from compy.inversion import auto, dram
from compy.inversion.cli import read_compliance
from compy.inversion.errors import default_covariance

#: A normal oceanic crust: rows [thickness m, density kg/m3, Vp m/s, Vs m/s].
REFERENCE = np.array([
    [1170.0, 2400.0, 4700.0, 2540.0],
    [1170.0, 2760.0, 6300.0, 3590.0],
    [4560.0, 2960.0, 6900.0, 3940.0],
    [0.0, 3190.0, 7730.0, 4310.0],
])  # fmt: skip


def synthetic():
    """Data from a known model, with noise drawn from the covariance used to invert."""
    depth, freqs = 2918.0, np.arange(6, 15) * 1.220703125e-3
    model = auto.SeafloorModel(REFERENCE)
    truth = model.reference_theta(sediment=120.0)
    first = model.names.index("vs_crust_1")
    truth[first : first + 4] *= [0.75, 0.8, 0.85, 0.9]
    clean = dram.forward(depth, freqs, model.layers(truth))
    cov_rel = default_covariance(freqs, clean, 0.0015 * clean)
    rng = np.random.default_rng(1)
    observed = clean + np.linalg.cholesky(
        cov_rel * np.outer(clean, clean)
    ) @ rng.standard_normal(len(freqs))
    data = auto.compliance_data(
        freqs, observed, depth, cov_rel=cov_rel, station="synthetic"
    )
    return data, model, truth


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("npz", nargs="?")
    parser.add_argument("--lat", type=float)
    parser.add_argument("--lon", type=float)
    parser.add_argument("--crust1")
    parser.add_argument("--output", default="vs_inversion")
    parser.add_argument(
        "--quick", action="store_true", help="short chains (for a first look)"
    )
    parser.add_argument("--processes", type=int, default=None)
    args = parser.parse_args()
    truth = None
    if args.npz is None:
        data, model, truth = synthetic()
    else:
        if args.crust1 is None or args.lat is None or args.lon is None:
            parser.error("a real station needs --lat, --lon and --crust1")
        freqs, eta, cov_rel, depth, _ = read_compliance(args.npz)
        data = auto.compliance_data(
            freqs, eta, depth, cov_rel=cov_rel, station=args.npz
        )
        model = auto.model_from_crust1(
            dram.crust1_column(args.lat, args.lon, args.crust1)
        )
    sampler = dict(n_iter=15000, n_burn=5000, max_rounds=1) if args.quick else {}
    result = auto.invert(
        data, model, gain_sigma=0.05, processes=args.processes, **sampler
    )

    quality = result.quality()
    print(f"quality: {quality['label']}; chi-square {quality['chi2']:.1f} for "
          f"{quality['n_data']} data; max R-hat {quality['max_rhat']:.3f}")  # fmt: skip
    for reason in quality["reasons"]:
        print("  note:", reason)
    summary = result.summary()

    def fmt(x):
        return f"{x:.0f}" if abs(x) >= 10 else f"{x:.3f}"

    header = "\nparameter               median   16 % to 84 %          width/prior"
    print(header + ("    true" if truth is not None else ""))
    for i, name in enumerate(model.names):
        p = summary["parameters"][name]
        line = f"  {name:20s} {fmt(p['median']):>8s}   {fmt(p['p16'])} to {fmt(p['p84'])}".ljust(
            58
        )
        line += f"{summary['resolution'][name]:5.2f}"
        print(line + (f"   {fmt(truth[i]):>8s}" if truth is not None else ""))
    profile = result.profile(np.array([100.0, 300.0, 700.0, 1500.0, 3000.0, 5000.0]))
    print("\ndepth below seafloor   Vs median (5-95 %)        resolved")
    for k, z in enumerate(profile["z"]):
        print(f"  {z:6.0f} m             {profile['vs'][2, k]:5.0f} m/s ({profile['vs'][0, k]:.0f}"
              f" to {profile['vs'][4, k]:.0f})".ljust(52) + ("yes" if profile["resolved"][k] else "no"))  # fmt: skip
    result.save(args.output)
    try:
        import matplotlib

        matplotlib.use("Agg")
        from compy.plotting.auto import station_figure

        station_figure(result, args.output + ".png")
        figure = f", {args.output}.png"
    except ImportError:
        figure = ""
    print(
        f"\nsaved {args.output}.npz, {args.output}.json, {args.output}_profile.csv{figure}"
    )


if __name__ == "__main__":
    main()
