"""Stage E, alternatives: other inversions of one station's compliance.

Usage: python e_alternative_inversions.py STATION [METHOD ...] [--processes N]

Methods, on the same stage C compliance and gain as e_invert.py (no gain is sampled):
  minimum_structure  after Crawford (2004), on CRUST1.0 Vp and density, from two
                     starts, keeping the smoother one if both fit, otherwise the better fit
  zha_webb_2016      Monte Carlo with the parameterization of Zha & Webb (2016)
  zheng_2024         Monte Carlo with the parameterization of Zheng et al. (2024)

Results go to results/<STA>/inversion_alternative. Monte Carlo errors are the standard
deviation of the monthly curves; minimum_structure uses the marginal errors of e_invert.py.
These methods test how much a result depends on the method: they do not replace stage E.
"""

import argparse
import json
import sys

import numpy as np
from config import COMPY, CRUST1, RESULTS
from e_invert import coordinates, gain_label, load_data

sys.path.insert(0, str(COMPY))

from compy.inversion.dram import crust1_column, crust1_model  # noqa: E402
from compy.inversion.literature_mcmc import (  # noqa: E402
    Zheng2024Model,
    ZhaWebb2016Model,
    monthly_standard_deviation,
    sample_literature_model,
)
from compy.inversion.minimum_structure import (  # noqa: E402
    crawford_thin_layers,
    invert_minimum_structure,
)

METHODS = ("minimum_structure", "zha_webb_2016", "zheng_2024")
FLAT_START = 2000.0  # m/s, the second starting model of the minimum-structure inversion
#: Grid of the saved Vs profiles [m].
DEPTHS = np.arange(0.0, 8001.0, 50.0)


def rock_column(column):
    """Moho, tops and Vp of the crystalline and sedimentary rows of a CRUST1.0 column."""
    thick = np.array([layer[1] for layer in column["layers"]])
    return dict(
        moho_depth=column["moho_depth"],
        vp_tops=np.r_[0.0, np.cumsum(thick)[:-1]],
        vp_values=np.array([layer[2] for layer in column["layers"]]),
        vp_mantle=column["mantle"][0],
    )


def monthly_data(station):
    """Frequency, compliance and monthly standard deviation of a station's stage C result."""
    label, summary = gain_label(station)
    stored = dict(np.load(RESULTS / station / f"{station}_compliance_{label}.npz"))
    keep = (stored["freqs"] >= stored["band"][0] - 1e-12) & (
        stored["freqs"] <= stored["band"][1] + 1e-12
    )
    curves = [
        m["compliance"]
        for m in summary["results"][label]["monthly"].values()
        if len(m["compliance"]) == int(keep.sum())
    ]
    return (
        float(stored["depth"]),
        stored["freqs"][keep],
        stored["compliance"][keep],
        monthly_standard_deviation(curves),
    )


def run_minimum_structure(station, column, out):
    """Minimum-structure inversion from two starts."""
    data, info = load_data(station, gain_rel=0.0)
    grid, edges = crawford_thin_layers(crust1_model(column))
    args = (data.eta, data.freqs, data.depth, grid, data.sigma)
    first = invert_minimum_structure(*args)
    second = invert_minimum_structure(*args, start=np.full(len(grid) - 1, FLAT_START))
    both_fit = first.reached_target and second.reached_target
    swap = (
        (second.roughness < first.roughness) if both_fit else (second.chi2 < first.chi2)
    )
    best, other = (second, first) if swap else (first, second)
    np.savez(
        out / f"{station}_minimum_structure.npz",
        model=best.model,
        other_start_vs=other.model[:-1, 3],
        edges=edges,
        prediction=best.prediction,
        freqs=data.freqs,
        eta=data.eta,
        sigma=data.sigma,
        depth=data.depth,
    )
    record = dict(
        chi2=best.chi2,
        n_data=len(data.eta),
        reached_target=best.reached_target,
        roughness=best.roughness,
        iterations=best.iterations,
        other_start_chi2=other.chi2,
        last_layer_vs=float(best.model[-2, 3]),
        half_space_vs=float(best.model[-1, 3]),
        gain_source=info["gain_source"],
    )
    (out / f"{station}_minimum_structure.json").write_text(json.dumps(record, indent=1))
    print(
        f"{station} minimum_structure: chi2 {best.chi2:.2f} / {len(data.eta)}, "
        f"{'reached' if best.reached_target else 'did NOT reach'} the target, "
        f"last layer {best.model[-2, 3]:.0f} m/s, half-space {best.model[-1, 3]:.0f} m/s",
        flush=True,
    )


def run_monte_carlo(station, method, column, out, processes):
    """Monte Carlo inversion with a published parameterization, and its prior."""
    depth, freqs, eta, sigma = monthly_data(station)
    if method == "zha_webb_2016":
        model = ZhaWebb2016Model(**rock_column(column))
    else:
        model = Zheng2024Model(column["moho_depth"])
    result = sample_literature_model(
        model, depth, freqs, eta, sigma, processes=processes
    )
    prior = sample_literature_model(
        model, depth, freqs, eta, sigma, processes=processes, prior_only=True, seed=12,
        n_iter=100000, n_burn=20000,
    )  # fmt: skip
    result.save(out / f"{station}_{method}")
    np.savez_compressed(
        out / f"{station}_{method}_profiles.npz",
        z=DEPTHS,
        vs_draws=result.vs_draws(DEPTHS),
        prior_quantiles=prior.vs_quantiles(DEPTHS, max_draws=4000),
    )
    s = result.summary()
    print(
        f"{station} {method}: chi2 {s['chi2_median_retained']:.2f} / {s['n_data']}, "
        f"retained {s['retained']} / {s['total']}, max R-hat {s['max_rhat']:.3f}, "
        f"min ESS {s['min_ess']:.0f}, Moho jump {s['moho_jump']['median']:.0f} m/s",
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("station")
    parser.add_argument("methods", nargs="*", default=list(METHODS))
    parser.add_argument("--processes", type=int, default=8)
    args = parser.parse_args()
    column = crust1_column(*coordinates(args.station), CRUST1)
    folder = RESULTS / args.station / "inversion_alternative"
    folder.mkdir(parents=True, exist_ok=True)
    for name in args.methods:
        if name == "minimum_structure":
            run_minimum_structure(args.station, column, folder)
        else:
            run_monte_carlo(args.station, name, column, folder, args.processes)
