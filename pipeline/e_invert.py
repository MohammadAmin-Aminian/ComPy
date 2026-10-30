"""Stage E: Bayesian Vs inversion of one station's compliance (ComPy 2.3/2.4).

Usage: python e_invert.py STATION [VARIANT ...] [--processes N]

Input is Stage C's compliance with the ComPy Rayleigh-wave gain ('compy',
the method of Aminian et al. 2025, section 4.1.5). Where that calibration
failed, the BRUIT-FM gain is used and the gain prior is widened (see
GAIN_BOUNDS). Monthly stacks give the month-to-month covariance; the Moho
depth below the seafloor comes from CRUST1.0 at the station.

Variants:
  final     free crustal Vp/Vs in [1.7, 2.4], pressure gain sampled within the
            gain bounds, parallel tempering (8 replicas); adds a prior-only
            run (final_prior). This is the model of the YV final run.
  main      Vp/Vs 1.85 (Zha & Webb 2016) and a 5 % gain covariance term
  final_h2  final with the H2 estimator S_zz/|S_pz| (noise on the pressure)
"""

import argparse
import json
import sys

import numpy as np
from obspy import read_inventory

from config import COMPY, CRUST1, INVENTORY, RESULTS, ROOT

sys.path.insert(0, str(COMPY))

import compy_inversion as ci  # noqa: E402

SAMPLER = dict(n_chains=8, n_iter=80000, n_burn=30000, thin=10)
MIN_MONTH_WINDOWS = 30
#: Gain prior (factor on the gain used) by gain source: the paper's +-5 % for
#: the Rayleigh calibration, wider for the microseism (BRUIT-FM) gain.
GAIN_BOUNDS = {"compy": (0.95, 1.05), "bruitfm": (0.85, 1.15)}
TEMPERATURES = 1.6 ** np.arange(8)


def coordinates(station):
    s = read_inventory(str(INVENTORY)).select(station=station)[0][0]
    return s.latitude, s.longitude


def gain_label(station):
    summary = json.loads((RESULTS / station / f"{station}_summary.json").read_text())
    for label in ("compy", "bruitfm"):
        if "n_windows" in summary["results"].get(label, {}):
            return label, summary
    raise SystemExit(f"{station}: no compliance result")


def h2_result(stored):
    """The stored result with compliance scaled to the H2 estimate."""
    from compy_processing import gravity_corrected_compliance

    d = dict(stored)
    f, depth = d["all_freqs"], float(d["depth"])
    spz, spp = np.abs(d["S_pz"]), d["S_pp"]
    ok = spz > 0
    h1 = np.where(
        ok,
        gravity_corrected_compliance(
            f, depth, np.where(ok, spz / spp, 0), np.sqrt(spp)
        ),
        1,
    )
    h2 = np.where(ok, gravity_corrected_compliance(f, depth, np.where(ok, d["S_zz"] / np.where(ok, spz, 1), 0),
                                                   np.sqrt(spp)), 1)  # fmt: skip
    scale = np.interp(d["freqs"], f, h2 / h1)
    d["compliance"] = d["compliance"] * scale
    d["uncertainty"] = d["uncertainty"] * scale
    return d, scale


def load_data(station, gain_rel=0.05, estimator="h1"):
    label, summary = gain_label(station)
    path = RESULTS / station / f"{station}_compliance_{label}.npz"
    stored = dict(np.load(path))
    keep = (stored["freqs"] >= stored["band"][0] - 1e-12) & (
        stored["freqs"] <= stored["band"][1] + 1e-12
    )
    months = summary["results"][label]["monthly"]
    subsets = np.array([m["compliance"] for m in months.values()
                        if len(m["compliance"]) == int(keep.sum()) and m["n"] >= MIN_MONTH_WINDOWS])  # fmt: skip
    if estimator == "h2":
        stored, scale = h2_result(stored)
        subsets = subsets * scale[keep] if len(subsets) else subsets
    data = ci.ComplianceData.from_result(stored, subsets=subsets if len(subsets) >= 3 else None,
                                         gain_rel=gain_rel, floor_rel=0.01)  # fmt: skip
    data.station = station
    info = dict(file=str(path.relative_to(ROOT)), gain_source=label, gain=float(stored["gain"]),
                months=len(subsets), band=np.asarray(stored["band"]).tolist(), estimator=estimator,
                n_windows=summary["results"][label]["n_windows"])  # fmt: skip
    return data, info


def build_model(station, **overrides):
    column = ci.crust1_column(*coordinates(station), CRUST1)
    settings = dict(moho_depth=column["moho_depth"], sediment_max=1000.0)
    settings.update(overrides)
    return ci.CrustalModel(**settings), column


VARIANTS = {
    "final": dict(model=dict(vp_vs_range=(1.7, 2.4)), data=dict(gain_rel=0.0),
                  sampler=dict(gain=True, temperatures=TEMPERATURES)),
    "main": dict(),
    "final_h2": dict(model=dict(vp_vs_range=(1.7, 2.4)), data=dict(gain_rel=0.0, estimator="h2"),
                     sampler=dict(gain=True, temperatures=TEMPERATURES)),
}  # fmt: skip


def run(station, variant, processes, seed=0):
    spec = VARIANTS[variant]
    data, info = load_data(station, **spec.get("data", {}))
    model, column = build_model(station, **spec.get("model", {}))
    sampler = dict(spec.get("sampler", {}))
    if sampler.pop("gain", False):
        sampler["gain_bounds"] = GAIN_BOUNDS[info["gain_source"]]
    out = RESULTS / station / "inversion"
    out.mkdir(parents=True, exist_ok=True)
    post = ci.sample(data, model, seed=seed, processes=processes, **SAMPLER, **sampler)
    post.save(out / f"{station}_{variant}")
    record = json.loads((out / f"{station}_{variant}.json").read_text())
    record.update(input=info, crust1=column)
    (out / f"{station}_{variant}.json").write_text(
        json.dumps(record, indent=1, default=float)
    )
    print(f"{station} {variant}: chi2 median {record['chi2']['median']:.1f} / N={record['n_data']}, "
          f"max R-hat {record['max_rhat']:.3f}, min ESS {record['min_ess']:.0f}, "
          f"acceptance {np.mean(record['acceptance']):.2f}", flush=True)  # fmt: skip
    if variant == "final":  # tempering is pointless without a likelihood
        settings = {k: v for k, v in sampler.items() if k != "temperatures"}
        prior = ci.sample(data, model, seed=seed + 1, processes=processes, prior_only=True,
                          **SAMPLER, **settings)  # fmt: skip
        prior.save(out / f"{station}_final_prior")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("station")
    parser.add_argument("variants", nargs="*", default=["final", "main", "final_h2"])
    parser.add_argument("--processes", type=int, default=8)
    args = parser.parse_args()
    for variant in args.variants:
        run(args.station, variant, args.processes)
