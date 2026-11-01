"""Stage A: Bayesian Vs inversion of the YV broadband compliance (ComPy 2.3).

Usage: python a_invert.py STATION [VARIANT ...]

Input is the ComPy 2.2 compliance of 2026-10-28 (gain 'compy', which matches
Aminian et al. 2025, Table 1). RR36, RR38, RR40 and RR50 use the stable
epoch before their pressure/vertical sensitivity step. Monthly stacks give
the month-to-month covariance. The Moho depth below the seafloor comes from
CRUST1.0 at the station.

Variants (each a full run; 'main' adds a prior-only run):
  main     Zha & Webb (2016) parameterization, Vp/Vs 1.85, Hamilton sediment
  brocher  Vp from Brocher (2005) instead of Vp/Vs 1.85
  sedscale sediment Vs multiplier free in [0.5, 2]
  bruitfm  pressure gain from gain_corrected_2026-10-28/gains.csv
  full     whole record instead of the stable epoch (epoch stations only)
  moho     Moho depth free within +-1 km (as in Aminian et al. 2025)
  fixedvp  Vp and density fixed at CRUST1.0 (Aminian et al. 2025; Crawford 1991)
  vpvs     one free crustal Vp/Vs in [1.7, 2.4]
  gain     main, with the pressure gain sampled within +-5 % (Aminian et al.
           2025, section 4.1.5) instead of a 5 % covariance term
  final    free crustal Vp/Vs, gain sampled per epoch, and for the four
           epoch stations the data after the sensitivity jump (Stage G, own
           Rayleigh gain) inverted together with the stable epoch; adds a
           prior-only run (final_prior)
  final_stable  final without the after-jump data
  final_h2 final with the H2 estimator S_zz/|S_pz| (noise on P; Stage H)
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COMPLIANCE = ROOT.parent / "compliance_yv_2026-10-28"
CRUST1 = COMPLIANCE / "literature" / "crust1"
sys.path.insert(0, str(ROOT.parent / "ComPy"))

import compy_inversion as ci  # noqa: E402

STATIONS = {  # latitude, longitude (Aminian et al. 2025, Table 1)
    "RR28": (-22.7152, 53.1595), "RR29": (-24.9657, 51.7488),
    "RR34": (-32.0783, 52.2114), "RR36": (-33.7018, 55.9578),
    "RR38": (-30.5650, 59.6858), "RR40": (-28.1461, 63.3020),
    "RR50": (-25.5182, 70.0222), "RR52": (-20.4723, 68.1094),
}  # fmt: skip
EPOCH = {"RR36", "RR38", "RR40", "RR50"}
SAMPLER = dict(n_chains=8, n_iter=80000, n_burn=30000, thin=10)
MIN_MONTH_WINDOWS = 30


def monthly_subsets(months, n_band):
    return np.array([m["compliance"] for m in months.values()
                     if m.get("compliance") is not None and len(m["compliance"]) == n_band
                     and m["n"] >= MIN_MONTH_WINDOWS])  # fmt: skip


def h2_result(npz):
    """The stored file with compliance replaced by the H2 estimate (Stage H)."""
    from compy_processing import gravity_corrected_compliance

    d = dict(np.load(npz))
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


def load_data(
    station, gain="compy", epoch=True, after=False, gain_rel=0.05, estimator="h1"
):
    """ComplianceData with random, monthly, gain and floor covariance terms.

    ``after`` joins the after-jump data of Stage G as a second group (epoch
    stations only); ``gain_rel=0`` drops the gain term for a sampled gain;
    ``estimator="h2"`` uses S_zz/|S_pz| (monthly subsets scaled alike).
    """
    sub = "epoch/" if epoch and station in EPOCH else ""
    npz = COMPLIANCE / "results" / station / f"{sub}{station}_compliance_{gain}.npz"
    summary = json.loads((npz.parent / f"{station}_summary.json").read_text())
    sources = [
        (npz, summary["results"][gain]["monthly"], "stable" if sub else "record")
    ]
    after_npz = ROOT / "compliance_after" / station / f"{station}_compliance_after.npz"
    if after and station in EPOCH and after_npz.exists():
        monthly = json.loads(
            (after_npz.parent / f"{station}_summary.json").read_text()
        )["monthly"]
        sources.append((after_npz, monthly, "after"))
    parts, info = [], dict(files=[], months=[], gains=[], bands=[], estimator=estimator)
    for path, months, name in sources:
        stored = dict(np.load(path))
        keep = (stored["freqs"] >= stored["band"][0] - 1e-12) & (
            stored["freqs"] <= stored["band"][1] + 1e-12
        )
        subsets = monthly_subsets(months, int(keep.sum()))
        if estimator == "h2":
            stored, scale = h2_result(path)
            subsets = subsets * scale[keep] if len(subsets) else subsets
        part = ci.ComplianceData.from_result(stored, subsets=subsets if len(subsets) >= 3 else None,
                                             gain_rel=gain_rel, floor_rel=0.01)  # fmt: skip
        part.station = station
        parts.append(part)
        info["files"].append(str(Path(path).relative_to(ROOT.parent)))
        info["months"].append(len(subsets))
        info["gains"].append(float(stored["gain"]))
        info["bands"].append(np.asarray(stored["band"]).tolist())
    data = (
        parts[0]
        if len(parts) == 1
        else ci.ComplianceData.join(parts, [s[2] for s in sources], station)
    )
    return data, info


def build_model(station, **overrides):
    column = ci.crust1_column(*STATIONS[station], CRUST1)
    settings = dict(moho_depth=column["moho_depth"], sediment_max=1000.0)
    settings.update(overrides)
    if isinstance(settings.get("reference"), str):  # "crust1": fixed Vp and density
        settings["reference"] = ci.crust1_model(column)
    return ci.CrustalModel(**settings), column


VARIANTS = {
    "main": dict(),
    "brocher": dict(model=dict(vp_vs=None)),
    "sedscale": dict(model=dict(sediment_scale=(0.5, 2.0))),
    "bruitfm": dict(data=dict(gain="bruitfm")),
    "full": dict(data=dict(epoch=False)),
    "moho": dict(model=dict(moho_shift=(-1000.0, 1000.0))),
    # Fixed Vp makes the posterior multimodal: parallel tempering, 8 replicas.
    "fixedvp": dict(
        model=dict(reference="crust1"), sampler=dict(temperatures=1.6 ** np.arange(8))
    ),
    "vpvs": dict(
        model=dict(vp_vs_range=(1.7, 2.4)),
        sampler=dict(temperatures=1.6 ** np.arange(8)),
    ),
    "gain": dict(data=dict(gain_rel=0.0), sampler=dict(gain_bounds=(0.95, 1.05))),
    "final": dict(
        model=dict(vp_vs_range=(1.7, 2.4)),
        data=dict(after=True, gain_rel=0.0),
        sampler=dict(gain_bounds=(0.95, 1.05), temperatures=1.6 ** np.arange(8)),
    ),
    "final_stable": dict(
        model=dict(vp_vs_range=(1.7, 2.4)),
        data=dict(gain_rel=0.0),
        sampler=dict(gain_bounds=(0.95, 1.05), temperatures=1.6 ** np.arange(8)),
    ),
    "final_h2": dict(
        model=dict(vp_vs_range=(1.7, 2.4)),
        data=dict(after=True, gain_rel=0.0, estimator="h2"),
        sampler=dict(gain_bounds=(0.95, 1.05), temperatures=1.6 ** np.arange(8)),
    ),
}


def run(station, variant, seed=0):
    spec = VARIANTS[variant]
    data, info = load_data(station, **spec.get("data", {}))
    model, column = build_model(station, **spec.get("model", {}))
    out = ROOT / "results" / station
    out.mkdir(parents=True, exist_ok=True)
    post = ci.sample(
        data, model, seed=seed, processes=8, **SAMPLER, **spec.get("sampler", {})
    )
    post.save(out / f"{station}_{variant}")
    record = json.loads((out / f"{station}_{variant}.json").read_text())
    record.update(input=info, crust1=column)
    (out / f"{station}_{variant}.json").write_text(
        json.dumps(record, indent=1, default=float)
    )
    s = record
    print(f"{station} {variant}: chi2 median {s['chi2']['median']:.1f} / N={s['n_data']}, "
          f"max R-hat {s['max_rhat']:.3f}, min ESS {s['min_ess']:.0f}, "
          f"acceptance {np.mean(s['acceptance']):.2f}", flush=True)  # fmt: skip
    if variant in ("main", "final"):  # tempering is pointless without a likelihood
        settings = {
            k: v for k, v in spec.get("sampler", {}).items() if k != "temperatures"
        }
        prior = ci.sample(
            data,
            model,
            seed=seed + 1,
            processes=8,
            prior_only=True,
            **SAMPLER,
            **settings,
        )
        prior.save(out / f"{station}_{'prior' if variant == 'main' else 'final_prior'}")


if __name__ == "__main__":
    station = sys.argv[1]
    for variant in sys.argv[2:] or ["main"]:
        if variant in ("full", "final_stable") and station not in EPOCH:
            continue
        run(station, variant)
