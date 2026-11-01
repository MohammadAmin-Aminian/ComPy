"""Stage G: numbers for the paper-style report (results/tables/report_numbers.json).

Collects, for every station, the compliance measurement (gain, band,
threshold, windows, coherence, monthly stability), the final inversion (fit,
convergence, Vp/Vs, sediment, Vs and serpentinization at fixed depths, width
relative to the prior, kernel sensitivity depths) and the Vp/Vs = 1.85
variant, plus network-pooled profiles. The LaTeX report reads this file.
"""

import json

import numpy as np

from f_figures import (
    MIN_WINDOWS,
    compliance_summary,
    coords,
    load_posterior,
    n_windows,
    usable,
)

import compy_inversion as ci  # noqa: E402  (on sys.path through f_figures/e_invert)
from config import RESULTS, STATIONS, WORK  # noqa: E402

DEPTHS = np.array([100.0, 300.0, 500.0, 1000.0, 2000.0, 3000.0, 4000.0, 6000.0])


def q3(x, axis=0):
    return np.quantile(x, [0.05, 0.5, 0.95], axis=axis)


def station_record(s, c):
    summ = compliance_summary(s)
    lat, lon, depth = c[s]
    rec = dict(station=s, lat=lat, lon=lon, depth=depth, n_windows=n_windows(s))
    cal_path = WORK / s / "calibration" / "compy_gain.json"
    cal = json.loads(cal_path.read_text()) if cal_path.exists() else {}
    rec.update(gain_compy=cal.get("gain"), gain_compy_std=cal.get("gain_std"),
               gain_events=len(cal.get("events_used", [])), calibration_error=cal.get("error"))  # fmt: skip
    if not summ:
        rec["status"] = "no compliance (stage B: no usable segments)"
        return rec
    rec["gain_bruitfm"] = summ.get("gains", {}).get("bruitfm")
    label = next(
        (
            lab
            for lab in ("compy", "bruitfm", "nominal")
            if "n_windows" in summ["results"].get(lab, {})
        ),
        None,
    )
    if label is None:
        rec["status"] = "no compliance (gates)"
        return rec
    r = summ["results"][label]
    d = np.load(RESULTS / s / f"{s}_compliance_{label}.npz")
    m = (d["freqs"] >= d["band"][0] - 1e-9) & (d["freqs"] <= d["band"][1] + 1e-9)
    full = d["compliance"][m]
    months = [np.array(v["compliance"]) for v in r["monthly"].values()
              if len(v["compliance"]) == m.sum() and v["n"] >= 30]  # fmt: skip
    dev = [float(np.median(x / full) - 1) for x in months]
    rec.update(gain_source=label, gain_used=r["gain"], band=r["band_hz"], threshold=r["threshold"],
               coherence=r["median_coherence_in_band"], rel_error=r["median_relative_uncertainty"],
               months=len(months), monthly_dev=[min(dev), max(dev)] if dev else None,
               eta=dict(freqs=d["freqs"][m].tolist(), eta=full.tolist()))  # fmt: skip
    if rec["n_windows"] < MIN_WINDOWS:
        rec["status"] = f"rejected: {rec['n_windows']} windows < {MIN_WINDOWS}"
        return rec
    if not usable(s):
        rec["status"] = "compliance only (inversion pending or failed)"
        return rec
    rec["status"] = "inverted" + (" (uncalibrated gain)" if label == "nominal" else "")
    post, column = load_posterior(s, "final")
    prior, _ = load_posterior(s, "final_prior")
    main, _ = load_posterior(s, "main")
    j = json.loads((RESULTS / s / "inversion" / f"{s}_final.json").read_text())
    flat, names = post.flat(), post.names
    vs = post.vs_draws(DEPTHS, 4000)
    rec.update(chi2=j["chi2"]["median"], n_data=j["n_data"], rhat=j["max_rhat"], ess=j["min_ess"],
               vp_vs=q3(flat[:, names.index("vp_vs")]).tolist(),
               sediment=q3(flat[:, names.index("sediment_thickness")]).tolist(),
               crust1_sediment=column["sediment_thickness"], crust1_moho=column["moho_depth"],
               vs=q3(vs).tolist(), serp=q3(ci.serpentinization(vs)).tolist())  # fmt: skip
    if "gain" in names:
        rec["gain_sampled"] = (q3(flat[:, names.index("gain")]) * r["gain"]).tolist()
    if prior is not None:
        qp = q3(prior.vs_draws(DEPTHS, 4000))
        qv = np.array(rec["vs"])
        rec["width_ratio"] = ((qv[2] - qv[0]) / (qp[2] - qp[0])).tolist()
    if main is not None:
        jm = json.loads((RESULTS / s / "inversion" / f"{s}_main.json").read_text())
        rec.update(
            vs_main=q3(main.vs_draws(DEPTHS, 4000)).tolist(),
            chi2_main=jm["chi2"]["median"],
        )
    model = post.model.layers(np.median(flat, axis=0))
    z, kern = ci.vs_kernels(post.data.depth, post.data.freqs[:1], model)
    rec["sens50"], rec["sens90"] = (
        float(ci.sensitivity_depth(z, kern[0], f)) for f in (0.5, 0.9)
    )
    return rec


def main():
    c = coords()
    recs = [station_record(s, c) for s in STATIONS]
    inv = [s for s, r in zip(STATIONS, recs) if r["status"].startswith("inverted")]
    calibrated = [s for s, r in zip(STATIONS, recs) if r["status"] == "inverted"]
    pooled = {}
    if calibrated:
        draws = np.concatenate(
            [load_posterior(s, "final")[0].vs_draws(DEPTHS, 1500) for s in calibrated]
        )
        pooled["vs"] = q3(draws).tolist()
        pooled["serp"] = q3(ci.serpentinization(draws)).tolist()
        mains = [load_posterior(s, "main")[0] for s in calibrated]
        mains = [m for m in mains if m is not None]
        if mains:
            pooled["vs_main"] = q3(
                np.concatenate([m.vs_draws(DEPTHS, 1500) for m in mains])
            ).tolist()
            pooled["n_main"] = len(mains)
        crust1 = []
        for s in calibrated:
            m = ci.crust1_model(load_posterior(s, "final")[1])
            tops = np.concatenate([[0.0], np.cumsum(m[:-1, 0])])
            crust1.append(m[np.searchsorted(tops, DEPTHS, side="right") - 1, 3])
        pooled["crust1_vs"] = np.mean(crust1, axis=0).tolist()
    out = dict(
        depths=DEPTHS.tolist(),
        stations=recs,
        inverted=inv,
        calibrated=calibrated,
        pooled=pooled,
    )
    (RESULTS / "tables").mkdir(parents=True, exist_ok=True)
    (RESULTS / "tables" / "report_numbers.json").write_text(
        json.dumps(out, indent=1, default=float)
    )
    for r in recs:
        print(r["station"], r["status"])


if __name__ == "__main__":
    main()
