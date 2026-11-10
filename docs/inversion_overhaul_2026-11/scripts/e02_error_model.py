"""E2: what the errors of a stacked compliance curve really are.

For every station of the all-station plan, the accepted windows of the wide-gate
stage C result are taken from the per-window spectra cache (13_audit) and restacked:

* ``bp``      Bendat & Piersol error stored by stage C;
* ``iid``     bootstrap over windows (windows resampled independently);
* ``block_d`` bootstrap over blocks of d days (1, 5, 30): adds whatever changes
              more slowly than d days (sea state, instrument, selection);
* correlation between frequency bins in each case, and the share of the block
  covariance carried by a common factor (all bins moving together) and by a
  tilt (a slope in log frequency);
* the H2/H1 ratio (1/coherence^2 of the stack).

Nothing in the archive is written. Output: one CSV row per station and a JSON
with the band-resolved curves.

Run: PYTHONPATH=<worktree> python e02_error_model.py <outdir> [--workers 10] [--stations NET.STA ...]
"""

import argparse
import csv
import json
import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np

ARCHIVE = Path.home() / "Desktop/ComPy"
CACHE = ARCHIVE / "13_audit/spectra_cache"
STAGE_C = ARCHIVE / "06_compliance_wide"
PLAN = ARCHIVE / "12_inversion_all_stations/tables/station_plan.csv"
N_BOOT = 400
BLOCK_DAYS = (1, 5, 30)


def _eta(f, k, attraction, zp_sum, pp_sum):
    from compy.numerics import FREE_AIR_GRADIENT

    omega2 = (2 * np.pi * f) ** 2
    ratio = np.abs(zp_sum) / pp_sum
    return k * (omega2 * ratio + attraction) / (omega2 + FREE_AIR_GRADIENT)


def station(args):
    net, sta, label, epoch_end = args
    try:
        return _station(net, sta, label, epoch_end)
    except Exception as exc:  # a failure is a row, not a stop
        return dict(sid=f"{net}.{sta}", error=f"{type(exc).__name__}: {exc}")


def _station(net, sta, label, epoch_end):
    from compy.numerics import GRAVITATIONAL_CONSTANT, GRAVITY, wavenumber

    folder = STAGE_C / net / "results" / sta
    if epoch_end:
        folder = folder / "epoch"
    res = dict(np.load(folder / f"{sta}_compliance_{label}.npz"))
    cache = np.load(CACHE / net / f"{sta}.npz")
    f_all, starts = cache["f"], cache["starts"]
    gain, depth = float(res["gain"]), float(res["depth"])
    lo, hi = res["band"]
    hi = min(hi, float(res["freqs"].max()))  # stage C stores at most 40 mHz (OUTPUT_BAND)
    band = (f_all >= lo - 1e-12) & (f_all <= hi + 1e-12)
    f = f_all[band]
    order = np.argsort(starts)
    idx = order[np.searchsorted(starts[order], res["window_starts"])]
    if not np.allclose(starts[idx], res["window_starts"], atol=1e-3):
        raise RuntimeError("accepted windows not found in the cache")
    t = starts[idx]
    zz = cache["zz"][idx][:, band]
    pp = cache["pp"][idx][:, band] * gain**2
    zp = cache["zp"][idx][:, band] * gain
    n = len(t)
    k = wavenumber(2 * np.pi * f, depth)
    kh = k * depth
    attraction = 2 * np.pi * (GRAVITATIONAL_CONSTANT / GRAVITY) * (1 + np.exp(-2 * kh)) / 2
    eta = _eta(f, k, attraction, zp.sum(0), pp.sum(0))
    stored = res["compliance"][(res["freqs"] >= lo - 1e-12) & (res["freqs"] <= hi + 1e-12)]
    check = float(np.max(np.abs(eta / stored - 1)))
    bp = res["uncertainty"][(res["freqs"] >= lo - 1e-12) & (res["freqs"] <= hi + 1e-12)] / stored
    coh2 = np.abs(zp.sum(0)) ** 2 / (zz.sum(0) * pp.sum(0))
    rng = np.random.default_rng(abs(hash(sta)) % 2**32)

    def boot(groups):
        """Bootstrap over groups of windows; returns relative deviations (B, nf)."""
        ids = np.unique(groups)
        # per-group sums
        inv = np.searchsorted(ids, groups)
        gz = np.zeros((len(ids), zp.shape[1]), complex)
        gp = np.zeros((len(ids), pp.shape[1]))
        np.add.at(gz, inv, zp)
        np.add.at(gp, inv, pp)
        w = rng.multinomial(len(ids), np.full(len(ids), 1 / len(ids)), size=N_BOOT)
        return _eta(f, k, attraction, w @ gz, w @ gp) / eta - 1, len(ids)

    out = dict(sid=f"{net}.{sta}", label=label, depth=depth, n_windows=n, n_freq=int(band.sum()),
               f_lo=float(f[0]), f_hi=float(f[-1]), restack_check=check,
               span_days=float((t[-1] - t[0]) / 86400),
               bp_med=float(np.median(bp)), h2_over_h1_med=float(np.median(1 / coh2)),
               h2_over_h1_lowf=float(1 / coh2[0]), h2_over_h1_highf=float(1 / coh2[-1]))
    curves = dict(f=f.tolist(), eta=eta.tolist(), bp=bp.tolist(), h2_over_h1=(1 / coh2).tolist())
    dev, _ = boot(np.arange(n))
    cov = np.cov(dev, rowvar=False)
    sd = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sd, sd)
    out.update(iid_med=float(np.median(sd)), iid_over_bp=float(np.median(sd / bp)),
               iid_corr_adjacent=float(np.mean(np.diag(corr, 1))) if len(f) > 1 else np.nan,
               iid_corr_lag2=float(np.mean(np.diag(corr, 2))) if len(f) > 2 else np.nan,
               iid_corr_ends=float(corr[0, -1]))
    curves["iid"] = sd.tolist()
    x = np.log(f / np.sqrt(f[0] * f[-1]))
    for d in BLOCK_DAYS:
        groups = np.floor((t - t[0]) / (86400 * d)).astype(int)
        n_groups = len(np.unique(groups))
        if n_groups < 6:
            for key in ("med", "over_iid", "corr_adjacent", "corr_ends", "common_share", "tilt_share",
                        "common_sd", "tilt_sd", "rest_sd", "n_blocks"):
                out[f"b{d}_{key}"] = np.nan
            continue
        dev, _ = boot(groups)
        cov = np.cov(dev, rowvar=False)
        sd_b = np.sqrt(np.diag(cov))
        corr = cov / np.outer(sd_b, sd_b)
        # decomposition of each bootstrap deviation: common level + tilt in log f + rest
        basis = np.column_stack([np.ones_like(x), x / (np.ptp(x) if np.ptp(x) > 0 else 1)])
        coef, *_ = np.linalg.lstsq(basis, dev.T, rcond=None)
        rest = dev - (basis @ coef).T
        total = np.sum(dev**2)
        out.update({
            f"b{d}_med": float(np.median(sd_b)), f"b{d}_over_iid": float(np.median(sd_b / sd)),
            f"b{d}_corr_adjacent": float(np.mean(np.diag(corr, 1))) if len(f) > 1 else np.nan,
            f"b{d}_corr_ends": float(corr[0, -1]),
            f"b{d}_common_share": float(np.sum((np.outer(coef[0], basis[:, 0])) ** 2) / total),
            f"b{d}_tilt_share": float(np.sum((np.outer(coef[1], basis[:, 1])) ** 2) / total),
            f"b{d}_common_sd": float(np.std(coef[0])), f"b{d}_tilt_sd": float(np.std(coef[1])),
            f"b{d}_rest_sd": float(np.sqrt(np.mean(rest**2))), f"b{d}_n_blocks": n_groups,
        })
        curves[f"b{d}"] = sd_b.tolist()
    # first half against second half of the record (non-stationarity over the deployment)
    half = t < np.median(t)
    a = _eta(f, k, attraction, zp[half].sum(0), pp[half].sum(0))
    b = _eta(f, k, attraction, zp[~half].sum(0), pp[~half].sum(0))
    d_half = np.log(a / b)
    out.update(half_level=float(np.mean(d_half)), half_tilt=float(np.polyfit(x, d_half, 1)[0] * np.ptp(x))
               if len(f) > 2 else np.nan, half_rms=float(np.sqrt(np.mean(d_half**2))))
    curves["half_log_ratio"] = d_half.tolist()
    # Dependence on the pressure level. With noise of constant power N on the pressure
    # channel, H1 of windows with pressure power S is H S/(S + N): eta rises with S.
    # Noise on the vertical leaves H1 unbiased at every S. Windows are grouped in five
    # classes of band-mean pressure power; eta of each class is fitted with
    # ln eta = a - b / S_rel, so b is N/S at the median level.
    level = pp.mean(1)
    edges = np.quantile(level, np.linspace(0, 1, 6))
    cls = np.clip(np.searchsorted(edges, level, side="right") - 1, 0, 4)
    eta_c, s_rel = [], []
    for c in range(5):
        m = cls == c
        eta_c.append(_eta(f, k, attraction, zp[m].sum(0), pp[m].sum(0)))
        s_rel.append(pp[m].mean() / np.median(level))
    eta_c, s_rel = np.array(eta_c), np.array(s_rel)
    slope = np.array([np.polyfit(1 / s_rel, np.log(eta_c[:, j]), 1)[0] for j in range(len(f))])
    out.update(level_ratio_top_bottom=float(np.median(eta_c[4] / eta_c[0])),
               level_ratio_lowf=float(eta_c[4, 0] / eta_c[0, 0]),
               level_ratio_highf=float(eta_c[4, -1] / eta_c[0, -1]),
               level_noise_fraction=float(np.median(-slope)),
               level_span_db=float(10 * np.log10(s_rel[4] / s_rel[0])))
    curves["level_classes"] = dict(s_rel=s_rel.tolist(), eta=eta_c.tolist())
    return dict(row=out, curves=curves)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out")
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--stations", nargs="*")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    jobs = []
    with open(PLAN) as fh:
        for r in csv.DictReader(fh):
            if args.stations and r["sid"] not in args.stations:
                continue
            jobs.append((r["network"], r["station"], r["wide_gain_used"], r["epoch_end"] or None))
    with mp.get_context("spawn").Pool(args.workers) as pool:
        results = pool.map(station, jobs, chunksize=1)
    rows = [r["row"] for r in results if "row" in r]
    failed = [r for r in results if "row" not in r]
    keys = list(rows[0])
    with open(out / "error_model.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    (out / "error_model_curves.json").write_text(
        json.dumps({r["row"]["sid"]: r["curves"] for r in results if "row" in r}))
    (out / "error_model_failed.json").write_text(json.dumps(failed, indent=1))
    print(f"{len(rows)} stations, {len(failed)} failed", flush=True)


if __name__ == "__main__":
    sys.exit(main())
