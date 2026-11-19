"""Station-balanced noise models of groups of pressure gauges.

Each gauge counts once: at each period its histogram of window levels (after its
gain) is normalised to a probability, weighted by the gauge's usable-band taper,
and summed over the gauges of the group. Percentiles of that sum give the
group's low (5 %), median and high (95 %) curves; a Peterson-style fit of
continuous straight lines in log period, with the fewest segments reaching
1 dB RMS, gives the model. Curves are drawn only where at least 8 gauges count
at full weight (3 for thin groups, marked as such).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from . import native
from .classify import season_months
from .config import db_edges

GRID = native.GRID


def band_taper(f_lo, f_hi, periods=GRID):
    """Weight 1 inside the usable band [f_lo, f_hi] Hz, falling as cos^2 over the outer octave, 0 outside.

    ``f_lo`` of 0 (or less) means no lower limit, as for an absolute gauge coupled
    down to DC: no taper is applied on the long-period side.
    """
    f = 1.0 / np.asarray(periods, float)
    u = np.ones(len(f))
    if not (np.isfinite(f_lo) and np.isfinite(f_hi) and f_hi > 0):
        return u * 0
    hi = np.log2(f_hi / f)
    lo = np.log2(f / f_lo) if f_lo > 0 else np.full(len(f), np.inf)
    u[lo < 0] = 0
    u[hi < 0] = 0
    a = (lo >= 0) & (lo < 1)
    u[a] *= np.sin(0.5 * np.pi * lo[a]) ** 2
    b = (hi >= 0) & (hi < 1)
    u[b] *= np.sin(0.5 * np.pi * hi[b]) ** 2
    return u


def shift_hist(h, db, step=1.0):
    """Shift histograms along their last (dB) axis by ``db`` decibels, keeping the counts."""
    db = db / step
    if db == 0 or not np.isfinite(db):
        return h
    k = int(np.floor(db))
    a = db - k
    out = np.zeros_like(h)
    n = h.shape[-1]
    for kk, wt in ((k, 1 - a), (k + 1, a)):
        if 0 <= kk < n:
            out[..., kk:] += wt * h[..., : n - kk]
        elif -n < kk < 0:
            out[..., : n + kk] += wt * h[..., -kk:]
    return out


def percentiles(P, cover, edges, pct):
    """Percentiles (linear inside bins) and smoothed mode of a (period, dB) probability array."""
    centres = 0.5 * (edges[:-1] + edges[1:])
    out = {p: np.full(P.shape[0], np.nan) for p in pct}
    mode = np.full(P.shape[0], np.nan)
    for j in range(P.shape[0]):
        if cover[j] <= 0:
            continue
        c = np.cumsum(P[j])
        if c[-1] <= 0:
            continue
        c = c / c[-1]
        for p in pct:
            q = p / 100.0
            i = min(int(np.searchsorted(c, q)), len(c) - 1)
            prev = c[i - 1] if i > 0 else 0.0
            out[p][j] = edges[i] + (q - prev) / max(c[i] - prev, 1e-12) * (
                edges[i + 1] - edges[i]
            )
        mode[j] = centres[
            int(np.argmax(np.convolve(P[j], np.ones(3) / 3, mode="same")))
        ]
    return out, mode


def _continuous_fit(x, v, knots):
    """Least-squares levels at the knots of a continuous piecewise-linear curve; returns knots x, levels, RMS."""
    kx = x[knots]
    A = np.zeros((len(x), len(knots)))
    for c_ in range(len(knots)):
        e = np.zeros(len(knots))
        e[c_] = 1
        A[:, c_] = np.interp(x, kx, e)
    lv, *_ = np.linalg.lstsq(A, v, rcond=None)
    return kx, lv, float(np.sqrt(np.mean((A @ lv - v) ** 2)))


def fit_piecewise(T, y, max_rms=1.0, max_seg=20):
    """Continuous straight lines in log10 period with the fewest segments reaching ``max_rms`` dB.

    Breakpoints by dynamic programming over the sample points, levels at the
    breakpoints by least squares. Returns the knots, the RMS and maximum misfit,
    and the Peterson table: P(T) = A + B log10(T) on each segment. None if fewer
    than six finite points.
    """
    ok = np.isfinite(y)
    x, v = np.log10(np.asarray(T, float)[ok]), np.asarray(y, float)[ok]
    n = len(x)
    if n < 6:
        return None
    cost = np.full((n, n), np.inf)
    for i in range(n - 1):
        sx = sy = sxx = sxy = syy = 0.0
        for j in range(i, n):
            sx += x[j]
            sy += v[j]
            sxx += x[j] ** 2
            sxy += x[j] * v[j]
            syy += v[j] ** 2
            m = j - i + 1
            if m >= 2:
                den = m * sxx - sx**2
                b = (m * sxy - sx * sy) / den if den > 0 else 0
                a = (sy - b * sx) / m
                cost[i, j] = max(syy - a * sy - b * sxy, 0)
    D = np.full((max_seg + 1, n), np.inf)
    B = np.zeros((max_seg + 1, n), dtype=int)
    D[1] = cost[0]
    best = None
    for k in range(1, max_seg + 1):
        if k > 1:
            for j in range(1, n):
                c = D[k - 1][:j] + cost[1 : j + 1, j]
                s = int(np.argmin(c))
                D[k][j], B[k][j] = c[s], s
        if not np.isfinite(D[k][n - 1]):
            break  # more segments than the points can carry: keep the best fit found so far
        idx, j = [n - 1], n - 1
        for q in range(k, 1, -1):
            j = B[q][j]
            idx.append(j)
        knots = sorted(set([0] + idx))
        knots = [0] + [q for q in knots[1:-1] if q not in (0, n - 1)] + [n - 1]
        # The dynamic programme gives disjoint segments; a breakpoint between samples s and s+1 is placed on
        # whichever of the two gives the better continuous fit.
        kx, lv, rms = _continuous_fit(x, v, knots)
        for q in range(1, len(knots) - 1):
            if knots[q] + 1 < knots[q + 1]:
                trial = knots[:q] + [knots[q] + 1] + knots[q + 1 :]
                tx, tl, tr = _continuous_fit(x, v, trial)
                if tr < rms:
                    knots, kx, lv, rms = trial, tx, tl, tr
        r = np.interp(x, kx, lv) - v
        if best is None or rms < best["rms_db"]:
            best = dict(period_knots_s=[float(10**q) for q in kx], level_knots_db=[float(q) for q in lv],
                        rms_db=rms, max_abs_db=float(np.abs(r).max()), segments=len(knots) - 1)  # fmt: skip
        if rms <= max_rms:
            break
    per, lev = np.array(best["period_knots_s"]), np.array(best["level_knots_db"])
    Bs = np.diff(lev) / np.diff(np.log10(per))
    As = lev[:-1] - Bs * np.log10(per[:-1])
    best["table"] = [dict(period_from_s=float(per[i]), period_to_s=float(per[i + 1]), A_db=float(As[i]),
                          B_db_per_decade=float(Bs[i])) for i in range(len(Bs))]  # fmt: skip
    return best


def slug(s):
    return re.sub(
        r"[^A-Za-z0-9]+", "_", str(s).replace("<", "lt").replace(">", "gt")
    ).strip("_")


class Atlas:
    """Groups of gauges and their models.

    Args:
        table: one row per gauge, in the order of ``cube``; columns ``key``, ``net``, ``lat``, ``gain_db``.
        cube: array (gauge, month, period, dB) of window counts, gains not applied.
        taper: array (gauge, period) of usable-band weights (0 removes a gauge or a band).
        cfg: atlas settings.
        out_dir: models are written to ``out_dir/<axis>/<group>.{csv,json}`` and ``_density.npz``.
    """

    def __init__(self, table, cube, taper, cfg, out_dir=None):
        self.d = table.reset_index(drop=True)
        self.m = cube
        self.U = np.asarray(taper, float)
        self.cfg = cfg
        self.edges = db_edges(cfg)
        self.centres = 0.5 * (self.edges[:-1] + self.edges[1:])
        self.out = Path(out_dir) if out_dir else None
        self.index = []
        self.densities = {}

    def gauge_probability(self, u, months=None, gain=True):
        """(P, total windows, median) of one gauge: probability over dB at each period, after its gain."""
        h = np.asarray(self.m[u])
        h = h[months].sum(axis=0) if months is not None else h.sum(axis=0)
        if gain:
            h = shift_hist(h, float(self.d.gain_db[u]), self.edges[1] - self.edges[0])
        tot = h.sum(axis=1)
        ok = tot >= self.cfg["windows"]["min_windows_per_period"]
        P = np.zeros_like(h)
        P[ok] = h[ok] / tot[ok, None]
        c = np.cumsum(P, axis=1)
        med = np.full(h.shape[0], np.nan)
        med[ok] = self.centres[np.argmax(c[ok] >= 0.5, axis=1)]
        return P, tot, ok, med

    def group(self, rows, months_of=None, gain=True):
        nP = self.m.shape[2]
        Pg, cover, used = (
            np.zeros((nP, len(self.centres))),
            np.zeros(nP),
            np.zeros(nP, int),
        )
        nets, meds, windows = [set() for _ in range(nP)], [], 0.0
        for u in rows:
            ms = months_of(u) if months_of else None
            if ms is not None and not len(ms):
                continue
            P, tot, ok, med = self.gauge_probability(u, ms, gain)
            w = self.U[u] * ok
            if not w.any():
                continue
            Pg += w[:, None] * P
            cover += w
            full = w > 0.5
            used += full
            for j in np.where(full)[0]:
                nets[j].add(self.d.net[u])
            med[~full] = np.nan
            meds.append(med)
            windows += float(tot.max())
        return (
            Pg,
            cover,
            np.array(meds),
            windows,
            used,
            np.array([len(s) for s in nets]),
        )

    def add(self, axis, name, rows, months_of=None, gain=True, save=True):
        """Model of one group; returns the per-period table and an index record."""
        mc = self.cfg["models"]
        rows = [u for u in rows if self.U[u].any()]
        P, cover, meds, windows, used, nets = self.group(rows, months_of, gain)
        n_units = len(meds)
        need = (
            mc["min_stations"]
            if n_units >= mc["min_stations"]
            else mc["min_stations_thin"]
        )
        pc, mode = percentiles(
            P, np.where(used >= need, cover, 0), self.edges, mc["percentiles"]
        )
        with np.errstate(all="ignore"):
            st = {
                q: np.nanpercentile(meds, q, axis=0)
                if n_units
                else np.full(len(GRID), np.nan)
                for q in (10, 50, 90)
            }
        tab = pd.DataFrame({"period_s": GRID, "frequency_hz": 1 / GRID, **{f"p{p}_db": pc[p] for p in mc["percentiles"]},
                            "mode_db": mode, **{f"gauge_median_p{q}_db": st[q] for q in (10, 50, 90)},
                            "gauges": used, "networks": nets})  # fmt: skip
        tab.loc[used < need, [f"gauge_median_p{q}_db" for q in (10, 50, 90)]] = np.nan
        rec = dict(axis=axis, group=str(name), gauges=n_units, networks=int(self.d.net[list(rows)].nunique()) if rows else 0,
                   window_count=windows, status="model" if n_units >= mc["min_stations"] else "thin_curves_only")  # fmt: skip
        rec["single_experiment"] = rec["networks"] < mc["robust_networks"]
        out = dict(axis=axis, group=str(name), units="dB re 1 Pa^2/Hz",
                   percentile_basis="station-balanced distribution of 3-h windows",
                   **{k: rec[k] for k in ("gauges", "networks", "window_count", "status", "single_experiment")})  # fmt: skip
        if rec["status"] == "model":
            mask = (used >= max(mc["min_stations"], mc["fit_share"] * used.max())) & (
                (GRID >= 0.5) | (used >= 0.5 * used.max())
            )
            for lab, col in (
                ("low", f"p{mc['low_percentile']}_db"),
                ("median", "p50_db"),
                ("high", f"p{mc['high_percentile']}_db"),
            ):
                f = fit_piecewise(
                    GRID,
                    np.where(mask, tab[col].values, np.nan),
                    mc["fit_max_rms_db"],
                    mc["fit_max_segments"],
                )
                out[lab] = f
                if f is not None:
                    rec.update(
                        {
                            f"{lab}_segments": f["segments"],
                            f"{lab}_rms_db": f["rms_db"],
                            f"{lab}_max_db": f["max_abs_db"],
                        }
                    )
                    tab[f"{lab}_fit_db"] = np.interp(np.log10(GRID), np.log10(f["period_knots_s"]), f["level_knots_db"],
                                                     left=np.nan, right=np.nan)  # fmt: skip
        Pn = P / np.maximum(P.sum(1, keepdims=True), 1e-12)
        Pn[used < need] = np.nan
        self.densities[(axis, str(name))] = (Pn.astype(np.float32), used)
        if save and self.out is not None:
            dd = self.out / axis
            dd.mkdir(parents=True, exist_ok=True)
            tab.round(3).to_csv(dd / f"{slug(name)}.csv", index=False)
            (dd / f"{slug(name)}.json").write_text(json.dumps(out, indent=1))
            np.savez_compressed(dd / f"{slug(name)}_density.npz", density=Pn.astype(np.float32), periods=GRID,
                                db_edges=self.edges, gauges_full_weight=used)  # fmt: skip
        self.index.append(rec)
        return tab, rec

    def by(self, axis, column, rows=None):
        """One model per value of a table column (within ``rows``)."""
        d = self.d if rows is None else self.d.iloc[list(rows)]
        for g, s in d.groupby(column):
            self.add(axis, g, list(s.index))

    def season_paired(self, depth_groups=((0, 200), (200, 1000), (1000, 1e9))):
        """Local winter minus summer median per gauge, then quartiles over gauges by hemisphere and depth group."""
        d = self.d
        diffs = {}
        for u in range(len(d)):
            if not self.U[u].any():
                continue
            med = {
                s: self.gauge_probability(u, season_months(s, d.lat[u]))[3]
                for s in ("winter", "summer")
            }
            x = med["winter"] - med["summer"]
            x[self.U[u] < 0.5] = np.nan
            diffs[u] = x
        rows = []
        for hemi, hm in (("north", d.lat >= 0), ("south", d.lat < 0)):
            for lo, hi in ((0, 1e9),) + tuple(depth_groups):
                us = [u for u in diffs if hm[u] and lo <= d.depth_m[u] < hi]
                if not us:
                    continue
                X = np.vstack([diffs[u] for u in us])
                n = np.isfinite(X).sum(axis=0)
                with np.errstate(all="ignore"):
                    q = np.nanpercentile(X, [25, 50, 75], axis=0)
                q[:, n < self.cfg["models"]["min_stations"]] = np.nan
                label = "all depths" if hi >= 1e9 and lo == 0 else f"{lo:g}-{hi:g} m"
                rows.append(pd.DataFrame({"hemisphere": hemi, "depth_group": label, "period_s": GRID,
                                          "p25_db": q[0], "median_db": q[1], "p75_db": q[2], "gauges": n}))  # fmt: skip
        return pd.concat(rows, ignore_index=True)

    def index_table(self):
        return pd.DataFrame(self.index)
