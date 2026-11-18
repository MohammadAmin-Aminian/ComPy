"""Screens that keep failing instruments out of the noise models.

Principle (atlas supplement, 9 November 2026): ocean noise is shared across a
network; a level that a gauge's network peers at similar depth do not see, or
that its own seismometer does not predict, is instrumental.

* ``Cascade``: channel-level selection with a reason for each removal.
* ``level_gate``: a gauge's level in a short-period band against its basin median.
* ``calibration_screen``: the P-Z calibration level against DPGs and APGs of the same depth class.
* ``day_screen``: day levels against same-network, same-depth peers on the same day.
* ``band_screen``: a gauge's band median against its peers, deep gauges only.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


class Cascade:
    """Ordered channel selection: each step removes the channels its mask rejects, with a reason."""

    def __init__(self, table):
        self.d = table.copy()
        self.d["drop_reason"] = ""
        self.steps = [dict(step="all channels", removed=0, remaining=len(self.d))]

    @property
    def live(self):
        return self.d.drop_reason.eq("")

    def step(self, keep_mask, reason):
        keep = pd.Series(keep_mask, index=self.d.index).fillna(False).astype(bool)
        out = self.live & ~keep
        self.d.loc[out, "drop_reason"] = reason
        self.steps.append(dict(step=reason, removed=int(out.sum()), remaining=int(self.live.sum())))
        return out

    def table(self):
        return pd.DataFrame(self.steps)


def level_gate(level_db, group, live, max_db):
    """Offset of each gauge's level from the median of the live gauges of its group; True where within ``max_db``.

    Gauges without a level pass (they are tested in the other gate band).
    """
    level_db = pd.Series(level_db, index=group.index)
    med = level_db[live].groupby(group[live]).median()
    ref = group.map(med).fillna(level_db[live].median())
    off = level_db - ref
    return off, ~(off.abs() > max_db)


def calibration_screen(table, pz_offset_db, cfg):
    """Gauges and network sensor types whose P-Z calibration level fails.

    Args:
        table: gauges with ``key``, ``net``, ``type``, ``depth_bin``, ``depth_m``.
        pz_offset_db: Series by key (see ``pz_level`` and ``pz_offsets``): the P-Z level
            R = 20 log10(|S_zp / S_zz| / w^2) - 20 log10(rho H) at 0.2 Hz, where a
            Rayleigh wave gives p = rho H a (An et al. 2017), minus the median R of the
            DPGs and APGs of the same depth class (gauges at ``deep_m`` or deeper).
    Returns:
        DataFrame by key with ``cal_offset_db``, ``cal_excluded``, ``cal_reason``.
    """
    s = cfg["screens"]
    d = table[["key", "net", "type"]].copy()
    d["cal_offset_db"] = d.key.map(pz_offset_db)
    fails = d.cal_offset_db.abs() > s["cal_db"]
    net = d[d.cal_offset_db.notna()].groupby(["net", "type"]).agg(n=("key", "size"), med=("cal_offset_db", "median"))
    bad = set(net[(net.n >= 3) & (net.med.abs() > s["cal_db"])].index)
    net_bad = pd.Series([(a, b) in bad for a, b in zip(d.net, d.type)], index=d.index)
    d["cal_excluded"] = net_bad | fails
    d["cal_reason"] = np.where(net_bad, "network fails P-Z level", np.where(fails, "gauge fails P-Z level", ""))
    return d.set_index("key")[["cal_offset_db", "cal_excluded", "cal_reason"]]


def pz_level(npz_path, depth_m, gain_db=0.0, f0=0.2, min_coherence=0.7, min_windows=20, rho=1028.0):
    """P-Z level R [dB] of one gauge from a ComPy stage B spectra file (``<NET>.<STA>_raw.npz``).

    Median over windows with squared coherence >= ``min_coherence`` at ``f0`` of
    20 log10(|S_zp / S_zz| / w^2) - 20 log10(rho H), with the pressure gain used
    in that file removed and ``gain_db`` applied. NaN with fewer than ``min_windows``.
    """
    import json

    z = np.load(npz_path)
    g = json.loads(str(z["settings"])).get("pressure_gain", 1.0)
    f = z["freqs"]
    j = int(np.argmin(abs(f - f0)))
    szo, szz, soo = z["ZP_szo"][:, j], z["ZP_szz"][:, j].astype(float), z["ZP_soo"][:, j].astype(float)
    with np.errstate(all="ignore"):
        sel = np.abs(szo) ** 2 / (szz * soo) >= min_coherence
    if sel.sum() < min_windows:
        return np.nan
    w = 2 * np.pi * f[j]
    return float(np.median(20 * np.log10(np.abs(szo[sel]) / szz[sel] / w**2)) - 20 * np.log10(rho * depth_m)
                 - 20 * np.log10(g) + gain_db)  # fmt: skip


def pz_offsets(levels, table):
    """P-Z calibration offset of each deep gauge: ``levels`` (Series by key of R at 0.2 Hz, dB) minus the median of
    the DPGs and APGs of its depth class (or of all DPGs and APGs when the class has none)."""
    d = table.set_index("key")
    R = pd.Series(levels).reindex(d.index)
    ctrl = d.type.isin(["DPG", "APG"])
    ref = R[ctrl].groupby(d.depth_bin[ctrl]).median()
    return R - d.depth_bin.map(ref).fillna(R[ctrl].median())


def usable_band(f_lo, f_hi, band):
    """True when a period band (s) lies inside a usable frequency band (Hz)."""
    lo, hi = band
    return (1.0 / hi >= f_lo) and (1.0 / lo <= f_hi)


def day_screen(daily, table, cfg):
    """Days to drop per gauge.

    A gauge's daily band levels minus the same-day median of the other gauges of
    its network and depth class (at least ``min_peers``), less the gauge's usual
    offset; a day departing by more than ``day_db`` in any usable band is
    dropped. Gauges at ``deep_m`` or deeper without enough peers are tested
    against their own median with ``self_db``. Shallower gauges without peers
    are not tested: storms raise their swell bands by tens of decibels.

    Args:
        daily: dict key -> DataFrame (index day, columns band names) of daily band levels.
        table: gauges with ``key``, ``net``, ``depth_bin``, ``depth_m``, ``f_lo``, ``f_hi``.
    Returns:
        (drops: dict key -> set of days, summary DataFrame)
    """
    s = cfg["screens"]
    bands = cfg["bands"]
    d = table.set_index("key")
    grp = d.net.astype(str) + "|" + d.depth_bin.astype(str)
    drops, rows = {}, []
    for _, keys in grp.groupby(grp):
        keys = [k for k in keys.index if daily.get(k) is not None and len(daily[k])]
        for k in keys:
            r, x = d.loc[k], daily[k]
            others = [daily[o] for o in keys if o != k]
            bad = pd.Series(False, index=x.index)
            tested = []
            for b, rng in bands.items():
                if b not in x or not usable_band(r.f_lo, r.f_hi, rng) or x[b].isna().all():
                    continue
                if len(others) >= s["min_peers"]:
                    M = pd.concat([o[b] for o in others if b in o], axis=1).reindex(x.index)
                    ok = M.notna().sum(axis=1) >= s["min_peers"]
                    res = x[b] - M.median(axis=1)
                    a = res - res[ok].median()
                    bad |= ok & (a.abs() > s["day_db"])
                    tested.append(f"{b}:peers")
                elif r.depth_m >= s["deep_m"]:
                    bad |= (x[b] - x[b].median()).abs() > s["self_db"]
                    tested.append(f"{b}:self")
            drops[k] = set(int(t) for t in x.index[bad.values])
            rows.append(dict(key=k, days=len(x), days_dropped=len(drops[k]),
                             share_dropped=len(drops[k]) / max(len(x), 1), tests=";".join(tested)))  # fmt: skip
    return drops, pd.DataFrame(rows)


def band_screen(daily, drops, table, cfg):
    """Gauge-bands whose median (after dropped days, with gain) departs from peers by more than ``band_db``.

    Deep gauges only (``deep_m``): shallower, levels at 10-50 s change steeply
    with depth inside a class. Peers: the other gauges of the same network and
    depth class (at least ``band_min_peers_network``), else the depth class in
    all networks (at least ``min_peers``).
    """
    s = cfg["screens"]
    d = table.set_index("key")
    med = {}
    for k in d.index:
        x = daily.get(k)
        if x is None or not len(x):
            continue
        med[k] = x[~x.index.isin(drops.get(k, set()))].median()
    M = pd.DataFrame(med).T.reindex(d.index).add(d.gain_db, axis=0)
    live = d[(~d.get("cal_excluded", pd.Series(False, index=d.index)).astype(bool)) & (d.depth_m >= s["deep_m"])]
    rows = []
    for k, r in live.iterrows():
        for b, rng in cfg["bands"].items():
            if b not in M or not usable_band(r.f_lo, r.f_hi, rng) or not np.isfinite(M.loc[k, b]):
                continue

            def ok(o, b=b, rng=rng):
                return o != k and np.isfinite(M.loc[o, b]) and usable_band(d.f_lo[o], d.f_hi[o], rng)

            peers = [o for o in live.index[(live.net == r.net) & (live.depth_bin == r.depth_bin)] if ok(o)]
            ref = "network and depth class"
            if len(peers) < s["band_min_peers_network"]:
                peers = [o for o in live.index[live.depth_bin == r.depth_bin] if ok(o)]
                ref = "depth class, all networks"
            if len(peers) < s["min_peers"]:
                continue
            off = M.loc[k, b] - M.loc[peers, b].median()
            rows.append(dict(key=k, band=b, offset_db=off, peers=len(peers), reference=ref, removed=abs(off) > s["band_db"]))
    return pd.DataFrame(rows, columns=["key", "band", "offset_db", "peers", "reference", "removed"])
