"""Per-channel products of the atlas from the native-rate spectra (``native``).

For one channel: the window screens (no level, clipped, dead, deployment
edges, dropped days), monthly histograms of the window levels on the dB grid
(one array per channel: month x period x dB), and daily band levels used by the
day and band screens. Gains are not applied here; the models shift each
channel's histogram by its gain.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from . import native
from .config import db_edges

FULL_SCALE_24 = 2.0**23
FULL_SCALE_32 = 2.0**31


def histogram(levels, months, edges):
    """Counts (12, n_periods, n_db) of window levels; NaN levels and levels outside the edges are not counted."""
    n_w, n_p = levels.shape
    nb = len(edges) - 1
    out = np.zeros((12, n_p, nb), np.float32)
    if n_w == 0:
        return out
    b = np.floor((levels - edges[0]) / (edges[1] - edges[0]))
    ok = np.isfinite(b) & (b >= 0) & (b < nb)
    m = np.broadcast_to(months[:, None], levels.shape)
    p = np.broadcast_to(np.arange(n_p)[None, :], levels.shape)
    flat = (m[ok] * n_p + p[ok]) * nb + b[ok].astype(np.int64)
    np.add.at(out.reshape(-1), flat, 1.0)
    return out


def edge_days(days, level, edge_db, n_days):
    """Days at either end of the record whose level departs from the record median by more than ``edge_db``.

    Searched over the first and last ``n_days`` distinct days: the screen removes
    every day up to the last departing one, so a descent or a recovery is cut whole.
    """
    ud = np.unique(days)
    if not len(ud):
        return set()
    dm = np.array([np.nanmedian(level[days == x]) for x in ud])
    out = np.abs(dm - np.nanmedian(level)) > edge_db
    drop = set()
    for seq in (
        np.arange(min(n_days, len(ud))),
        np.arange(len(ud) - 1, max(len(ud) - 1 - n_days, -1), -1),
    ):
        hit = [k for k, q in enumerate(seq) if out[q]]
        if hit:
            drop.update(ud[seq[: hit[-1] + 1]].tolist())
    return drop


def screen_windows(z, cfg, drop_days=None):
    """Boolean masks of the window screens of one channel (all True = kept) and their counts."""
    lv = z["levels_db"]
    w = cfg["windows"]
    finite = np.isfinite(lv).any(axis=1)
    mx = np.asarray(z["maxabs_counts"], float)
    fs_24 = np.mean(mx >= w["clip_fraction"] * FULL_SCALE_24) <= 0.5
    clipped = mx >= w["clip_fraction"] * (FULL_SCALE_24 if fs_24 else FULL_SCALE_32)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lev = np.nanmedian(np.where(finite[:, None], lv, 0.0), axis=1)
    live = finite & ~clipped
    dead = np.zeros(len(lv), bool)
    if live.any():
        dead = live & (lev < np.median(lev[live]) - w["dead_window_db"])
    keep = live & ~dead
    days = np.floor(np.asarray(z["starts"], float) / 86400.0).astype(int)
    a, b = w["edge_band_s"]
    band = (native.GRID >= a) & (native.GRID <= b)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lev_band = np.nanmean(lv[:, band], axis=1)
    edge = np.zeros(len(lv), bool)
    if keep.any() and np.isfinite(lev_band[keep]).any():
        drop = edge_days(days[keep], lev_band[keep], w["edge_db"], w["edge_days"])
        edge = keep & np.isin(days, list(drop))
    keep &= ~edge
    dropped = np.zeros(len(lv), bool)
    if drop_days is not None and len(drop_days):
        dropped = keep & np.isin(days, np.asarray(list(drop_days)))
        keep &= ~dropped
    counts = dict(n_total=len(lv), n_no_level=int((~finite).sum()), n_clipped=int((finite & clipped).sum()),
                  n_dead=int(dead.sum()), n_edge=int(edge.sum()), n_day_screen=int(dropped.sum()),
                  n_kept=int(keep.sum()), full_scale_bits=24 if fs_24 else 32)  # fmt: skip
    return keep, counts


def months_of(starts):
    """Calendar month index (0 = January) of window start times [s since 1970]."""
    t = (np.asarray(starts, float) + 5400.0).astype("datetime64[s]")
    return (t.astype("datetime64[M]").astype(int) % 12).astype(int)


def daily_band_levels(z, keep, bands, min_windows=4):
    """Median over each day of the window band means [dB], kept windows only; days with fewer windows are NaN."""
    lv = z["levels_db"][keep]
    days = np.floor(np.asarray(z["starts"], float)[keep] / 86400.0).astype(int)
    cols = {}
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for name, (a, b) in bands.items():
            sel = (native.GRID >= a) & (native.GRID <= b)
            cols[name] = np.nanmean(lv[:, sel], axis=1)
    df = pd.DataFrame(cols)
    df["day"] = days
    g = df.groupby("day")
    out = g.median()
    out[g.size() < min_windows] = np.nan
    return out


def channel_products(path, cfg, drop_days=None):
    """Window screens, histograms (all windows with a level, and kept windows) and daily band levels.

    Args:
        path: native-rate product file of the channel.
        cfg: atlas settings.
        drop_days: days (since 1970) removed by the day screen.
    Returns:
        (counts dict, cube_all, cube_kept, daily DataFrame)
    """
    z = native.load(path)
    edges = db_edges(cfg)
    keep, counts = screen_windows(z, cfg, drop_days)
    months = months_of(z["starts"])
    finite = np.isfinite(z["levels_db"]).any(axis=1)
    cube_all = histogram(z["levels_db"][finite], months[finite], edges)
    cube_kept = histogram(z["levels_db"][keep], months[keep], edges)
    daily = daily_band_levels(z, keep, cfg["bands"])
    if keep.any():
        t = np.asarray(z["starts"], float)[keep]
        counts.update(t_start=float(t.min()), t_end=float(t.max() + native.WINDOW_S),
                      months_present=int(len(np.unique(months[keep]))))  # fmt: skip
    counts.update(fs=z["fs"], input_units=z["input_units"], sensitivity=float(z["sensitivity"][0]),
                  n_epochs=int(z["response_db"].shape[0]))  # fmt: skip
    return counts, cube_all, cube_kept, daily
