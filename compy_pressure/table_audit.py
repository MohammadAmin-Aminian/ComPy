"""Anomaly and correlation scan of the pressure QC/calibration master table.

ObsPy-free. Reads the one-row-per-pressure-epoch table that BRUIT-FM writes
(``pressure_qc_calibration_table.csv``) and, optionally, the per-channel
metadata offsets of the provisional decisions, and flags with
``compy.anomaly.Anomaly`` records:

* ``gain_outlier``: a channel whose gain differs from the other channels of
  its network and sensor type (robust z of 20 log10 g);
* ``method_disagreement``: Rayleigh, microseism, local-P and peer-offset
  multipliers of one channel that differ by >= 6 dB;
* ``polarity``: pressure signs that differ between methods on one channel,
  or within a network for one method;
* ``network_gain_structure``: a gain trend with depth across the table that
  vanishes within network x sensor groups (confounding, not physics);
* ``response_rolloff``: hydrophone gains that follow the low edge of the
  usable band (the calibration sits on a mis-modelled roll-off);
* ``infragravity_outlier`` / ``infragravity_basin``: Godin residuals far
  from their basin and depth peers, and basin-wide excesses;
* ``nearby_hidden``: channels that pass the broad expected-level limits but
  differ from their neighbours by >= 6 dB (a relative gain error);
* ``offset_outlier``: a channel offset far from its own network's offset;
* ``offset_invalid``: offsets that cannot be a measurement (dead channels
  give log10 of ~0 power, e.g. -3088 dB), which the peer-offset multiplier
  then turns into 1e154.

:func:`scan_table` runs them; ``compy.anomaly.rank`` orders the result.
"""

import csv

import numpy as np
from scipy import stats

from compy.anomaly import Anomaly, mad_sigma, p_to_sigma

FLOAT_COLUMNS = (
    "Lat", "Lon", "Profondeur", "mean_score", "score_expected_levels", "score_nearby",
    "score_godin", "pressure_expected_median_in_range_fraction", "nearby_median_difference_db",
    "infragravity_median_residual_db", "pressure_expected_band_min_hz",
    "pressure_expected_band_max_hz", "rayleigh_g", "rayleigh_events", "rayleigh_sign",
    "microseism_g", "microseism_coherence_squared", "local_p_g", "local_p_polarity",
    "score_rayleigh", "score_microseism", "score_local_p", "mean_score_all",
)  # fmt: skip


#: Larger peer offsets are not measurements (dead or constant data).
MAX_OFFSET_DB = 120.0


def _float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def load_table(path, offsets=None):
    """Columns of the master table as arrays (floats where numeric).

    Args:
        path (path): ``pressure_qc_calibration_table.csv``.
        offsets (path or None): ``metadata_offsets_channels.csv``; adds
            ``offset_db`` and ``peer_multiplier`` matched on Seed_id and start.

    Returns:
        dict: ``{column: ndarray}`` plus ``network`` and ``station``.
    """
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    table = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    for k in FLOAT_COLUMNS:  # absent columns read as not measured
        table[k] = (
            np.array([_float(v) for v in table[k]])
            if k in table
            else np.full(len(rows), np.nan)
        )
    for k in ("type", "basin", "depth_bin", "final_class", "channel_start_time"):
        table.setdefault(k, np.full(len(rows), "", dtype=object))
    table["network"] = np.array([s.split(".")[0] for s in table["Seed_id"]])
    table["station"] = np.array([".".join(s.split(".")[:2]) for s in table["Seed_id"]])
    table["offset_db"] = np.full(len(rows), np.nan)
    table["peer_multiplier"] = np.full(len(rows), np.nan)
    if offsets is not None:
        with open(offsets, newline="") as fh:
            index = {
                (r["Seed_id"], r["channel_start_time"]): r for r in csv.DictReader(fh)
            }
        for i, (sid, start) in enumerate(
            zip(table["Seed_id"], table["channel_start_time"])
        ):
            r = index.get((sid, start))
            if r:
                table["offset_db"][i] = _float(r["offset_db"])
                table["peer_multiplier"][i] = _float(r["multiplier_to_match_peers"])
    return table


def gains(table, *, min_coherence_squared=0.5):
    """``{method: array of g}`` with unusable estimates as NaN.

    Microseism estimates below ``min_coherence_squared`` are dropped, as in
    the provisional scoring.
    """
    micro = np.where(table["microseism_coherence_squared"] >= min_coherence_squared,
                     table["microseism_g"], np.nan)  # fmt: skip
    peer = np.where(
        np.abs(table["offset_db"]) <= MAX_OFFSET_DB, table["peer_multiplier"], np.nan
    )
    out = {"rayleigh": table["rayleigh_g"], "microseism": micro, "local_p": table["local_p_g"],
           "peer_offset": peer}  # fmt: skip
    return {k: np.where(np.isfinite(v) & (v > 0), v, np.nan) for k, v in out.items()}


def _db(g):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 20 * np.log10(g)


def gain_outliers(table, g, *, min_group=4, z_min=4.0, db_min=3.0):
    out = []
    groups = {}
    for method, values in g.items():
        if method == "peer_offset":
            continue
        for i in np.flatnonzero(np.isfinite(values)):
            groups.setdefault(
                (method, table["network"][i], table["type"][i]), []
            ).append(i)
    for (method, net, kind), idx in sorted(groups.items()):
        if len(idx) < min_group:
            continue
        db = _db(g[method][idx])
        med, scale = np.median(db), max(mad_sigma(db), 1.0)
        for i, v in zip(idx, db):
            z = (v - med) / scale
            if abs(z) >= z_min and abs(v - med) >= db_min:
                out.append(Anomaly(table["Seed_id"][i], "instrument", "gain_outlier",
                                   f"{method} g = {g[method][i]:.3g} is {v - med:+.1f} dB from the median of "
                                   f"{len(idx)} {net} {kind} channels ({10 ** (med / 20):.3g})",
                                   min(abs(z), 10.0), 100 * (10 ** ((v - med) / 20) - 1),
                                   table["channel_start_time"][i][:10],
                                   {"method": method, "group_median_g": float(10 ** (med / 20))}))  # fmt: skip
    return out


def method_disagreement(table, g, *, db_min=6.0):
    out = []
    pairs = [("rayleigh", "microseism"), ("rayleigh", "local_p"), ("microseism", "local_p"),
             ("rayleigh", "peer_offset"), ("microseism", "peer_offset")]  # fmt: skip
    for a, b in pairs:
        both = np.isfinite(g[a]) & np.isfinite(g[b])
        for i in np.flatnonzero(both):
            d = _db(g[a][i] / g[b][i])
            if abs(d) >= db_min:
                out.append(Anomaly(table["Seed_id"][i], "instrument", "method_disagreement",
                                   f"{a} g = {g[a][i]:.3g} vs {b} g = {g[b][i]:.3g} ({d:+.1f} dB)",
                                   min(abs(d) / 1.5, 10.0), 100 * (g[a][i] / g[b][i] - 1),
                                   table["channel_start_time"][i][:10], {"pair": f"{a}/{b}"}))  # fmt: skip
    return out


def polarity(table):
    out, info = [], {}
    # Signs only count where the method produced an estimate.
    rs = np.where(np.isfinite(table["rayleigh_g"]), table["rayleigh_sign"], np.nan)
    lp = np.where(np.isfinite(table["local_p_g"]), table["local_p_polarity"], np.nan)
    both = np.isfinite(rs) & np.isfinite(lp)
    for i in np.flatnonzero(both & (np.sign(rs) != np.sign(lp))):
        out.append(Anomaly(table["Seed_id"][i], "instrument", "polarity",
                           f"Rayleigh sign {rs[i]:+.0f} but local-P polarity {lp[i]:+.0f} on the same channel",
                           3.0, float("nan"), table["channel_start_time"][i][:10]))  # fmt: skip
    for name, values in (("Rayleigh sign", rs), ("local-P polarity", lp)):
        for net in np.unique(table["network"]):
            sel = (table["network"] == net) & np.isfinite(values)
            if sel.sum() < 3:
                continue
            pos, neg = int(np.sum(values[sel] > 0)), int(np.sum(values[sel] < 0))
            info.setdefault(name, {})[net] = (pos, neg)
            if pos and neg:
                minority = min(pos, neg)
                out.append(Anomaly(net, "instrument", "polarity",
                                   f"{name} mixed within the network: {pos} positive, {neg} negative",
                                   p_to_sigma(stats.binomtest(minority, pos + neg, 0.05, "greater").pvalue),
                                   float("nan")))  # fmt: skip
    return out, info


def depth_structure(table, g, *, min_group=4):
    """Depth trend of log gains across the table and within network x type groups."""
    out, info = [], {}
    depth = table["Profondeur"]
    for method in ("rayleigh", "microseism"):
        ok = np.isfinite(g[method]) & np.isfinite(depth)
        across = stats.spearmanr(depth[ok], _db(g[method][ok]))
        xs, ys = [], []
        for net in np.unique(table["network"][ok]):
            for kind in np.unique(table["type"][ok]):
                sel = ok & (table["network"] == net) & (table["type"] == kind)
                if sel.sum() >= min_group:
                    xs += list(depth[sel] - depth[sel].mean())
                    ys += list(_db(g[method][sel]) - _db(g[method][sel]).mean())
        within = stats.spearmanr(xs, ys) if len(xs) >= 8 else None
        slope = (
            float(np.polyfit(np.array(xs) / 1000, ys, 1)[0]) if len(xs) >= 8 else np.nan
        )
        info[method] = dict(n=int(ok.sum()), across_rho=float(across.statistic),
                            across_p=float(across.pvalue), n_within=len(xs),
                            within_rho=None if within is None else float(within.statistic),
                            within_p=None if within is None else float(within.pvalue),
                            within_slope_db_per_km=slope)  # fmt: skip
        if across.pvalue < 0.01 and within is not None and within.pvalue > 0.05:
            out.append(Anomaly("all", "network", "network_gain_structure",
                               f"{method} gains fall with depth across the table (Spearman {across.statistic:+.2f}, "
                               f"p = {across.pvalue:.1e}) but not within network × sensor groups "
                               f"(Spearman {within.statistic:+.2f}, p = {within.pvalue:.2f}, "
                               f"{slope:+.1f} dB/km): gains are network properties, not a depth effect",
                               p_to_sigma(across.pvalue), float("nan"), evidence=info[method]))  # fmt: skip
    return out, info


def response_rolloff(table, g):
    """Do hydrophone gains follow the low edge of the usable band?"""
    out, info = [], {}
    low = table["pressure_expected_band_min_hz"]
    for method in ("rayleigh", "microseism"):
        sel = (
            (table["type"] == "Hydrophone")
            & np.isfinite(g[method])
            & np.isfinite(low)
            & (low > 0)
        )
        if sel.sum() < 8:
            continue
        rho, p = stats.spearmanr(np.log(low[sel]), _db(g[method][sel]))
        info[method] = dict(n=int(sel.sum()), spearman=float(rho), p=float(p),
                            median_g=float(np.median(g[method][sel])),
                            median_band_min_hz=float(np.median(low[sel])))  # fmt: skip
        if p < 0.01:
            out.append(Anomaly("Hydrophone", "instrument", "response_rolloff",
                               f"hydrophone {method} g correlates with the low edge of the usable band "
                               f"(Spearman {rho:+.2f}, p = {p:.1e}, n = {sel.sum()}; median g = "
                               f"{np.median(g[method][sel]):.2f}): the calibration band sits on a roll-off "
                               "that the response does not describe",
                               p_to_sigma(p), 100 * (np.median(g[method][sel]) - 1), evidence=info[method]))  # fmt: skip
    return out, info


def infragravity(table, *, z_min=4.0, db_min=10.0):
    out, info = [], {}
    res = table["infragravity_median_residual_db"]
    ok = np.isfinite(res)
    for basin in np.unique(table["basin"][ok]):
        sel = ok & (table["basin"] == basin)
        info[basin] = dict(n=int(sel.sum()), median_db=float(np.median(res[sel])),
                           p10=float(np.percentile(res[sel], 10)), p90=float(np.percentile(res[sel], 90)))  # fmt: skip
        if sel.sum() >= 10 and abs(np.median(res[sel])) >= 6:
            out.append(Anomaly(basin, "network", "infragravity_basin",
                               f"{basin} infragravity levels {np.median(res[sel]):+.1f} dB from Godin et al. (2013) "
                               f"Eq. 6 (median of {sel.sum()} channels)",
                               min(10.0, abs(np.median(res[sel])) / 3), float("nan"), evidence=info[basin]))  # fmt: skip
        for dbin in np.unique(table["depth_bin"][sel]):
            grp = sel & (table["depth_bin"] == dbin)
            if grp.sum() < 6:
                continue
            med, scale = np.median(res[grp]), max(mad_sigma(res[grp]), 2.0)
            for i in np.flatnonzero(grp):
                z = (res[i] - med) / scale
                if abs(z) >= z_min and abs(res[i] - med) >= db_min:
                    out.append(Anomaly(table["Seed_id"][i], "noise", "infragravity_outlier",
                                       f"Godin residual {res[i]:+.0f} dB vs {med:+.0f} dB for {grp.sum()} {basin} "
                                       f"channels at {dbin}",
                                       min(abs(z), 10.0), float("nan"), table["channel_start_time"][i][:10]))  # fmt: skip
    return out, info


def nearby_hidden(table, *, db_min=6.0):
    out = []
    d, s1 = table["nearby_median_difference_db"], table["score_expected_levels"]
    for i in np.flatnonzero(
        np.isfinite(d) & np.isfinite(s1) & (s1 >= 1.0) & (np.abs(d) >= db_min)
    ):
        out.append(Anomaly(table["Seed_id"][i], "instrument", "nearby_hidden",
                           f"passes the expected-level limits but sits {d[i]:+.1f} dB from its neighbours",
                           min(abs(d[i]) / 2, 10.0), 100 * (10 ** (d[i] / 20) - 1),
                           table["channel_start_time"][i][:10]))  # fmt: skip
    return out


def offset_outliers(table, *, db_min=10.0, min_network=4):
    out = []
    off = np.where(
        np.abs(table["offset_db"]) <= MAX_OFFSET_DB, table["offset_db"], np.nan
    )
    for net in np.unique(table["network"]):
        sel = (table["network"] == net) & np.isfinite(off)
        if sel.sum() < min_network:
            continue
        med = np.median(off[sel])
        scale = max(mad_sigma(off[sel]), 3.0)
        for i in np.flatnonzero(sel):
            if abs(off[i] - med) >= db_min:
                out.append(Anomaly(table["Seed_id"][i], "instrument", "offset_outlier",
                                   f"level {off[i]:+.0f} dB from its basin's good gauges, {off[i] - med:+.0f} dB from "
                                   f"its network ({net} median {med:+.0f} dB)",
                                   min(abs(off[i] - med) / scale, 10.0), float("nan"),
                                   table["channel_start_time"][i][:10]))  # fmt: skip
    return out


def offset_invalid(table):
    out = []
    off = table["offset_db"]
    for i in np.flatnonzero(np.isfinite(off) & (np.abs(off) > MAX_OFFSET_DB)):
        out.append(Anomaly(table["Seed_id"][i], "instrument", "offset_invalid",
                           f"peer offset {off[i]:+.0f} dB stored as measured (class {table['final_class'][i]}); "
                           f"multiplier {table['peer_multiplier'][i]:.2g}",
                           5.0, float("nan"), table["channel_start_time"][i][:10]))  # fmt: skip
    return out


def scan_table(table):
    """All detectors; returns (anomalies, info)."""
    g = gains(table)
    found, info = [], {}
    found += gain_outliers(table, g)
    found += method_disagreement(table, g)
    a, info["polarity"] = polarity(table)
    found += a
    a, info["depth"] = depth_structure(table, g)
    found += a
    a, info["rolloff"] = response_rolloff(table, g)
    found += a
    a, info["infragravity"] = infragravity(table)
    found += a
    found += nearby_hidden(table)
    found += offset_outliers(table)
    found += offset_invalid(table)
    info["n_channels"] = len(table["Seed_id"])
    info["n_with_gain"] = {k: int(np.isfinite(v).sum()) for k, v in g.items()}
    return found, info
