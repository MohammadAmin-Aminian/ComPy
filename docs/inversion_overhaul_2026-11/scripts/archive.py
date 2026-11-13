"""Read-only access to the 2026 archive (~/Desktop/ComPy) for the overhaul scripts.

``station_inputs`` returns, for one station of the all-station plan, the stacked
compliance of the wide-gate stage C result with a covariance measured from the
accepted windows in the per-window spectra cache (13_audit/spectra_cache), and
the CRUST1.0 column at the station.
"""

import csv
import json
from pathlib import Path

import numpy as np

ARCHIVE = Path.home() / "Desktop/ComPy"
CACHE = ARCHIVE / "13_audit/spectra_cache"
STAGE_C = ARCHIVE / "06_compliance_wide"
PLAN = ARCHIVE / "12_inversion_all_stations/tables/station_plan.csv"
CRUST1 = Path.home() / "Desktop/BRUIT-FM/compliance_yv_2026-10-28/literature/crust1"
OUT = ARCHIVE / "14_inversion_overhaul"


def plan():
    """Rows of the all-station plan as dictionaries keyed by 'NET.STA'."""
    with open(PLAN) as fh:
        return {r["sid"]: r for r in csv.DictReader(fh)}


def station_windows(net, sta, label, epoch_end=None):
    """Accepted windows of the stage C result: band, f, pp, zp, starts, gain, depth, stored."""
    folder = STAGE_C / net / "results" / sta
    if epoch_end:
        folder = folder / "epoch"
    res = dict(np.load(folder / f"{sta}_compliance_{label}.npz"))
    cache = np.load(CACHE / net / f"{sta}.npz")
    f_all, starts = cache["f"], cache["starts"]
    lo, hi = (float(v) for v in res["band"])
    band = (f_all >= lo - 1e-12) & (f_all <= hi + 1e-12)
    order = np.argsort(starts)
    idx = order[np.searchsorted(starts[order], res["window_starts"])]
    if not np.allclose(starts[idx], res["window_starts"], atol=1e-3):
        raise RuntimeError("accepted windows not found in the cache")
    return dict(f=f_all[band], pp=cache["pp"][idx][:, band], zp=cache["zp"][idx][:, band],
                zz=cache["zz"][idx][:, band], starts=starts[idx], gain=float(res["gain"]),
                depth=float(res["depth"]), band=(lo, hi), stored=res,
                fc=float(cache["fc"]))  # fmt: skip


LEAK_TABLE = (
    Path(__file__).resolve().parents[2] / "audit_2026-11-data/data/dpg_leak_factors.csv"
)
PHASE_TABLE = ARCHIVE / "13_audit/phase/phase_metadata.csv"


def response_bias(sid, freqs):
    """Compliance bias [%] of the pressure-gauge leak corner measured from the tide.

    From the per-station leak factor of the 2.9 audit
    (docs/audit_2026-11-data/data/dpg_leak_factors.csv, column k_recommended) and
    the channel's StationXML stages, through ``compy.response_check.compliance_bias``.
    Returns (bias, record) or (None, None) for a station without a factor.
    """
    from obspy import read_inventory

    from compy import response_check as rc

    with open(LEAK_TABLE) as fh:
        leak = {r["sid"]: r for r in csv.DictReader(fh)}
    if sid not in leak or not leak[sid]["k_recommended"]:
        return None, None
    with open(PHASE_TABLE) as fh:
        phase = {r["sid"]: r for r in csv.DictReader(fh)}
    net, sta = sid.split(".")
    inv = read_inventory(
        str(ARCHIVE / "00_inventory/stationxml" / f"{net}.xml")
    ).select(station=sta)
    chans = [c for n in inv for s in n for c in s if c.code == phase[sid]["p_channel"]]
    if not chans:
        return None, None
    k = float(leak[sid]["k_recommended"])
    bias = rc.compliance_bias(rc.paz_stages(chans[0]), np.asarray(freqs, float), k)
    record = dict(
        k=k,
        source=leak[sid]["k_source"],
        channel=phase[sid]["p_channel"],
        bias_low_pct=float(bias[0]),
        bias_high_pct=float(bias[-1]),
    )
    # Uncertainty of the correction: half the difference between the tide and the
    # phase estimate of the bias at the band bottom, where both exist.
    try:
        record["tide_minus_phase_pct"] = float(leak[sid]["bias_lo_tide_pct"]) - float(
            leak[sid]["bias_lo_phase_pct"]
        )
    except ValueError:
        record["tide_minus_phase_pct"] = None
    return bias, record


def station_inputs(sid, *, max_fraction_fc=0.8, response="auto", seed=0, **data_kwargs):
    """ComplianceData, CRUST1.0 column and a record of the inputs for one station.

    ``response="auto"`` corrects the compliance for the pressure-gauge leak corner
    where the audit measured one (``response_bias``) and marks the response
    "corrected"; other stations stay "unchecked". "unchecked" forces no correction.
    """
    from compy.inversion import auto, dram
    from compy.inversion.errors import stack_covariance, usable_bins

    row = plan()[sid]
    net, sta = row["network"], row["station"]
    w = station_windows(net, sta, row["wide_gain_used"], row["epoch_end"] or None)
    stat = stack_covariance(w["f"], w["pp"], w["zp"], w["starts"], w["depth"], gain=w["gain"],
                            seed=seed)  # fmt: skip
    keep = usable_bins(stat["sd_rel"])
    trimmed = (
        [float(w["f"][~keep].min()), float(w["f"][~keep].max())]
        if not keep.all()
        else None
    )
    correction = None
    if response == "auto":
        bias, correction = response_bias(sid, w["f"])
        response = "unchecked" if bias is None else "corrected"
        if bias is not None:
            stat["eta"] = stat["eta"] / (1.0 + bias / 100.0)
    if correction is not None and "tilt_rel" not in data_kwargs:
        diff = correction["tide_minus_phase_pct"]
        spread = 0.0 if diff is None else 0.5 * abs(diff) / 100.0
        data_kwargs["tilt_rel"] = float(
            np.hypot(spread, auto.RESPONSE_TILT["corrected"])
        )
    data = auto.compliance_data(w["f"][keep], stat["eta"][keep], w["depth"],
                                cov_rel=stat["cov_rel"][np.ix_(keep, keep)],
                                max_fraction_fc=max_fraction_fc, response=response,
                                station=sid, **data_kwargs)  # fmt: skip
    column = dram.crust1_column(float(row["latitude"]), float(row["longitude"]), CRUST1)
    info = dict(sid=sid, gain_source=row["wide_gain_used"], gain=w["gain"], depth=w["depth"],
                n_windows=int(len(w["starts"])), stage_c_band=list(w["band"]),
                trimmed_noisy_bins=trimmed,
                band=[float(data.freqs[0]), float(data.freqs[-1])], n_data=int(len(data.freqs)),
                f_c=w["fc"], slow_level=stat["slow_level"], slow_tilt=stat["slow_tilt"],
                block_days=stat["block_days"], n_blocks=stat["n_blocks"],
                sd_rel_median=float(np.median(stat["sd_rel"])),
                half_level=stat["half_level"], half_tilt=stat["half_tilt"],
                latitude=float(row["latitude"]), longitude=float(row["longitude"]),
                tier=int(row["tier"]), response=response, response_correction=correction,
                tilt_rel=data_kwargs.get("tilt_rel", auto.RESPONSE_TILT[response]),
                crust1=json.loads(json.dumps(column, default=float)))  # fmt: skip
    return data, column, info
