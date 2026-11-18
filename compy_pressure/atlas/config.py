"""Settings of the pressure atlas. Every threshold used by the atlas lives here.

``load(path)`` reads a JSON or YAML file and overrides the defaults key by key,
so a run file only names what it changes. The values are those of the atlas of
6 November 2026 and its supplement of 9 November 2026 (ComPy branch
``analysis/pressure-atlas``), unless a comment says otherwise.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

DEFAULTS = {
    "histogram": {"db_min": -150.0, "db_max": 130.0, "db_step": 1.0},
    "windows": {
        "dead_window_db": 30.0,  # windows this far below the channel median are dropouts
        "edge_db": 15.0,  # deployment-edge days: 1-10 s day level this far from the channel median
        "edge_days": 20,  # searched at each end of the record
        "edge_band_s": [1.0, 10.0],
        "clip_fraction": 0.99,  # of the digitiser full scale (24 bit; 32 bit when most windows exceed 24 bit)
        "min_windows_per_period": 28,  # 42 h of 3-h windows at 50 % overlap
    },
    "population": {
        "min_windows": 112,  # 7 days of 3-h windows at 50 % overlap
        "level_gate_db": 20.0,  # gauge level within this of its basin median, in each gate band
        "level_gate_bands_s": [[1.0, 2.0], [2.5, 4.0]],
        "one_hz_min_period_s": 4.0,  # 1-Hz channels: decimation filter not in StationXML (A3, 18 Nov 2026)
        "nyquist_fraction": 0.4,  # highest usable frequency as a fraction of the sampling rate
    },
    "screens": {
        "cal_db": 20.0,  # P-Z calibration level against same-depth DPG/APG at 0.2 Hz
        "day_db": 20.0,  # day level against network peers of the same depth class
        "self_db": 30.0,  # day level against the gauge's own median, deep gauges without peers
        "band_db": 25.0,  # gauge-band median against peers
        "min_peers": 3,
        "band_min_peers_network": 5,
        "deep_m": 1000.0,
    },
    "models": {
        "min_stations": 8,  # fewer: percentile curves only, drawn thin
        "min_stations_thin": 3,  # fewer: nothing drawn
        "robust_networks": 3,
        "low_percentile": 5,
        "high_percentile": 95,
        "percentiles": [1, 5, 10, 50, 90, 95, 99],
        "fit_max_rms_db": 1.0,
        "fit_max_segments": 20,
        "fit_share": 0.25,  # line fits only where the group keeps this share of its gauges (and min_stations)
    },
    "bands": {
        "short_1_3s": [1.0, 3.0],
        "sec_microseism_3_10s": [3.0, 10.0],
        "prim_microseism_10_15s": [10.0, 15.0],
        "notch_15_50s": [15.0, 50.0],
        "compliance_50_100s": [50.0, 100.0],
        "ig_100_300s": [100.0, 300.0],
        "long_300_1000s": [300.0, 1000.0],
    },
    "depth_bins_m": [0, 200, 1000, 2000, 3000, 4000, 5000, 100000],
    "depth_labels": [
        "<200 m",
        "200-1000 m",
        "1-2 km",
        "2-3 km",
        "3-4 km",
        "4-5 km",
        ">5 km",
    ],
    "distance_bins_km": [0, 50, 200, 500, 100000],
    "distance_labels": ["<50 km", "50-200 km", "200-500 km", ">500 km"],
    "shelf_break_bins_km": [-1000, 0, 50, 200, 500, 100000],
    "shelf_break_labels": [
        "on the shelf",
        "0-50 km",
        "50-200 km",
        "200-500 km",
        ">500 km",
    ],
    "stats": {"seed": 20261118, "n_bootstrap": 1000},
}


def _merge(base, over):
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        out[k] = (
            _merge(out[k], v)
            if isinstance(v, dict) and isinstance(out.get(k), dict)
            else v
        )
    return out


def load(path=None, **overrides):
    """Settings: the defaults, then a JSON/YAML file, then keyword overrides (nested dicts)."""
    cfg = copy.deepcopy(DEFAULTS)
    if path:
        text = Path(path).read_text()
        if str(path).endswith((".yaml", ".yml")):
            import yaml

            data = yaml.safe_load(text)
        else:
            data = json.loads(text)
        cfg = _merge(cfg, data)
    return _merge(cfg, overrides)


def db_edges(cfg):
    import numpy as np

    h = cfg["histogram"]
    return np.arange(h["db_min"], h["db_max"] + 0.5 * h["db_step"], h["db_step"])
