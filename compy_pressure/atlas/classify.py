"""Classification rules of the pressure atlas: sensor type and brand, basin, depth and distance classes, season.

Sensor type and basin follow the rules of the ComPy whole-archive run
(``scripts/build_inventory.py``), which are those of the PressureQC report of
28 October 2026; the brand rule was written for the atlas of 6 November 2026.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def sensor_type(description, channel=""):
    """APG, DPG, Hydrophone or Other from the StationXML sensor description and the channel code."""
    s = str(description).lower() if pd.notna(description) else ""
    if "parosci" in s or "apg" in s:
        return "APG"
    if any(k in s for k in ("differential", "dpg", "pascals gauge", "cox-webb")):
        return "DPG"
    if any(k in s for k in ("hydrophone", "hti", "hitech", "hitec", "oas", "itc")):
        return "Hydrophone"
    if channel == "HDH" and not s.strip():
        return "Hydrophone"
    return "Other"


def brand(description, kind):
    """Sensor brand or model from the description and the sensor type."""
    s = str(description).lower() if pd.notna(description) else ""
    if kind == "APG":
        return "Paroscientific APG" if "parosci" in s else "Unknown APG"
    if kind == "DPG":
        if "dpg5004" in s:
            return "INSU DPG5004"
        if "ldeo" in s:
            return "LDEO DPG"
        if "sio" in s or "cox-webb" in s:
            return "SIO Cox-Webb DPG"
        return "Generic DPG"
    if kind == "Hydrophone":
        for key, name in (("oas", "OAS hydrophone"), ("itc", "ITC hydrophone"), ("edo", "EDO hydrophone"),
                          ("ssq", "HTI-94-SSQ hydrophone (2-30 kHz)")):  # fmt: skip
            if key in s:
                return name
        if any(
            k in s for k in ("hti", "ht1", "hitech", "hitec", "high tech", "hightech")
        ):
            return "HTI hydrophone"
        return "Generic hydrophone"
    return "Q330 internal" if "q330" in s else "Other"


def basin(lat, lon):
    """Ocean basin of a position (coarse rule, enough for station grouping)."""
    if lat > 66:
        return "Arctic"
    if lat < -50:
        return "Southern"
    if 30 < lat < 46 and -6 < lon < 37:
        return "Mediterranean"
    if -70 < lon < 20 or (lat > 10 and -98 < lon < 20):
        return "Atlantic"
    if 20 <= lon < 147 and lat < 30:
        return "Indian"
    return "Pacific"


def binned(values, edges, labels):
    """Class label of each value, intervals closed on the left."""
    return pd.cut(
        np.asarray(values, float), bins=edges, labels=labels, right=False
    ).astype(object)


SEASONS = ("winter", "spring", "summer", "autumn")


def local_season(month, lat):
    """Season of a calendar month (1-12) at a latitude: DJF is winter in the north, JJA in the south."""
    m = month if lat >= 0 else (month + 5) % 12 + 1
    return {12: "winter", 1: "winter", 2: "winter", 3: "spring", 4: "spring", 5: "spring",
            6: "summer", 7: "summer", 8: "summer", 9: "autumn", 10: "autumn", 11: "autumn"}[m]  # fmt: skip


def season_months(season, lat):
    """Month indices (0 = January) of a local season at a latitude."""
    return np.array([m for m in range(12) if local_season(m + 1, lat) == season])


def add_classes(table, cfg):
    """Add ``basin``, ``depth_bin``, ``distance_bin`` and ``shelf_bin`` columns to a channel table (copy)."""
    t = table.copy()
    if "basin" not in t:
        t["basin"] = [basin(a, b) for a, b in zip(t.lat, t.lon)]
    t["depth_bin"] = binned(t.depth_m, cfg["depth_bins_m"], cfg["depth_labels"])
    if "coast_km" in t:
        t["distance_bin"] = binned(
            t.coast_km, cfg["distance_bins_km"], cfg["distance_labels"]
        )
    if "shelf_break_km" in t:
        t["shelf_bin"] = binned(
            t.shelf_break_km, cfg["shelf_break_bins_km"], cfg["shelf_break_labels"]
        )
    return t
