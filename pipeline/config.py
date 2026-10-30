"""Run configuration for the ComPy end-to-end pipeline.

Every stage script imports this module, so a new network needs only new
values here. The defaults are the validated XE run (Pacific Array Young ORCA,
2018-2019, docs/xe_young_orca_2026-10-29); the YV values of Aminian et al.
(2025) are given in comments.

Two environment variables avoid editing this file:

* ``COMPY_RUN_DIR``: the run folder (work/, results/, figures/, logs/,
  metadata/); default the current directory.
* ``COMPY_PIPELINE_CONFIG``: a Python file executed after these defaults,
  which may override any name below (e.g. NETWORK, STATIONS, INVENTORY).
"""

import os
from pathlib import Path

ROOT = Path(os.environ.get("COMPY_RUN_DIR", Path.cwd())).resolve()
COMPY = Path(__file__).resolve().parents[1]  # the ComPy repository
WORK = ROOT / "work"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

# --- data -------------------------------------------------------------------
NETWORK = "XE"  # YV
# SDS archive; COMPY_SDS points stage A at a local copy (a_fast.py).
SDS = Path(os.environ.get("COMPY_SDS", "/Volumes/DATA_bruitfm/SDS"))
SDS_LOCATION = "00"  # location code in the SDS file names
YEARS = (2018, 2019)  # YV: (2012, 2013)
INVENTORY = ROOT / "metadata" / f"{NETWORK}.xml"  # seismometer and DPG responses
Z, H1, H2, P = "BHZ", "BH1", "BH2", "BDH"
CHANNELS = (Z, H1, H2, P)
STATIONS = (
    "CC01", "CC03", "CC04", "CC05", "CC06", "CC07", "CC08", "CC09", "CC11",
    "EC01", "EC02", "EC03", "EC04", "EC05",
    "EE01", "EE02", "EE03", "EE04",
    "WC01", "WC02", "WC03", "WC04", "WC05",
    "WW01", "WW02", "WW03", "WW04",
)  # fmt: skip

# --- stage A: decimation and responses --------------------------------------
RAW_RATE = 50.0  # YV: 62.5
FACTORS = (5, 4, 2)  # TiSKitPy FIR decimation to 1.25 Hz; YV: (5, 5, 2)
OUT_RATE = 1.25

# --- stage B: cleaning --------------------------------------------------------
# Periodic-transient (hourly glitch) removal, needed for the INSU-IPGP OBS of
# YV. The Scripps ABALONE OBS of XE show none (fold strength below 2.0 in all
# 584 15-day blocks tested), and detection costs about 1 h per station.
DEGLITCH = False  # YV: True

# --- catalogue (stage 0, fetch_catalog.py) -----------------------------------
CATALOG_SPAN = ("2018-04-01", "2019-06-30")
CATALOG = WORK / f"catalog_usgs_M5.85_{CATALOG_SPAN[0]}_{CATALOG_SPAN[1]}.qml"

# --- gains ----------------------------------------------------------------------
# Optional BRUIT-FM gain table (columns Seed_id, gain_status, gain_applied),
# compared with ComPy's own Rayleigh-wave calibration. None to skip.
BRUITFM_GAINS = None

# --- inversion ------------------------------------------------------------------
CRUST1 = ROOT / "metadata" / "crust1"  # crust1.bnds, .vp, .vs, .rho


# --- figures (f_figures.py) -----------------------------------------------------
# Station-code prefix -> sub-array label for Figs. 8 and 10. Optional files:
# metadata/bathy/etopo_<net>.nc (map) and metadata/bathy/etopo_region.nc (inset).
GROUPS = {"WW": "West (WW)", "WC": "West-central (WC)", "CC": "Central (CC)",
          "EC": "East-central (EC)", "EE": "East (EE)"}  # fmt: skip


def seed_id(station, channel=P):
    """SEED id as written in the BRUIT-FM tables."""
    return f"{NETWORK}.{station}..{channel}"


_override = os.environ.get("COMPY_PIPELINE_CONFIG")
if _override:
    exec(compile(Path(_override).read_text(), _override, "exec"))
