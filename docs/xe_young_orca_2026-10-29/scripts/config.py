"""Run configuration: XE (Pacific Array Young ORCA, 2018-2019) broadband compliance.

Every stage script imports this module, so a new network only needs a new
config. The values mirror the YV run of 2026-10-28 (Aminian et al. 2025)
except where the instruments differ (50 Hz sampling, empty location code).
"""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPY = ROOT.parent / "ComPy"
WORK = ROOT / "work"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

# --- data -------------------------------------------------------------------
NETWORK = "XE"
# COMPY_SDS points stage A at a local copy (a_fast.py) of a slow archive.
SDS = Path(os.environ.get("COMPY_SDS", "/Volumes/DATA_bruitfm/SDS"))
SDS_LOCATION = "00"  # location code in the SDS file names (traces carry "")
YEARS = (2018, 2019)
INVENTORY = ROOT / "metadata" / "XE.xml"  # IRIS seismic + pressure responses
Z, H1, H2, P = "BHZ", "BH1", "BH2", "BDH"
CHANNELS = (Z, H1, H2, P)

# Stations with more than 3 months of data and a working DPG (BRUIT-FM QC):
# CC02 (pressure non-working), CC10 (50 d) and CC12 (17 d) are left out.
STATIONS = (
    "CC01", "CC03", "CC04", "CC05", "CC06", "CC07", "CC08", "CC09", "CC11",
    "EC01", "EC02", "EC03", "EC04", "EC05",
    "EE01", "EE02", "EE03", "EE04",
    "WC01", "WC02", "WC03", "WC04", "WC05",
    "WW01", "WW02", "WW03", "WW04",
)  # fmt: skip

# --- stage A: decimation and responses --------------------------------------
RAW_RATE = 50.0
FACTORS = (5, 4, 2)  # 50 Hz -> 1.25 Hz with TiSKitPy SAC FIR filters
OUT_RATE = 1.25

# --- stage B: cleaning --------------------------------------------------------
# Periodic-transient (glitch) removal. The YV INSU-IPGP OBS have an hourly
# mass-position glitch; the Scripps ABALONE OBS of XE do not: 584 15-day blocks
# at CC04, CC05, CC06, CC11, EE04 and WC03 gave fold strengths of 1.15 (median)
# and 1.76 (max), below the 2.0 needed to subtract a template. Detection costs
# about 1 h per station, so it is switched off for the remaining XE stations.
DEGLITCH = False

# --- catalogue --------------------------------------------------------------
CATALOG = WORK / "catalog_usgs_M5.85_20180401_20190630.qml"
CATALOG_SPAN = ("2018-04-01", "2019-06-30")

# --- gains ------------------------------------------------------------------
BRUITFM_GAINS = ROOT.parent / "gain_corrected_2026-10-28" / "gains.csv"


def seed_id(station, channel=P):
    """SEED id as written in the BRUIT-FM tables."""
    return f"{NETWORK}.{station}..{channel}"


# --- inversion ----------------------------------------------------------------
CRUST1 = ROOT.parent / "compliance_yv_2026-10-28" / "literature" / "crust1"

# --- figures --------------------------------------------------------------------
GROUPS = {"WW": "West (WW)", "WC": "West-central (WC)", "CC": "Central (CC)",
          "EC": "East-central (EC)", "EE": "East (EE)"}  # fmt: skip
