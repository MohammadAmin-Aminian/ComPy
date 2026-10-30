"""compy_pressure: ocean-bottom pressure-gauge QC, calibration and noise models.

Moved from BRUIT-FM 1.3.0 (https://github.com/MohammadAmin-Aminian/bruit-fm,
commit ba735b90af6ec3b9556a60527387679ef4697ae5) by its integration/export_to_compy.py. BRUIT-FM is kept as a
backup; this package is the maintained copy.
"""

from importlib.metadata import PackageNotFoundError, version

BRUITFM_VERSION = "1.3.0"
BRUITFM_COMMIT = "ba735b90af6ec3b9556a60527387679ef4697ae5"

try:
    __version__ = version("seafloor-compy")
except PackageNotFoundError:  # running from a source tree that is not installed
    __version__ = "0+bruitfm1.3.0"
