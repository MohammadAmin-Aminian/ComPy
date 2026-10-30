"""compy_pressure: ocean-bottom pressure-gauge QC, calibration and noise models.

Exported from BRUIT-FM 1.3.0 (https://github.com/MohammadAmin-Aminian/bruit-fm)
by its integration/export_to_compy.py. Until BRUIT-FM is frozen, change the
code there and export it again: a change made only here is overwritten.
"""

from importlib.metadata import PackageNotFoundError, version

BRUITFM_VERSION = "1.3.0"
BRUITFM_COMMIT = "356f840a8d263d1ccf8a683302d0d16a7c15dca4"

try:
    __version__ = version("seafloor-compy")
except PackageNotFoundError:  # running from a source tree that is not installed
    __version__ = "0+bruitfm1.3.0"
