"""Deprecated: ComPy 2.x module, split into compy.plotting.spectra.

Every old name still resolves to its new home; renamed ones warn with the new
name. See ``compy._legacy.LEGACY["ffplot"]`` for the full table.
"""

from compy._deprecation import legacy_getattr, legacy_module
from compy._legacy import LEGACY

legacy_module("ffplot", "compy.plotting.spectra")

_TABLE = LEGACY["ffplot"]
__getattr__ = legacy_getattr(
    __name__,
    _TABLE,
    quiet={old for old, new in _TABLE.items() if new.endswith(":" + old)},
)


def __dir__():
    return sorted(_TABLE)
