"""Deprecated: ComPy 2.x module, split into compy.inversion (and compy.plotting.inversion).

Every old name still resolves to its new home; renamed ones warn with the new
name. See ``compy._legacy.LEGACY["inv_compy"]`` for the full table.
"""

from compy._deprecation import legacy_getattr, legacy_module
from compy._legacy import LEGACY

legacy_module("inv_compy", "compy.inversion (and compy.plotting.inversion)")

_TABLE = LEGACY["inv_compy"]
__getattr__ = legacy_getattr(
    __name__,
    _TABLE,
    quiet={old for old, new in _TABLE.items() if new.endswith(":" + old)},
)


def __dir__():
    return sorted(_TABLE)
