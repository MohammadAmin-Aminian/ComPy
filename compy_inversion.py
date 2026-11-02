"""Deprecated alias of :mod:`compy.inversion.dram` (ComPy 2.x module name)."""

import sys

from compy._deprecation import legacy_module

legacy_module("compy_inversion", "compy.inversion.dram")

import compy.inversion.dram as _module  # noqa: E402

sys.modules[__name__] = _module
