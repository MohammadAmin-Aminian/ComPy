"""Deprecated alias of :mod:`compy.numerics` (ComPy 2.x module name)."""

import sys

from compy._deprecation import legacy_module

legacy_module("compy_numerics", "compy.numerics")

import compy.numerics as _module  # noqa: E402

sys.modules[__name__] = _module
