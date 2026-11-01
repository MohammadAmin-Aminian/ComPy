"""Deprecated alias of :mod:`compy.diagnostics` (ComPy 2.x module name)."""

import sys

from compy._deprecation import legacy_module

legacy_module("compy_diagnostics", "compy.diagnostics")

import compy.diagnostics as _module  # noqa: E402

sys.modules[__name__] = _module
