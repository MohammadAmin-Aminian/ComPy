"""Deprecated alias of :mod:`compy.cleaning` (ComPy 2.x module name)."""

import sys

from compy._deprecation import legacy_module

legacy_module("compy_cleaning", "compy.cleaning")

import compy.cleaning as _module  # noqa: E402

sys.modules[__name__] = _module
