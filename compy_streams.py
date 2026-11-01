"""Deprecated alias of :mod:`compy.streams` (ComPy 2.x module name)."""

import sys

from compy._deprecation import legacy_module

legacy_module("compy_streams", "compy.streams")

import compy.streams as _module  # noqa: E402

sys.modules[__name__] = _module
