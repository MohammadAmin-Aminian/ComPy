"""Deprecated alias of :mod:`compy.compliance` (ComPy 2.x module name)."""

import sys

from compy._deprecation import legacy_module

legacy_module("compy_processing", "compy.compliance")

import compy.compliance as _module  # noqa: E402

sys.modules[__name__] = _module
