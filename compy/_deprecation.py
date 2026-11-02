"""Helpers that keep pre-package names importable, with a DeprecationWarning."""

import importlib
import warnings

REMOVAL = "ComPy 3.0"


def _warn(old, new):
    warnings.warn(
        f"{old} is deprecated and will be removed in {REMOVAL}; use {new}",
        DeprecationWarning,
        stacklevel=3,
    )


def resolve(target):
    """Return the object named ``"package.module:name"``."""
    module, _, name = target.partition(":")
    return getattr(importlib.import_module(module), name)


def legacy_getattr(module_name, table, *, quiet=(), submodules=()):
    """Build a module ``__getattr__`` serving the old names in ``table``.

    ``table`` maps an old name to ``"new.module:new_name"``. Access warns,
    except for names in ``quiet`` (current public API re-exported unchanged).
    Names in ``submodules`` are imported lazily as ``module_name.name``.
    """

    def __getattr__(name):
        if name in submodules:
            return importlib.import_module(f"{module_name}.{name}")
        try:
            target = table[name]
        except KeyError:
            raise AttributeError(
                f"module {module_name!r} has no attribute {name!r}"
            ) from None
        if name not in quiet:
            module, _, new = target.partition(":")
            _warn(f"{module_name}.{name}", f"{module}.{new}")
        return resolve(target)

    return __getattr__


def renamed(module_name, names):
    """Module ``__getattr__`` for names renamed within ``module_name``."""
    return legacy_getattr(
        module_name, {old: f"{module_name}:{new}" for old, new in names.items()}
    )


def legacy_module(old, new):
    """Warn that the top-level module ``old`` moved to ``new``."""
    warnings.warn(
        f"the {old} module is deprecated and will be removed in {REMOVAL}; "
        f"import {new} instead",
        DeprecationWarning,
        stacklevel=3,
    )
