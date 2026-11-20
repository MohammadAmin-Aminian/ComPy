"""Entry points of the ``compy-pressure*`` commands.

The pressure tools need pandas and openpyxl (the ``pressure`` extra) and the
maps need cartopy (the ``maps`` extra). A command started without them ends
with one line that names the extra instead of an import traceback.
"""

from __future__ import annotations

import importlib

#: Optional module -> the extra that installs it.
EXTRAS = {"openpyxl": "pressure", "pandas": "pressure", "cartopy": "maps"}


def run(module, argv=None):
    """Import ``module`` and run its ``main``; name the missing extra if an import fails."""
    try:
        target = importlib.import_module(module)
        return target.main() if argv is None else target.main(argv)
    except ModuleNotFoundError as exc:
        missing = (exc.name or "").split(".")[0]
        if missing not in EXTRAS:
            raise
        raise SystemExit(
            f"This command needs '{missing}', which the '{EXTRAS[missing]}' extra installs: "
            f"python -m pip install 'seafloor-compy[{EXTRAS[missing]}]'"
        ) from None


def network():
    """``compy-pressure``."""
    return run("compy_pressure.network")


def inventory():
    """``compy-pressure-inventory``."""
    return run("compy_pressure.inventory")


def tilt():
    """``compy-pressure-tilt``."""
    return run("compy_pressure.tilt")


def local_p():
    """``compy-pressure-local-p``."""
    return run("compy_pressure.local_p_check")


def microseism():
    """``compy-pressure-microseism``."""
    return run("compy_pressure.microseism_check")


def paired():
    """``compy-pressure-paired``."""
    return run("compy_pressure.paired_spectra")


def noise_model():
    """``compy-pressure-noise-model``."""
    return run("compy_pressure.group_noise_model")


def approved_noise_model():
    """``compy-pressure-approved-noise-model``."""
    return run("compy_pressure.noise_model")


def atlas():
    """``compy-pressure-atlas``."""
    return run("compy_pressure.atlas.cli")
