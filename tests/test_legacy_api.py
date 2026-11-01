"""ComPy 2.x names keep working after the move to the compy package."""

import importlib
import warnings

import pytest

import compy
from compy._deprecation import resolve
from compy._legacy import LEGACY

ALIAS_MODULES = {
    "compy_numerics": "compy.numerics",
    "compy_streams": "compy.streams",
    "compy_cleaning": "compy.cleaning",
    "compy_diagnostics": "compy.diagnostics",
    "compy_inversion": "compy.inversion.dram",
    "compy_processing": "compy.compliance",
}

CASES = [
    (module, old, new) for module, table in LEGACY.items() for old, new in table.items()
]


@pytest.mark.parametrize("module, old, target", CASES)
def test_every_old_name_resolves_to_its_new_home(module, old, target):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        legacy = importlib.import_module(module)
        assert getattr(legacy, old) is resolve(target)


RENAMED = [(m, o, n) for m, o, n in CASES if n.partition(":")[2] != o]


@pytest.mark.parametrize("module, old, target", RENAMED)
def test_renamed_names_warn_with_the_new_name(module, old, target):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        legacy = importlib.import_module(module)
    with pytest.warns(DeprecationWarning, match=target.partition(":")[2]):
        getattr(legacy, old)


@pytest.mark.parametrize("old, new", sorted(ALIAS_MODULES.items()))
def test_old_module_names_are_the_new_modules(old, new):
    # The same module object, so private names and monkeypatching still work.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        assert importlib.import_module(old) is importlib.import_module(new)


def test_paper_api_on_the_package():
    with pytest.warns(DeprecationWarning, match="calculate_compliance"):
        assert compy.Calculate_Compliance is compy.calculate_compliance
    with pytest.warns(DeprecationWarning, match="correct_tilt"):
        assert compy.Rotate is compy.correct_tilt
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # current names are silent
        assert compy.remove_coherent_noise is compy.preprocessing.remove_coherent_noise
        assert compy.wavenumber is compy.numerics.wavenumber


def test_compliance_module_keeps_estimate():
    from compy import compliance

    with pytest.warns(DeprecationWarning, match="estimate_compliance"):
        assert compliance.estimate is compliance.estimate_compliance


def test_unknown_names_still_raise():
    with pytest.raises(AttributeError):
        compy.no_such_function  # noqa: B018
