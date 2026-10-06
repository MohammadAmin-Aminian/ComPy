import ast
from pathlib import Path

import numpy as np


def _functions(path, names):
    tree = ast.parse(Path(path).read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    module = ast.Module(body=nodes, type_ignores=[])
    ns = {"np": np}
    exec(compile(module, path, "exec"), ns)
    return ns


def test_gaussian_likelihood_uses_squared_normalized_residual():
    ns = _functions("inv_compy.py", {"liklihood"})
    d = np.array([1.0, 2.0])
    m = np.array([0.0, 0.0])
    s = np.array([1.0, 2.0])
    expected = np.exp(-0.5 * np.sum(((d - m) / s) ** 2))
    assert np.isclose(ns["liklihood"](d, m, s=s), expected)


def test_roughness_is_nonnegative_and_zero_for_constant_profile():
    ns = _functions("inv_compy.py", {"Roughness"})
    assert np.isclose(ns["Roughness"](np.ones(20), 1), 0.0)
    assert ns["Roughness"](np.linspace(0, 1, 20) ** 2, 2) >= 0.0


def test_pressure_calibration_retains_intended_band():
    source = Path("Pressure_calibration.py").read_text()
    hp = 'stream22.filter("highpass", freq=freq1)'
    lp = 'stream22.filter("lowpass", freq=freq2)'
    assert hp in source and lp in source
    assert source.index(hp) < source.index(lp)
