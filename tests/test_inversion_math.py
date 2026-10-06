import ast
from pathlib import Path

import numpy as np


def _load_functions(path, names):
    tree = ast.parse(Path(path).read_text())
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    ns = {"np": np}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), path, "exec"), ns)
    return ns


def test_likelihood_is_one_for_perfect_fit():
    ns = _load_functions("inv_compy.py", {"liklihood"})
    data = np.array([1.0, 2.0, 3.0])
    assert np.isclose(ns["liklihood"](data, data, s=0.1), 1.0)


def test_likelihood_decreases_with_residual():
    ns = _load_functions("inv_compy.py", {"liklihood"})
    data = np.zeros(3)
    near = ns["liklihood"](data, np.array([0.1, 0.0, 0.0]), s=1.0)
    far = ns["liklihood"](data, np.array([1.0, 0.0, 0.0]), s=1.0)
    assert 0.0 < far < near <= 1.0


def test_roughness_constant_profile_zero():
    ns = _load_functions("inv_compy.py", {"Roughness"})
    assert np.isclose(ns["Roughness"](np.full(32, 2500.0), 2), 0.0)


def test_roughness_nonnegative():
    ns = _load_functions("inv_compy.py", {"Roughness"})
    profile = np.array([1500.0, 1700.0, 2200.0, 2800.0, 3000.0])
    assert ns["Roughness"](profile, 1) >= 0.0
    assert ns["Roughness"](profile, 2) >= 0.0


def test_velocity_to_vp_poisson_quarter():
    ns = _load_functions("inv_compy.py", {"velp"})
    vs = np.array([1000.0, 2000.0])
    expected = vs * np.sqrt(3.0)
    np.testing.assert_allclose(ns["velp"](vs, p=0.25), expected)


def test_density_positive_for_positive_velocity():
    ns = _load_functions("inv_compy.py", {"density", "velp"})
    rho = ns["density"](np.array([500.0, 2000.0, 4000.0]))
    assert np.all(np.isfinite(rho))
    assert np.all(rho > 0)
