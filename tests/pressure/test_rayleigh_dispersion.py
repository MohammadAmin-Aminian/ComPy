"""Independent half-space dispersion and source/provenance guards."""

import hashlib

import numpy as np
import pytest
from scipy.optimize import brentq
from test_rayleigh_calibration import synthetic_station as synthetic_station

from compy_pressure.rayleigh_dispersion import LAYER_COLUMNS, load_dispersion


def write_model(path, layers):
    np.savetxt(path, layers, delimiter=",", header=LAYER_COLUMNS, comments="")


def test_disba_homogeneous_model_matches_rayleigh_secular_equation(tmp_path):
    pytest.importorskip("disba")
    path = tmp_path / "layers.csv"
    write_model(path, [[1, 6, 3, 2.7], [0, 6, 3, 2.7]])
    original = path.read_bytes()
    table, provenance = load_dispersion(earth_model_file=path, band=(0.04, 0.1))
    ratio = brentq(
        lambda x: (2 - x * x) ** 2 - 4 * np.sqrt(1 - x * x) * np.sqrt(1 - 0.25 * x * x),
        0.5,
        0.999,
    )
    # Nondispersive homogeneous half-space: group = phase, in m/s.
    np.testing.assert_allclose(table[:, 1:], ratio * 3000, rtol=2e-4)
    assert table[0, 0] == 0.02 and table[-1, 0] == 0.1
    assert np.all(np.diff(table[:, 0]) > 0)
    assert provenance["mode"] == 0
    assert provenance["source_sha256"] == hashlib.sha256(original).hexdigest()
    assert path.read_bytes() == original


@pytest.mark.parametrize(
    "layers",
    [
        [[1, 6, 3, 2.7], [1, 6, 3, 2.7]],
        [[-1, 6, 3, 2.7], [0, 6, 3, 2.7]],
        [[1, 6, 3, 0], [0, 6, 3, 2.7]],
        [[1, 6, 3, 2.7], [0, 6, 0, 2.7]],
    ],
)
def test_invalid_layers_are_not_silently_repaired(tmp_path, layers):
    path = tmp_path / "layers.csv"
    write_model(path, layers)
    with pytest.raises(ValueError):
        load_dispersion(earth_model_file=path, band=(0.04, 0.1))


def test_model_source_is_exclusive_and_table_is_not_extrapolated(tmp_path):
    path = tmp_path / "table.csv"
    np.savetxt(
        path,
        [[0.02, 4000, 3800], [0.085, 4000, 3800]],
        delimiter=",",
        header="f,c,u",
        comments="",
    )
    with pytest.raises(ValueError, match="exactly one"):
        load_dispersion(dispersion_file=path, earth_model_file=path, band=(0.04, 0.1))
    with pytest.raises(ValueError, match="extrapolated"):
        load_dispersion(dispersion_file=path, band=(0.04, 0.1))


def test_network_model_input_changes_cache_dependency(synthetic_station, tmp_path):
    from compy_pressure import network
    from compy_pressure.inventory import load_master_dataset
    from compy_pressure.rayleigh_network import dependency_stamp

    data, exports, catalog, _ = synthetic_station
    model = tmp_path / "layers.csv"
    write_model(model, [[1, 6, 3, 2.7], [0, 6, 3, 2.7]])
    master = load_master_dataset(exports)
    row = master[master.input_units.eq("PA")].iloc[0]
    files = network.scan_files(data)
    options = dict(catalog=str(catalog), earth_model_file=str(model))
    before = dependency_stamp(row, master, files, options)
    write_model(model, [[1, 6, 3.1, 2.7], [0, 6, 3.1, 2.7]])
    assert before != dependency_stamp(row, master, files, options)
    assert dependency_stamp(
        row, master, files, options, {"disba_version": "a"}
    ) != dependency_stamp(row, master, files, options, {"disba_version": "b"})


def test_layer_model_calibration_saves_generated_table_and_preserves_master(
    synthetic_station, tmp_path
):
    pytest.importorskip("disba")
    from compy_pressure.rayleigh_calibration import run_calibration

    data, exports, catalog, _ = synthetic_station
    before = {p: p.read_bytes() for p in exports.glob("master_inventory.*")}
    ratio = brentq(
        lambda x: (2 - x * x) ** 2 - 4 * np.sqrt(1 - x * x) * np.sqrt(1 - 0.25 * x * x),
        0.5,
        0.999,
    )
    vs = 4 / ratio
    model = tmp_path / "layers.csv"
    write_model(model, [[1, 2 * vs, vs, 2.7], [0, 2 * vs, vs, 2.7]])
    result = run_calibration(
        catalog=catalog,
        earth_model_file=model,
        master_dir=exports,
        data_dir=data,
        network="XX",
        station="TEST",
        pressure_sign=1,
    )
    assert result["events"][0]["status"] == "candidate"
    assert result["events"][0]["pressure_multiplier"] == pytest.approx(0.67, rel=0.02)
    assert result["dispersion_model"]["source_kind"] == "layer_model_disba"
    assert str(model) in result["inputs"]
    generated = np.loadtxt(result["files"]["dispersion"], delimiter=",", skiprows=1)
    assert generated[-1, 0] == 0.1
    assert before == {p: p.read_bytes() for p in before}
    from compy_pressure import network, rayleigh_network

    report = rayleigh_network.run_network(
        exports,
        tmp_path / "network-model-results",
        network.scan_files(data),
        update_master=False,
        options=dict(
            catalog=catalog, earth_model_file=model, pressure_sign=1, band=[0.02, 0.06]
        ),
    )
    assert report["dispersion_model"]["source_kind"] == "layer_model_disba"
    entry = next(iter(report["entries"].values()))
    assert entry["status"] == "unresolved"  # Single event, not a model failure.
    assert any(p.endswith("model_dispersion.csv") for p in entry["files"])
    assert before == {p: p.read_bytes() for p in before}
