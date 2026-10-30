"""Validated tabulated dispersion or explicit layered-model disba generation."""

import hashlib
from importlib.metadata import version
from io import StringIO
from pathlib import Path

import numpy as np

LAYER_COLUMNS = "thickness_km,vp_km_s,vs_km_s,density_g_cm3"
DISPERSION_COLUMNS = "frequency_hz,phase_velocity_m_s,group_velocity_m_s"


def validate_dispersion(table, band):
    table = np.asarray(table, dtype=float)
    if table.ndim != 2 or table.shape[1] != 3 or len(table) < 2:
        raise ValueError("Expected frequency, phase velocity, group velocity CSV")
    if not np.isfinite(table).all() or np.any(table <= 0):
        raise ValueError("Dispersion must be finite and positive")
    if np.any(np.diff(table[:, 0]) <= 0):
        raise ValueError("Dispersion frequency must be strictly increasing")
    if band[0] < table[0, 0] or band[1] > table[-1, 0]:
        raise ValueError("Dispersion cannot be extrapolated into the fit band")
    if not table[0, 0] <= 0.04 <= table[-1, 0]:
        raise ValueError("Dispersion must bracket the 40 mHz packet")
    return table


def dispersion_from_layers(path, band):
    path = Path(path)
    content = path.read_bytes()
    lines = content.decode().splitlines()
    if not lines or lines[0].strip() != LAYER_COLUMNS:
        raise ValueError(f"Earth model CSV header must be {LAYER_COLUMNS}")
    layers = np.loadtxt(StringIO(content.decode()), delimiter=",", skiprows=1, ndmin=2)
    if layers.shape[1] != 4 or len(layers) < 2 or not np.isfinite(layers).all():
        raise ValueError("Earth model requires at least two finite four-column layers")
    if np.any(layers[:-1, 0] <= 0) or layers[-1, 0] != 0:
        raise ValueError(
            "Layer thicknesses must be positive, with a zero-thickness half-space"
        )
    if (
        np.any(layers[:, 1] <= 0)
        or np.any(layers[:, 3] <= 0)
        or np.any(layers[:, 2] < 0)
    ):
        raise ValueError("Earth model velocities/densities must be physical and finite")
    if np.any(layers[1:, 2] == 0) or np.any(layers[:, 1] <= layers[:, 2]):
        raise ValueError("Only the first layer may be fluid; require Vp > Vs")
    try:
        from disba import GroupDispersion, PhaseDispersion
    except ImportError as exc:
        raise ImportError(
            "Layer-model option requires installation of .[dispersion]"
        ) from exc
    fmin, fmax = min(0.02, band[0]), max(0.1, band[1])
    periods = np.geomspace(1 / fmax, 1 / fmin, 160)
    phase = PhaseDispersion(*layers.T)(periods, mode=0, wave="rayleigh")
    group = GroupDispersion(*layers.T)(periods, mode=0, wave="rayleigh")
    if (
        len(phase.period) != len(periods)
        or len(group.period) != len(periods)
        or not np.allclose(phase.period, periods, rtol=1e-12, atol=0)
        or not np.allclose(group.period, periods, rtol=1e-12, atol=0)
    ):
        raise ValueError("disba did not solve the full requested period range")
    frequencies = 1 / periods[::-1]
    frequencies[[0, -1]] = [fmin, fmax]
    table = np.column_stack(
        (frequencies, phase.velocity[::-1] * 1000, group.velocity[::-1] * 1000)
    )
    if path.read_bytes() != content:
        raise ValueError("Earth model changed during dispersion calculation")
    provenance = dict(
        source_kind="layer_model_disba",
        source_file=str(path),
        source_sha256=hashlib.sha256(content).hexdigest(),
        disba_version=version("disba"),
        numba_version=version("numba"),
        algorithm="dunkin",
        mode=0,
        wave="rayleigh",
        layers_units="km, km/s, km/s, g/cm3",
        period_seconds=[float(periods[0]), float(periods[-1])],
        period_count=160,
        field_validation="not_established_by_solver",
    )
    return validate_dispersion(table, band), provenance


def load_dispersion(*, dispersion_file=None, earth_model_file=None, band):
    if len(band) != 2 or not np.isfinite(band).all() or not 0 < band[0] < band[1]:
        raise ValueError("Invalid dispersion band")
    if bool(dispersion_file) == bool(earth_model_file):
        raise ValueError("Supply exactly one dispersion_file or earth_model_file")
    if earth_model_file:
        return dispersion_from_layers(earth_model_file, band)
    path = Path(dispersion_file)
    table = validate_dispersion(np.loadtxt(path, delimiter=",", skiprows=1), band)
    return table, dict(
        source_kind="dispersion_table",
        source_file=str(path),
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    )
