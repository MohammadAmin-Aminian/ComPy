"""Multiple-chain diagnostics and portable run records for ComPy.

Diagnostics use ArviZ's rank-normalized, folded split-Rhat and bulk/tail ESS.
Install ``seafloor-compy[diagnostics]`` to use ``summarize_chains``.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import os
from tempfile import TemporaryDirectory
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


def summarize_chains(draws, parameter_names, *, burnin=0):
    """Summarize ``(chains, iterations, parameters)`` draws after burn-in.

    Use independent chains with the same posterior target. Undefined diagnostics
    are returned as None with a warning, never as successful convergence.
    Quantiles describe the retained sample, even when it has not converged.
    """
    values = np.asarray(draws, dtype=float)
    names = list(parameter_names)
    if (
        values.ndim != 3
        or values.shape[0] < 2
        or values.shape[2] == 0
        or not np.isfinite(values).all()
        or len(names) != values.shape[2]
        or any(not isinstance(name, str) or not name for name in names)
        or len(set(names)) != len(names)
    ):
        raise ValueError(
            "Need finite (chains, iterations, parameters) draws and unique names"
        )
    if (
        isinstance(burnin, (bool, np.bool_))
        or not isinstance(burnin, (int, np.integer))
        or not 0 <= burnin <= values.shape[1] - 4
    ):
        raise ValueError("burnin must leave at least four iterations per chain")
    try:
        import arviz as az
    except ImportError as exc:
        raise ImportError(
            "Install ComPy with '.[diagnostics]' for chain diagnostics"
        ) from exc
    retained = values[:, burnin:, :]
    messages = []
    if len(values) < 4:
        messages.append("Use at least four independent chains for a final analysis.")
    if any(
        np.array_equal(retained[i], retained[j])
        for i in range(len(values))
        for j in range(i)
    ):
        messages.append("Identical chains detected; check that seeds are independent.")
    summaries = {}
    for index, name in enumerate(names):
        samples = retained[:, :, index]
        quantiles = np.quantile(samples, [0.05, 0.5, 0.95])
        entry = dict(zip(("q05", "median", "q95"), map(float, quantiles)))
        if np.any(np.ptp(samples, axis=1) == 0):
            entry.update(rhat=None, ess_bulk=None, ess_tail=None)
            messages.append(f"{name}: a chain is constant; diagnostics are undefined.")
        else:
            diagnostics = {
                "rhat": az.rhat(samples, method="rank"),
                "ess_bulk": az.ess(samples, method="bulk"),
                "ess_tail": az.ess(samples, method="tail"),
            }
            entry.update(
                {
                    key: float(value) if np.isfinite(value) else None
                    for key, value in diagnostics.items()
                }
            )
            if any(entry[key] is None for key in diagnostics):
                messages.append(f"{name}: a diagnostic is undefined.")
            if entry["rhat"] is not None and entry["rhat"] >= 1.01:
                messages.append(f"{name}: Rhat >= 1.01; chains may not have mixed.")
            if any(
                entry[key] is not None and entry[key] < 100 * len(values)
                for key in ("ess_bulk", "ess_tail")
            ):
                messages.append(f"{name}: effective sample size below 100 per chain.")
        summaries[name] = entry
    return {
        "chains": len(values),
        "retained_iterations_per_chain": retained.shape[1],
        "burnin": int(burnin),
        "parameters": summaries,
        "warnings": messages,
        "interpretation": "These diagnostics do not prove convergence or model correctness.",
    }


def save_run(output, arrays, *, settings, diagnostics=None):
    """Save numeric arrays, strict JSON metadata and a SHA-256 artifact manifest.

    Existing output paths are refused. This records caller-supplied settings;
    include seeds, units, uncertainties, data identifiers and sampler parameters.
    Load ``arrays.npz`` with ``allow_pickle=False``. No network access is required.
    """
    import compy

    output = Path(output)
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    converted = {}
    for key, array in arrays.items():
        if (
            not isinstance(key, str)
            or not key.isidentifier()
            or key in ("file", "allow_pickle")
        ):
            raise ValueError(
                "Array names must be identifiers other than file/allow_pickle"
            )
        value = np.asarray(array)
        if value.dtype.kind not in "biufc" or not np.isfinite(value).all():
            raise ValueError(f"{key}: only finite numeric arrays can be saved")
        converted[key] = value
    if not converted:
        raise ValueError("Supply at least one numeric array")
    dependencies = {}
    for package in (
        "numpy",
        "scipy",
        "matplotlib",
        "obspy",
        "tiskitpy",
        "disba",
        "arviz",
    ):
        try:
            dependencies[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            dependencies[package] = None
    code_hashes = {}
    root = Path(compy.__file__).parent
    for name in (
        "compy",
        "inv_compy",
        "compy_numerics",
        "compy_streams",
        "compy_processing",
        "Pressure_calibration",
        "ffplot",
        "compy_diagnostics",
    ):
        source = root / f"{name}.py"
        if source.is_file():
            code_hashes[source.name] = hashlib.sha256(source.read_bytes()).hexdigest()
    metadata = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "compy_version": compy.__version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "dependencies": dependencies,
        "source_sha256": code_hashes,
        "settings": settings,
        "diagnostics": diagnostics,
        "arrays": {
            key: {"shape": list(value.shape), "dtype": str(value.dtype)}
            for key, value in converted.items()
        },
    }
    # Reject unserializable settings/NaNs before creating any output.
    content = json.dumps(metadata, indent=2, allow_nan=False) + "\n"
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".compy-run-", dir=output.parent) as temporary:
        staging = Path(temporary)
        np.savez_compressed(staging / "arrays.npz", **converted)
        (staging / "metadata.json").write_text(content, encoding="utf-8")
        manifest = {
            name: hashlib.sha256((staging / name).read_bytes()).hexdigest()
            for name in ("arrays.npz", "metadata.json")
        }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        # Reserve the destination exclusively, including against concurrent runs.
        # Publish the complete record with a same-filesystem directory rename.
        output.mkdir(exist_ok=False)
        try:
            os.replace(staging, output)
        except OSError:
            output.rmdir()
            raise
    return output
