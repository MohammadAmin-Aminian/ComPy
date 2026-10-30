"""Small helpers shared by the spectra, calibration and pilot modules."""

import csv
from pathlib import Path

import numpy as np

# Home of the default exports/ and results/ folders, and the base of relative
# paths in a master table. An installed package must not write next to its
# code, so this is the current directory.
PROJECT_DIRECTORY = Path.cwd()


def complete_samples(stream, start, sample_rate, npts):
    """Accept one continuous, finite window; never interpolate gaps or pad data."""
    delta = 1 / sample_rate
    for trace in stream:
        if trace.stats.sampling_rate != sample_rate:
            raise ValueError("Waveform sample rate differs from inventory")
        offset = max(
            0, int(np.ceil((start - trace.stats.starttime) * sample_rate - 1e-5))
        )
        actual = trace.stats.starttime + offset * delta
        if actual < start - 1e-6 or actual - start > delta + 1e-6:
            continue
        data = trace.data[offset : offset + npts]
        if len(data) != npts or np.ma.is_masked(data) or not np.isfinite(data).all():
            continue
        if np.ptp(data) == 0:
            continue
        return np.asarray(data, dtype=float), actual
    return None, None


def write_csv(path, rows, *, fieldnames=None):
    if rows or fieldnames:
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=list(rows[0]) if rows else fieldnames,
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
