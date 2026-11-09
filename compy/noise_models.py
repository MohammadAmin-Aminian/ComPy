"""Reference noise levels for OBS diagnostics, without ObsPy.

* Peterson (1993) new low and high noise models (NLNM, NHNM) for vertical
  acceleration, in dB re 1 (m/s²)²/Hz. The tables are Peterson's piecewise
  ``A + B log10(T)`` coefficients (USGS Open-File Report 93-322, Table 1),
  the same values ObsPy's ``get_nlnm``/``get_nhnm`` evaluate.
* Webb (1998, Fig. 2) seafloor pressure bounds, in dB re 1 Pa²/Hz, read
  from ``compy_pressure/reference_data/pressure_references.json``.

Land-based Peterson models are only a yardstick for OBS: seafloor verticals
sit above the NHNM at infragravity frequencies whenever compliance or tilt
noise is present (Webb 1998; Janiszewski et al. 2023, GJI 233, 297–315).
"""

import json
from pathlib import Path

import numpy as np

#: Peterson (1993) NLNM: period lower bound [s], A, B (last row ends at 1e5 s).
NLNM_TABLE = np.array([
    [0.10, -162.36, 5.64], [0.17, -166.70, 0.00], [0.40, -170.00, -8.30],
    [0.80, -166.40, 28.90], [1.24, -168.60, 52.48], [2.40, -159.98, 29.81],
    [4.30, -141.10, 0.00], [5.00, -71.36, -99.77], [6.00, -97.26, -66.49],
    [10.00, -132.18, -31.57], [12.00, -205.27, 36.16], [15.60, -37.65, -104.33],
    [21.90, -114.37, -47.10], [31.60, -160.58, -16.28], [45.00, -187.50, 0.00],
    [70.00, -216.47, 15.70], [101.00, -185.00, 0.00], [154.00, -168.34, -7.61],
    [328.00, -217.43, 11.90], [600.00, -258.28, 26.60], [10000.00, -346.88, 48.75],
])  # fmt: skip

#: Peterson (1993) NHNM: period lower bound [s], A, B (last row ends at 1e5 s).
NHNM_TABLE = np.array([
    [0.10, -108.73, -17.23], [0.22, -150.34, -80.50], [0.32, -122.31, -23.87],
    [0.80, -116.85, 32.51], [3.80, -108.48, 18.08], [4.60, -74.66, -32.95],
    [6.30, 0.66, -127.18], [7.90, -93.37, -22.42], [15.40, 73.54, -162.98],
    [20.00, -151.52, 10.01], [354.80, -206.66, 31.63],
])  # fmt: skip

_PERIOD_MAX = 1e5


def _peterson(table, freqs):
    f = np.asarray(freqs, float)
    out = np.full(f.shape, np.nan)
    with np.errstate(divide="ignore"):
        period = 1.0 / f
    ok = (period >= table[0, 0]) & (period <= _PERIOD_MAX)
    row = np.searchsorted(table[:, 0], period[ok], side="right") - 1
    out[ok] = table[row, 1] + table[row, 2] * np.log10(period[ok])
    return out


def nlnm(freqs):
    """Peterson NLNM acceleration PSD [dB re 1 (m/s²)²/Hz]; NaN outside 0.1–1e5 s."""
    return _peterson(NLNM_TABLE, freqs)


def nhnm(freqs):
    """Peterson NHNM acceleration PSD [dB re 1 (m/s²)²/Hz]; NaN outside 0.1–1e5 s."""
    return _peterson(NHNM_TABLE, freqs)


def displacement_to_acceleration_db(freqs, psd_displacement):
    """Displacement PSD [m²/Hz] to acceleration PSD [dB re 1 (m/s²)²/Hz]."""
    f = np.asarray(freqs, float)
    acc = np.asarray(psd_displacement, float) * (2 * np.pi * f) ** 4
    with np.errstate(divide="ignore"):
        return 10 * np.log10(acc)


def _reference_file():
    return (
        Path(__file__).resolve().parents[1]
        / "compy_pressure"
        / "reference_data"
        / "pressure_references.json"
    )


def webb_pressure_bounds(path=None):
    """Webb (1998) seafloor pressure low/high bounds.

    Args:
        path (path or None): ``pressure_references.json``; defaults to the
            copy shipped with ``compy_pressure``.

    Returns:
        tuple: (frequency [Hz], low [dB re 1 Pa²/Hz], high [dB re 1 Pa²/Hz]),
        sorted by frequency, or None if the file or entry is missing.
    """
    path = Path(path) if path is not None else _reference_file()
    if not path.exists():
        return None
    entry = json.loads(path.read_text()).get("webb")
    if not entry:
        return None
    points = entry.get("points") or entry.get("data") or entry.get("curve")
    if not points:
        return None
    rows = np.array(
        [
            (p["frequency_hz"], p["low_db_pa2_hz"], p["high_db_pa2_hz"])
            if isinstance(p, dict)
            else tuple(p[:3])
            for p in points
        ],
        float,
    )
    rows = rows[np.argsort(rows[:, 0])]
    return rows[:, 0], rows[:, 1], rows[:, 2]
