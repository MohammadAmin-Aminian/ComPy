"""Read the stage C and E products of processed stations without ObsPy.

The pipeline writes, per station, ``<STA>_compliance_<gain>.npz`` (stacked
spectra, compliance, accepted window starts), ``<STA>_summary.json`` (gates,
monthly stacks) and, optionally, ``<STA>_window_levels.npy`` (per-window
transfer-function and PSD levels) and ``inversion/<STA>_<variant>.json``.
:func:`load_station` gathers them in one :class:`StationProducts`, and
:func:`discover` lists the stations of a results folder, so that the
diagnostics (``compy.anomaly``, ``compy.plotting.diagnostics``) can run on
archived results, e.g. ``docs/results/*/results``.

All quantities are SI: Hz, Pa, m, 1/Pa; times are POSIX seconds (UTC).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

#: Order in which gain labels are preferred, as in ``pipeline/e_invert.py``.
GAIN_LABELS = ("compy", "bruitfm", "nominal")


@dataclass
class StationProducts:
    """Archived compliance products of one station.

    Attributes:
        station (str): Station code.
        label (str): Gain label of the result ("compy", "bruitfm", "nominal").
        depth (float): Water depth [m].
        gain (float): Pressure gain applied.
        band (tuple): Compliance band (f_min, f_max) [Hz].
        freqs (ndarray): Output frequencies [Hz].
        compliance (ndarray): Normalized compliance [1/Pa].
        uncertainty (ndarray): Bendat & Piersol one-sigma error [1/Pa].
        coherence (ndarray): Amplitude coherence of the stack.
        phase_deg (ndarray): Phase of Z/P [deg].
        window_mean, window_std (ndarray): Per-window compliance statistics.
        all_freqs (ndarray): Full frequency axis [Hz].
        S_zz, S_pp (ndarray): Stacked displacement [m²/Hz] and calibrated
            pressure [Pa²/Hz] PSDs on ``all_freqs``.
        S_pz (ndarray): Stacked cross-spectrum E[conj(P) Z] [m Pa/Hz].
        mean_coherence_all (ndarray): Mean per-window coherence.
        window_starts (ndarray): POSIX start times of accepted windows [s].
        summary (dict): The stage C summary of this gain label.
        monthly (dict): ``{"YYYY-MM": {"n", "compliance", "uncertainty"}}``.
        window_levels (ndarray or None): Rows (time, log10 |S_pz|/S_pp,
            pressure dB, acceleration dB) of accepted windows, if saved.
        inversion (dict): Stage E JSON records by variant name.
        n_windows_total (int or None): Candidate windows before gating.
        window_s (float or None): Window length [s].
    """

    station: str
    label: str
    depth: float
    gain: float
    band: tuple
    freqs: np.ndarray
    compliance: np.ndarray
    uncertainty: np.ndarray
    coherence: np.ndarray
    phase_deg: np.ndarray
    window_mean: np.ndarray
    window_std: np.ndarray
    all_freqs: np.ndarray
    S_zz: np.ndarray
    S_pp: np.ndarray
    S_pz: np.ndarray
    mean_coherence_all: np.ndarray
    window_starts: np.ndarray
    summary: dict = field(default_factory=dict)
    monthly: dict = field(default_factory=dict)
    window_levels: np.ndarray | None = None
    inversion: dict = field(default_factory=dict)
    n_windows_total: int | None = None
    window_s: float | None = None

    @property
    def in_band(self):
        """Boolean mask of ``freqs`` inside ``band``."""
        lo, hi = self.band
        return (self.freqs >= lo - 1e-12) & (self.freqs <= hi + 1e-12)

    @property
    def acceptance(self):
        """Accepted over candidate windows, or None if unknown."""
        if not self.n_windows_total:
            return None
        return len(self.window_starts) / self.n_windows_total

    def monthly_matrix(self, min_windows=1):
        """Months, window counts and an (n_months, n_band) compliance array.

        Only months whose stack covers the whole band and has at least
        ``min_windows`` windows are returned.
        """
        n_band = int(self.in_band.sum())
        rows = [
            (key, m["n"], m["compliance"])
            for key, m in sorted(self.monthly.items())
            if len(m["compliance"]) == n_band and m["n"] >= min_windows
        ]
        if not rows:
            return [], np.array([], int), np.empty((0, n_band))
        keys, counts, values = zip(*rows)
        return list(keys), np.array(counts), np.array(values, float)

    def sediment_thickness(self, variant="final"):
        """Posterior median sediment thickness [m] of an inversion variant."""
        record = self.inversion.get(variant)
        if not record:
            return None
        par = record.get("parameters", {}).get("sediment_thickness")
        return None if par is None else float(par["median"])


def month_start(key):
    """POSIX time of the first instant of a ``YYYY-MM`` month (UTC)."""
    year, month = map(int, key.split("-"))
    return datetime(year, month, 1, tzinfo=timezone.utc).timestamp()


def _choose_label(directory, station, label):
    if label is not None:
        return label
    for candidate in GAIN_LABELS:
        if (directory / f"{station}_compliance_{candidate}.npz").exists():
            return candidate
    raise FileNotFoundError(f"no compliance file for {station} in {directory}")


def load_station(directory, station=None, label=None, *, inversion_dirs=()):
    """Load the archived products of one station.

    Args:
        directory (path): Folder with ``<STA>_compliance_<label>.npz`` and
            ``<STA>_summary.json`` (e.g. ``results/<STA>``).
        station (str or None): Station code; defaults to the folder name.
        label (str or None): Gain label; defaults to the first available
            in :data:`GAIN_LABELS`.
        inversion_dirs (iterable of path): Extra folders searched for
            ``<STA>_<variant>.json`` stage E records, besides
            ``directory/inversion``.

    Returns:
        StationProducts: The station's products.
    """
    directory = Path(directory)
    station = station or directory.name
    label = _choose_label(directory, station, label)
    with np.load(directory / f"{station}_compliance_{label}.npz") as d:
        arrays = {k: d[k] for k in d.files}
    summary_path = directory / f"{station}_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    entry = summary.get("results", {}).get(label, {})
    levels_path = directory / f"{station}_window_levels.npy"
    levels = np.load(levels_path) if levels_path.exists() else None
    inversion = {}
    for folder in (directory / "inversion", *map(Path, inversion_dirs)):
        for path in sorted(folder.glob(f"{station}_*.json")) if folder.is_dir() else []:
            variant = path.stem[len(station) + 1 :]
            try:
                inversion[variant] = json.loads(path.read_text())
            except json.JSONDecodeError:
                continue
    return StationProducts(
        station=station,
        label=label,
        depth=float(arrays["depth"]),
        gain=float(arrays["gain"]),
        band=tuple(float(v) for v in arrays["band"]),
        freqs=arrays["freqs"],
        compliance=arrays["compliance"],
        uncertainty=arrays["uncertainty"],
        coherence=arrays["coherence"],
        phase_deg=arrays["phase_deg"],
        window_mean=arrays["window_mean"],
        window_std=arrays["window_std"],
        all_freqs=arrays["all_freqs"],
        S_zz=arrays["S_zz"],
        S_pp=arrays["S_pp"],
        S_pz=arrays["S_pz"],
        mean_coherence_all=arrays["mean_coherence_all"],
        window_starts=np.sort(np.asarray(arrays["window_starts"], float)),
        summary=summary,
        monthly=entry.get("monthly", {}),
        window_levels=levels,
        inversion=inversion,
        n_windows_total=summary.get("windows_total"),
        window_s=summary.get("window_s"),
    )


def discover(results_root, label=None, *, inversion_root=None):
    """Load every station folder of a stage C results folder.

    Args:
        results_root (path): Folder holding one sub-folder per station.
        label (str or None): Gain label (see :func:`load_station`).
        inversion_root (path or None): Folder holding one sub-folder per
            station with stage E JSON records, if not under each station.

    Returns:
        dict: ``{station: StationProducts}`` sorted by station code.
    """
    out = {}
    for folder in sorted(Path(results_root).iterdir()):
        if not folder.is_dir() or not list(
            folder.glob(f"{folder.name}_compliance_*.npz")
        ):
            continue
        extra = [Path(inversion_root) / folder.name] if inversion_root else []
        try:
            out[folder.name] = load_station(folder, label=label, inversion_dirs=extra)
        except FileNotFoundError:
            continue
    return out
