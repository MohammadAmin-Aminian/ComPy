"""Native-rate pressure spectra with one estimator from 0.02 s to about 1260 s.

Every channel is read at its own sampling rate and cut into 3-h windows with
50 % overlap (starts on a fixed 1.5-h grid from 00:00 UTC). A window is used
only when it has no gap. Its power spectral density is the Welch average of
2700-s Hann segments with 75 % overlap (the ObsPy PPSD convention for a 3-h
window), divided by the squared modulus of the full StationXML response, and
averaged as power over 1/4-octave bands centred on a 1/8-octave grid.

Two sample streams of the same window feed the one estimator:

* periods shorter than ``SPLIT_PERIOD_S`` (4 s) come from the native samples;
* periods of 4 s and longer come from the same samples decimated to 1 Hz with
  a zero-phase polyphase filter (flat to 0.4 Hz) and high-passed at 0.2 mHz
  (4-pole Butterworth, forward and backward) over a context of 6 h on each
  side of the window, which removes the tide and slower signals that would
  otherwise leak into the lowest Welch bins of absolute gauges. The high-pass
  changes the level at 1260 s by less than 1e-8 dB; it is divided out anyway.

Grid nodes are ``1060 / 2**(j/8)`` s, the nodes of the pressure atlas, extended
to 1260.6 s. Above about 470 s a 1/4-octave band is narrower than the Welch
resolution (1/2700 Hz); there the level is interpolated between the nearest
Welch bins, and the smoothing is the Welch resolution rather than 1/4 octave.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy import signal

WINDOW_S = 10800.0
STEP_S = 5400.0
SEGMENT_S = 2700.0
SEGMENT_OVERLAP = 0.75
SPLIT_PERIOD_S = 4.0
LONG_RATE_HZ = 1.0
HIGHPASS_HZ = 2.0e-4
HIGHPASS_ORDER = 4
CONTEXT_S = 21600.0
MIN_CONTEXT_S = (
    10800.0  # less context on a side: periods above LONG_CONTEXT_PERIOD_S are not kept
)
LONG_CONTEXT_PERIOD_S = 300.0
EDGE_DECIMATION_S = 120.0  # 1-Hz samples this close to a gap are not used
PRODUCT_VERSION = "native-1/4oct-v1"


def grid_periods(p_min=0.02, p_max=1300.0, anchor=1060.0, per_octave=8):
    """Grid periods [s], ascending: ``anchor / 2**(j/per_octave)`` within [p_min, p_max]."""
    j_lo = math.ceil(-per_octave * math.log2(p_max / anchor))
    j_hi = math.floor(per_octave * math.log2(anchor / p_min))
    j = np.arange(j_lo, j_hi + 1)
    return np.sort(anchor / 2.0 ** (j / per_octave))


GRID = grid_periods()


def band_weights(freqs, periods, width_octaves=0.25):
    """Weights (n_periods, n_freqs) for a power mean over 1/4-octave bands.

    A band that holds no frequency bin is interpolated linearly in log frequency
    between the two nearest bins. Returns the weights and a boolean array that
    is True where the band holds at least one bin.
    """
    from scipy import sparse

    freqs = np.asarray(freqs, float)
    fc = 1.0 / np.asarray(periods, float)
    half = 2.0 ** (width_octaves / 2.0)
    rows, cols, vals = [], [], []
    resolved = np.zeros(len(fc), bool)
    for i, f in enumerate(fc):
        lo, hi = (
            np.searchsorted(freqs, f / half, "left"),
            np.searchsorted(freqs, f * half, "left"),
        )
        if hi > lo:
            rows += [i] * (hi - lo)
            cols += list(range(lo, hi))
            vals += [1.0 / (hi - lo)] * (hi - lo)
            resolved[i] = True
            continue
        k = int(np.searchsorted(freqs, f))
        if k <= 1 or k >= len(freqs):
            continue
        a = (math.log(f) - math.log(freqs[k - 1])) / (
            math.log(freqs[k]) - math.log(freqs[k - 1])
        )
        rows += [i, i]
        cols += [k - 1, k]
        vals += [1.0 - a, a]
    w = sparse.csr_matrix((vals, (rows, cols)), shape=(len(fc), len(freqs)))
    return w, resolved


_WELCH_CACHE = {}


def welch_psd(x, fs):
    """One-sided Welch PSD of one window (2700-s Hann segments, 75 % overlap, linear detrend).

    Same result as ``scipy.signal.welch(x, fs, "hann", nperseg, noverlap, detrend="linear")``
    (tested to 1e-10), with the per-segment detrend done in closed form for speed.
    """
    from numpy.lib.stride_tricks import sliding_window_view
    from scipy import fft

    nper = int(round(SEGMENT_S * fs))
    step = nper - int(round(nper * SEGMENT_OVERLAP))
    key = (nper, fs)
    if key not in _WELCH_CACHE:
        win = signal.get_window("hann", nper)
        t = np.arange(nper) - (nper - 1) / 2.0
        _WELCH_CACHE[key] = (
            win,
            t,
            float(t @ t),
            1.0 / (fs * float(win @ win)),
            np.fft.rfftfreq(nper, 1.0 / fs),
        )
    win, t, tt, scale, freqs = _WELCH_CACHE[key]
    seg = sliding_window_view(np.asarray(x, float), nper)[::step]
    mean = seg.mean(axis=1, keepdims=True)
    slope = (seg @ t)[:, None] / tt
    spec = fft.rfft((seg - mean - slope * t[None, :]) * win[None, :], axis=1)
    p = np.mean(spec.real**2 + spec.imag**2, axis=0) * scale
    p[1:] *= 2.0
    if nper % 2 == 0:
        p[-1] /= 2.0
    return freqs, p


def decimation_factors(fs, target=LONG_RATE_HZ):
    """Integer up/down factors that take ``fs`` to ``target`` exactly."""
    r = Fraction(target).limit_denominator(1000) / Fraction(fs).limit_denominator(
        100000
    )
    return r.numerator, r.denominator


def highpass_sos(fs=LONG_RATE_HZ):
    return signal.butter(
        HIGHPASS_ORDER, HIGHPASS_HZ, btype="highpass", fs=fs, output="sos"
    )


@dataclass
class ResponseTable:
    """|response|^2 of one channel at the native and 1-Hz Welch frequencies."""

    native_freqs: np.ndarray
    native_power: np.ndarray
    long_freqs: np.ndarray
    long_power: np.ndarray
    sensitivity: float
    input_units: str
    response_db_grid: np.ndarray  # 10 log10 |R|^2 / sensitivity^2 at the grid periods


def response_table(response, fs):
    """Evaluate a full ObsPy ``Response`` (output in its input units) for one channel."""
    nper_n = int(round(SEGMENT_S * fs))
    fn = np.fft.rfftfreq(nper_n, 1.0 / fs)
    fl = np.fft.rfftfreq(int(round(SEGMENT_S * LONG_RATE_HZ)), 1.0 / LONG_RATE_HZ)

    def power(f):
        out = np.full(len(f), np.nan)
        pos = f > 0
        r = response.get_evalresp_response_for_frequencies(f[pos], output="DEF")
        out[pos] = np.abs(r) ** 2
        return out

    sens = response.instrument_sensitivity
    s = float(sens.value)
    units = str(sens.input_units or "")
    rg = power(1.0 / GRID)
    return ResponseTable(
        fn, power(fn), fl, power(fl), s, units, 10 * np.log10(rg / s**2)
    )


@dataclass
class ChannelSpectra:
    """Per-window levels of one channel on ``GRID`` [dB re 1 unit^2/Hz]."""

    seed_id: str
    fs: float
    starts: list = field(default_factory=list)
    levels: list = field(default_factory=list)
    maxabs: list = field(default_factory=list)
    std: list = field(default_factory=list)
    long_context_ok: list = field(default_factory=list)

    def save(self, path, response: ResponseTable, meta: dict):
        lv = np.array(self.levels, np.float32).reshape(-1, len(GRID))
        np.savez_compressed(
            path, periods=GRID, starts=np.array(self.starts, float), levels_db=lv,
            maxabs_counts=np.array(self.maxabs, float), std_counts=np.array(self.std, float),
            long_context_ok=np.array(self.long_context_ok, bool), response_db=response.response_db_grid,
            sensitivity=response.sensitivity, input_units=response.input_units, fs=self.fs,
            seed_id=self.seed_id, product=PRODUCT_VERSION, meta=json.dumps(meta),
        )  # fmt: skip


class _Estimator:
    """Holds the per-channel tables so that each window costs two Welch calls."""

    def __init__(self, response_tab: ResponseTable, fs):
        from scipy import sparse

        self.fs = fs
        self.r = response_tab
        short = GRID < SPLIT_PERIOD_S
        self.short_idx = np.where(short)[0]
        self.long_idx = np.where(~short)[0]
        self.up, self.down = decimation_factors(fs)
        self.sos = highpass_sos()
        wn, _ = band_weights(response_tab.native_freqs, GRID[short])
        wl, self.long_resolved = band_weights(response_tab.long_freqs, GRID[~short])
        _, hp = signal.sosfreqz(self.sos, worN=response_tab.long_freqs, fs=LONG_RATE_HZ)
        with np.errstate(divide="ignore", invalid="ignore"):
            inv_n = np.where(
                response_tab.native_power > 0, 1.0 / response_tab.native_power, 0.0
            )
            inv_l = np.where(response_tab.long_power * np.abs(hp) ** 2 > 0,
                             1.0 / (response_tab.long_power * np.abs(hp) ** 2), 0.0)  # fmt: skip
        # Fold the response into the band weights; bands above the native Nyquist stay NaN.
        self.wn = (wn @ sparse.diags(inv_n)).tocsr()
        self.wn_valid = (np.asarray(wn.sum(axis=1)).ravel() > 0) & (
            1.0 / GRID[short] < 0.45 * fs
        )
        self.wl = (wl @ sparse.diags(inv_l)).tocsr()
        self.wl_valid = (np.asarray(wl.sum(axis=1)).ravel() > 0) & (
            1.0 / GRID[~short] < 0.4 * LONG_RATE_HZ
        )

    def levels(self, native, long_ctx, long_off, long_ctx_ok):
        out = np.full(len(GRID), np.nan)
        if len(self.short_idx) and self.wn_valid.any():
            _, p = welch_psd(native, self.fs)
            v = self.wn @ p
            out[self.short_idx[self.wn_valid]] = 10 * np.log10(
                np.maximum(v[self.wn_valid], 1e-300)
            )
        x = long_ctx - np.polyval(
            np.polyfit(np.arange(len(long_ctx)), long_ctx, 1), np.arange(len(long_ctx))
        )
        y = signal.sosfiltfilt(
            self.sos, x, padtype="odd", padlen=min(len(x) - 1, 3 * int(1 / HIGHPASS_HZ))
        )
        n = int(round(WINDOW_S * LONG_RATE_HZ))
        _, p = welch_psd(y[long_off : long_off + n], LONG_RATE_HZ)
        v = self.wl @ p
        lv = np.full(len(self.long_idx), np.nan)
        lv[self.wl_valid] = 10 * np.log10(np.maximum(v[self.wl_valid], 1e-300))
        if not long_ctx_ok:
            lv[GRID[self.long_idx] > LONG_CONTEXT_PERIOD_S] = np.nan
        out[self.long_idx] = lv
        return out


_DOY = re.compile(r"(\d{4})\.(\d{3})(?:\.mseed)?$")


def file_day(path):
    """UTC day number (days since 1970-01-01) from an SDS file name."""
    m = _DOY.search(str(path))
    if not m:
        raise ValueError(f"not an SDS day file: {path}")
    y, d = int(m.group(1)), int(m.group(2))
    return (
        int((np.datetime64(f"{y}-01-01") - np.datetime64("1970-01-01")).astype(int))
        + d
        - 1
    )


def _one_location(stream, loc):
    """Traces of the channel's own location code.

    Data and metadata often disagree on the code of a single sensor, so a stream
    with one location code is taken as it is. With several codes only the
    requested one is kept: two sensors under one channel code must not be mixed.
    """
    codes = {tr.stats.location for tr in stream}
    if len(codes) <= 1:
        return stream
    own = stream.select(location=loc)
    if not len(own):
        raise ValueError(
            f"location codes {sorted(codes)} in the data and none is '{loc}'"
        )
    return own


def _contiguous(stream, fs):
    """Gap-free runs as (t0 [s since epoch], float64 samples), merged across files, and the number of traces left out.

    Traces whose sampling rate is not ``fs`` (beyond a clock-drift tolerance of
    2e-5) are left out and counted: with another rate the window length and the
    frequency axis would be wrong without any error.
    """
    other = 0
    for tr in list(stream):
        if abs(tr.stats.sampling_rate - fs) / fs < 2e-5:
            tr.stats.sampling_rate = fs
        else:
            stream.remove(tr)
            other += 1
    stream.merge(method=1, fill_value=None)
    runs = []
    for tr in stream.split():
        if tr.stats.npts:
            runs.append(
                (float(tr.stats.starttime.timestamp), np.asarray(tr.data, np.float64))
            )
    runs.sort(key=lambda r: r[0])
    return runs, other


def process_channel(
    files,
    response,
    seed_id,
    fs,
    out_path,
    chunk_days=None,
    read=None,
    meta=None,
    log=None,
):
    """Compute the native-rate levels of one channel from its day files and save them.

    Args:
        files: SDS day files of the channel (any order).
        response: ObsPy ``Response`` of the channel (one epoch).
        seed_id: NET.STA.LOC.CHA kept from the files.
        fs: nominal sampling rate [Hz].
        out_path: ``.npz`` written at the end.
        chunk_days: days computed per read (default from the rate); one more day is read on each side.
        read: function path -> ObsPy Stream (default ``obspy.read``).
    """
    from obspy import Stream

    if read is None:
        from obspy import read as read
    net, sta, loc, cha = seed_id.split(".")
    if (
        chunk_days is None
    ):  # about 2e7 samples per chunk at most, plus one day of context on each side
        chunk_days = 4 if fs <= 50 else (2 if fs <= 125 else 1)
    rtab = response_table(response, fs)
    est = _Estimator(rtab, fs)
    by_day = {}
    for f in files:
        by_day.setdefault(file_day(f), []).append(f)
    days = sorted(by_day)
    spec = ChannelSpectra(seed_id, fs)
    n_long = int(round(WINDOW_S * LONG_RATE_HZ))
    ctx = int(round(CONTEXT_S * LONG_RATE_HZ))
    stats = dict(
        windows_tried=0, windows_gap=0, chunks=0, read_errors=0, traces_other_rate=0
    )
    i = 0
    while i < len(days):
        block = [d for d in days[i:] if d < days[i] + chunk_days]
        i += len(block)
        d0, d1 = block[0], block[-1] + 1
        st = Stream()
        for d in range(d0 - 1, d1 + 1):
            for f in by_day.get(d, []):
                try:
                    s = read(str(f))
                except Exception:
                    stats["read_errors"] += 1
                    continue
                st += s.select(network=net, station=sta, channel=cha)
        if not len(st):
            continue
        st = _one_location(st, loc)
        runs, other = _contiguous(st, fs)
        stats["traces_other_rate"] += other
        if not runs:
            continue
        stats["chunks"] += 1
        for t0, x in runs:
            t1 = t0 + len(x) / fs
            # 1-Hz copy of this run, on whole seconds.
            off = (math.ceil(t0) - t0) * fs
            k0 = int(round(off))
            mean = float(x.mean())
            x -= mean  # in place; Welch detrends each segment, and the clip test adds the mean back
            if est.up == 1 and est.down == 1:
                lx, lt0 = x, t0
            else:
                lx = signal.resample_poly(
                    x[k0:], est.up, est.down, window=("kaiser", 8.0), padtype="line"
                )
                lt0 = t0 + k0 / fs
            edge = int(EDGE_DECIMATION_S * LONG_RATE_HZ)
            l_lo, l_hi = lt0 + edge, lt0 + (len(lx) - edge) / LONG_RATE_HZ
            first = math.ceil((max(t0, d0 * 86400.0)) / STEP_S) * STEP_S
            ws = first
            while ws < d1 * 86400.0 and ws + WINDOW_S <= t1 + 1e-6:
                stats["windows_tried"] += 1
                a = int(round((ws - t0) * fs))
                b = a + int(round(WINDOW_S * fs))
                if a < 0 or b > len(x) or ws < l_lo or ws + WINDOW_S > l_hi:
                    stats["windows_gap"] += 1
                    ws += STEP_S
                    continue
                la = int(round((ws - lt0) * LONG_RATE_HZ))
                ca = max(int(round((l_lo - lt0) * LONG_RATE_HZ)), la - ctx)
                cb = min(int(round((l_hi - lt0) * LONG_RATE_HZ)), la + n_long + ctx)
                ok = (la - ca) >= MIN_CONTEXT_S * LONG_RATE_HZ and (
                    cb - la - n_long
                ) >= MIN_CONTEXT_S * LONG_RATE_HZ
                w = x[a:b]
                lev = est.levels(w, lx[ca:cb], la - ca, ok)
                spec.starts.append(ws)
                spec.levels.append(lev.astype(np.float32))
                spec.maxabs.append(float(max(abs(w.max() + mean), abs(w.min() + mean))))
                spec.std.append(float(np.std(w)))
                spec.long_context_ok.append(bool(ok))
                ws += STEP_S
        if log:
            log(f"{seed_id} days {d0}-{d1 - 1}: {len(spec.starts)} windows")
    meta = dict(meta or {}, **stats, files=len(files), days=len(days), product=PRODUCT_VERSION,
                window_s=WINDOW_S, step_s=STEP_S, segment_s=SEGMENT_S, split_period_s=SPLIT_PERIOD_S,
                highpass_hz=HIGHPASS_HZ, context_s=CONTEXT_S)  # fmt: skip
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(out_path) + ".part.npz")
    spec.save(tmp, rtab, meta)
    tmp.replace(out_path)
    return meta


def pick_epoch_windows(starts, limits):
    """Which epoch each window comes from, when the day files at an epoch limit were computed with both responses.

    Args:
        starts: one array of window start times [s since epoch] per response epoch.
        limits: (start, end) of each epoch [s since epoch]; ``inf`` for an open end.

    Returns:
        Sorted list of (start, epoch index, index inside that epoch). A window is
        taken from the epoch that contains its start; a window outside every
        epoch is taken from the first epoch that computed it.
    """
    best = {}
    for k, st in enumerate(starts):
        for i, s in enumerate(st):
            rank = (0 if limits[k][0] <= s < limits[k][1] else 1, k)
            if s not in best or rank < best[s][0]:
                best[s] = (rank, k, i)
    return sorted((s, k, i) for s, (_, k, i) in best.items())


def load(path):
    """Read a native-rate product file: dict with periods, starts, levels_db and metadata."""
    z = np.load(path, allow_pickle=False)
    out = {k: z[k] for k in z.files}
    out["meta"] = json.loads(str(out["meta"]))
    for k in ("seed_id", "input_units", "product"):
        out[k] = str(out[k])
    out["fs"] = float(out["fs"])
    out["sensitivity"] = np.ravel(out["sensitivity"]).astype(float)
    out["response_db"] = np.atleast_2d(out["response_db"])
    return out
