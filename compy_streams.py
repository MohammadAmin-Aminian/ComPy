"""Sample-aligned, non-mutating ObsPy window helpers."""

import numpy as np
from obspy import Stream

from compy_numerics import positive


def validate_stream(stream):
    if not isinstance(stream, Stream) or not stream:
        raise ValueError("a nonempty ObsPy Stream is required")
    fs = positive(stream[0].stats.sampling_rate, "sampling rate")
    ids = [tr.id for tr in stream]
    if len(set(ids)) != len(ids):
        raise ValueError("merge trace fragments and resolve gaps before processing")
    origin = stream[0].stats.starttime
    for tr in stream:
        if tr.stats.sampling_rate != fs:
            raise ValueError("all traces must have the same sampling rate")
        offset = (tr.stats.starttime - origin) * fs
        if abs(offset - round(offset)) > 1e-5:
            raise ValueError("trace start times must lie on a common sample grid")
        if (
            not len(tr)
            or np.ma.isMaskedArray(tr.data)
            or not np.all(np.isfinite(tr.data))
        ):
            raise ValueError("traces must contain finite, unmasked samples")
    return fs


def sample_count(seconds, fs, name):
    count = seconds * fs
    if not np.isfinite(count) or abs(count - round(count)) > 1e-5 or count < 1:
        raise ValueError(f"{name} must correspond to a positive integer sample count")
    return int(round(count))


def cut_stream_with_overlap(stream, window_length, overlap):
    """Return copied complete windows, using half-open sample intervals.

    Duration N/fs contains N samples; ObsPy's inclusive endpoint is adjusted
    by one sample. Incomplete final windows are omitted.
    """
    fs = validate_stream(stream)
    window_length = positive(window_length, "window_length")
    overlap = float(overlap)
    if not np.isfinite(overlap) or not 0 <= overlap < window_length:
        raise ValueError("overlap must satisfy 0 <= overlap < window_length")
    width = sample_count(window_length, fs, "window_length")
    step = sample_count(window_length - overlap, fs, "step")
    start = max(tr.stats.starttime for tr in stream)
    end = min(tr.stats.endtime for tr in stream)
    n = int(round((end - start) * fs)) + 1
    if n < width:
        return []
    return [
        stream.slice(
            start + i / fs, start + (i + width - 1) / fs, nearest_sample=False
        ).copy()
        for i in range(0, n - width + 1, step)
    ]


def split_stream(stream, duration):
    return cut_stream_with_overlap(stream, duration, 0)


def trim_streams_to_same_length(stream, channels=("BH1", "BH2", "BDH", "BHZ")):
    traces = []
    for channel in channels:
        matches = stream.select(channel=channel)
        if len(matches) != 1:
            raise ValueError(
                f"expected exactly one trace for {channel}, found {len(matches)}"
            )
        traces.append(matches[0].copy())
    result = Stream(traces)
    validate_stream(result)
    start = max(tr.stats.starttime for tr in result)
    end = min(tr.stats.endtime for tr in result)
    if end < start:
        raise ValueError("channels have no common time interval")
    return result.trim(start, end, nearest_sample=False)
