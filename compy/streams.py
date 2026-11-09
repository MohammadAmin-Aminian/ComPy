"""Sample-aligned, non-mutating ObsPy window helpers.

ObsPy is imported when a function needs it, so that the array-level parts
of ``compy.compliance`` (estimators, gates, corrections) import without it.
"""

import numpy as np

from compy.numerics import positive


def validate_stream(stream):
    """Check a stream is ready for spectral work and return its sampling rate [Hz].

    Requires a nonempty ObsPy Stream with unique trace IDs, one sampling rate
    and start times on a common sample grid.
    """
    from obspy import Stream

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
    """Convert a duration ``seconds`` at ``fs`` Hz to a whole, positive number of samples."""
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
    """Cut a stream into consecutive, non-overlapping windows of ``duration`` seconds."""
    return cut_stream_with_overlap(stream, duration, 0)


def trim_streams_to_same_length(stream, channels=("BH1", "BH2", "BDH", "BHZ")):
    """Select one trace per channel in ``channels`` and trim all to their common time span."""
    from obspy import Stream

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


def split_and_save_stream(stream, interval_minutes, output_dir):
    """Save consecutive copied windows, including the final partial window."""
    from pathlib import Path

    fs = validate_stream(stream)
    width = sample_count(
        positive(interval_minutes, "interval_minutes") * 60, fs, "interval"
    )
    start = max(tr.stats.starttime for tr in stream)
    end = min(tr.stats.endtime for tr in stream)
    if end < start:
        raise ValueError("channels have no common time interval")
    n = int(round((end - start) * fs)) + 1
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    paths = []
    for i in range(0, n, width):
        first = start + i / fs
        last = start + (min(i + width, n) - 1) / fs
        window = stream.slice(first, last, nearest_sample=False).copy()
        filename = (
            f"{stream[0].stats.station}_{first.strftime('%Y%m%dT%H%M%S%f')}.mseed"
        )
        path = output / filename
        if path.exists():
            raise FileExistsError(path)
        window.write(str(path), format="MSEED")
        paths.append(path)
    return paths
