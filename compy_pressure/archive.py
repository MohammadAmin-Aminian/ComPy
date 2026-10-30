"""Read daily SDS MiniSEED files from a full archive or a station folder."""

from dataclasses import dataclass
from datetime import date
from pathlib import Path

from obspy import Stream, UTCDateTime, read
from obspy.io.mseed.util import get_record_information

# (filename, header) location codes accepted as one channel, as reviewed in
# docs/coverage_sources.json: some archives name empty-location files "00".
LOCATION_ALIASES = {("00", "")}
SELECTION_MARGIN_SECONDS = 1.0


def alias_seed_id(seed_id):
    """Return the header identifier a filename alias stands for, or None."""
    parts = seed_id.split(".")
    for name, header in LOCATION_ALIASES:
        if len(parts) == 4 and parts[2] == name:
            return ".".join([*parts[:2], header, parts[3]])
    return None


def same_channel(filename_seed, header_seed):
    return header_seed in {filename_seed, alias_seed_id(filename_seed)}


class SampleGridError(ValueError):
    """A channel contains records on different native sample grids."""

    def __init__(self, seed_id, residual_seconds, sample_rate):
        self.seed_id = seed_id
        self.residual_seconds = residual_seconds
        self.sample_rate = sample_rate
        super().__init__(
            f"Within-channel sample-grid offset: {seed_id}; "
            f"{residual_seconds:.9g} s; merge refused"
        )


@dataclass(frozen=True)
class WaveformFile:
    path: Path
    seed_id: str
    start: UTCDateTime
    end: UTCDateTime
    sample_rate: float


def waveform_seed_id(path):
    """Recognize SDS or daily NET.STA.LOC.CHA.YYYYMMDD.mseed names."""
    parts = Path(path).name.split(".")
    if (
        len(parts) == 7
        and parts[4] == "D"
        and parts[5].isdigit()
        and parts[6].isdigit()
    ):
        return ".".join(parts[:4])
    if (
        len(parts) == 6
        and parts[5].lower() == "mseed"
        and len(parts[4]) == 8
        and parts[4].isdigit()
    ):
        stamp = parts[4]
        try:
            date.fromisoformat(f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:]}")
        except ValueError:
            return None
        parts[2] = "" if parts[2] == "--" else parts[2]
        return ".".join(parts[:4])
    return None


def discover_waveforms(directory, seed_ids=None):
    """Index actual trace headers; never load waveform samples during discovery.

    Accept SDS and downloaded daily MiniSEED filenames anywhere below directory.
    Header times, rather than filename dates, handle files crossing midnight.
    """
    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        raise FileNotFoundError(f"Waveform directory does not exist: {directory}")
    paths = []
    for path in sorted(directory.rglob("*")):
        seed_id = waveform_seed_id(path)
        if not path.is_file() or seed_id is None:
            continue
        if seed_ids is not None and seed_id not in seed_ids:
            continue
        paths.append(path)
    return index_waveform_files(paths)


def _record_span(path, seed_id):
    """Index a file from its first and last records only, or return None.

    Daily files hold fixed-length records of one channel, so these two headers
    give the same identifier, rate, start and last-sample time as reading every
    header. Gaps inside the file stay in the data and are handled when read.
    Irregular files return None and are indexed from all headers instead.
    """
    try:
        first = get_record_information(str(path))
        length, count = first["record_length"], first["number_of_records"]
        if first["excess_bytes"] or count * length != first["filesize"]:
            return None
        last = (
            first
            if count == 1
            else get_record_information(str(path), offset=(count - 1) * length)
        )
    except Exception:
        return None
    ids = {
        ".".join((info["network"], info["station"], info["location"], info["channel"]))
        for info in (first, last)
    }
    rates = {float(info["samp_rate"]) for info in (first, last)}
    if (
        len(ids) != 1
        or len(rates) != 1
        or last["record_length"] != length
        or last["endtime"] < first["starttime"]
    ):
        return None
    (header_id,) = ids
    if not same_channel(seed_id, header_id):
        raise ValueError(f"Filename/header SEED ID mismatch: {path}")
    # Drifting clocks can leave the true last sample a few ms from the last
    # record's end; read_window selects files with a margin for this.
    end = max(first["endtime"], last["endtime"])
    return WaveformFile(path, header_id, first["starttime"], end, rates.pop())


def index_waveform_files(paths):
    """Index a preselected channel's files without rescanning a large archive."""
    records = []
    for path in sorted(map(Path, paths)):
        seed_id = waveform_seed_id(path)
        if seed_id is None:
            raise ValueError(f"Unrecognized daily MiniSEED filename: {path}")
        # An empty file (e.g. a failed download) holds no samples to index.
        if path.stat().st_size == 0:
            continue
        span = _record_span(path, seed_id)
        if span is not None:
            records.append(span)
            continue
        try:
            headers = read(str(path), format="MSEED", headonly=True)
        except Exception as exc:
            raise RuntimeError(f"Cannot read MiniSEED header: {path}") from exc
        for trace in headers:
            if not same_channel(seed_id, trace.id):
                raise ValueError(f"Filename/header SEED ID mismatch: {path}")
            records.append(
                WaveformFile(
                    path,
                    trace.id,
                    trace.stats.starttime,
                    trace.stats.endtime,
                    float(trace.stats.sampling_rate),
                )
            )
    if not records:
        raise FileNotFoundError("No matching daily SDS MiniSEED files")
    return records


def read_window(records, start, end, *, strict_sample_grid=False):
    """Read a half-open interval [start, end), preserving gaps as missing data.

    Adjacent files are merged without interpolation or zero filling. Conflicting
    overlaps become missing samples. Masked gaps are split into separate traces.
    With strict_sample_grid, reject fractional-sample offsets before ObsPy merge
    can align overlapping records to a common grid. Calibration uses this guard.
    """
    stream = Stream()
    # The margin covers indexed ends a few ms early (see _record_span); extra
    # files are trimmed to [start, end) by read, so the data are unchanged.
    paths = {
        record.path
        for record in records
        if record.start < end and record.end + SELECTION_MARGIN_SECONDS >= start
    }
    for path in sorted(paths):
        stream += read(str(path), format="MSEED", starttime=start, endtime=end)
    if not stream:
        return stream
    if strict_sample_grid:
        grids = {}
        for trace in stream:
            rate = float(trace.stats.sampling_rate)
            reference, reference_rate = grids.setdefault(
                trace.id, (trace.stats.starttime, rate)
            )
            if rate != reference_rate:
                raise ValueError(f"Unequal within-channel sample rates: {trace.id}")
            offset = float(trace.stats.starttime - reference)
            residual = abs(offset - round(offset * rate) / rate)
            if residual > 1e-6:
                raise SampleGridError(trace.id, residual, rate)
    stream.merge(method=0, fill_value=None)
    for trace in stream:
        trace.trim(start, end - trace.stats.delta, nearest_sample=False)
    return stream.split()
