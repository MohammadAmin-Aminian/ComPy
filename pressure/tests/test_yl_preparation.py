"""Do not skip new download inputs just because a PPSD was previously complete."""

import pytest
from obspy import UTCDateTime
from test_ppsd import dataset as dataset

from compy_pressure.ppsd import process_epoch
from scripts.prepare_yl_test import unchanged_complete_ppsd


@pytest.mark.parametrize("change", ["new_file", "inventory"])
def test_completed_ppsd_reuse_requires_unchanged_inputs(dataset, change):
    row, records, _, output = dataset
    result = process_epoch(
        row,
        records,
        output,
        ppsd_length=512,
        minimum_frequency=0.02,
        chunk_seconds=1024,
    )
    row = row.copy()
    for key, value in result.items():
        row[key] = value
    files = {row.Seed_id: sorted({r.path for r in records})}
    assert unchanged_complete_ppsd(row, files)
    if change == "new_file":
        new = files[row.Seed_id][0].parent / "XX.TEST.00.BDH.D.2020.003"
        new.write_bytes(b"new source requiring header verification during processing")
        files[row.Seed_id].append(new)
    else:
        from pathlib import Path

        inventory = Path(row.inventory_file)
        inventory.write_bytes(inventory.read_bytes() + b"\n")
    assert not unchanged_complete_ppsd(row, files)


@pytest.mark.parametrize("bound", ["start", "end"])
def test_bounded_complete_pilot_is_not_full_archive_ppsd(dataset, bound):
    row, records, _, output = dataset
    origin = min(r.start for r in records)
    result = process_epoch(
        row,
        records,
        output,
        ppsd_length=512,
        minimum_frequency=0.02,
        chunk_seconds=1024,
        **{bound: str(origin + 1024)},
    )
    assert result["ppsd_status"] == "complete"
    assert UTCDateTime(result[f"ppsd_requested_{bound}_time"]) == origin + 1024
    row = row.copy()
    for key, value in result.items():
        row[key] = value
    files = {row.Seed_id: sorted({r.path for r in records})}
    assert not unchanged_complete_ppsd(row, files)
