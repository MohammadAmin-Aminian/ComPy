"""Migration check: what the coverage planner needs to keep reusing products."""

import json
import os

from scripts.check_migration import check_ppsds, source_state


def record(path):
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def test_source_state_needs_same_path_size_and_mtime(tmp_path):
    day = tmp_path / "YV.RR52.00.BDH.D.2013.001"
    day.write_bytes(b"x" * 10)
    item = record(day)
    assert source_state(item) == "unchanged"
    # A copy that does not preserve modification times blocks reuse.
    os.utime(day, ns=(item["mtime_ns"], item["mtime_ns"] + 1_000_000_000))
    assert source_state(item) == "modified"
    day.unlink()
    assert source_state(item) == "missing"


def test_check_ppsds_counts_products_and_sources(tmp_path):
    data = tmp_path / "day"
    data.write_bytes(b"abc")
    product = tmp_path / "RR52" / "ppsd.npz"
    product.parent.mkdir()
    product.write_bytes(b"npz")
    (product.with_name("processing.json")).write_text(
        json.dumps({"source_files": [record(data)]})
    )
    old = tmp_path / "old" / "ppsd.npz"
    old.parent.mkdir()
    old.write_bytes(b"npz")
    registry = tmp_path / "master.csv"
    registry.write_text(
        "Seed_id,ppsd_status,ppsd_npz_file\n"
        f"YV.RR52.00.BDH,complete,{product}\n"
        f"YV.RR52.00.BDH,complete,{product}\n"  # duplicate epoch row counted once
        f"YV.RR01.00.BDH,complete,{old}\n"
        f"YV.RR02.00.BDH,complete,{tmp_path / 'gone.npz'}\n"
        "YV.RR03.00.BDH,not_processed,\n"
    )
    products, sources, examples = check_ppsds(
        [registry, tmp_path / "absent.csv"], tmp_path
    )
    assert products[("YV", "present")] == 1
    assert products[("YV", "no_processing_json")] == 1
    assert products[("YV", "missing_npz")] == 1
    assert products[("-", "missing_registry")] == 1
    assert sources[("YV", "unchanged")] == 1
    assert any("gone.npz" in line for line in examples)
