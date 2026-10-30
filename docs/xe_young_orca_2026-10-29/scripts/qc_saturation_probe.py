"""Probe raw XE counts for clipped or dead seismometer channels on sample days."""

import sys
from concurrent.futures import ProcessPoolExecutor
from obspy import UTCDateTime, read
import numpy as np
from config import NETWORK, SDS, SDS_LOCATION, STATIONS

CLIP = 2**23 - 1


def probe(sta):
    out = []
    for y, j in [
        (2018, 130),
        (2018, 160),
        (2018, 200),
        (2018, 260),
        (2018, 320),
        (2019, 20),
        (2019, 80),
        (2019, 120),
    ]:
        row = []
        for ch in ("BHZ", "BH1", "BH2"):
            f = (
                SDS
                / str(y)
                / NETWORK
                / sta
                / f"{ch}.D"
                / f"{NETWORK}.{sta}.{SDS_LOCATION}.{ch}.D.{y}.{j:03d}"
            )
            if not f.exists():
                row.append("  --  ")
                continue
            t = UTCDateTime(year=y, julday=j) + 43200
            tr = read(str(f), starttime=t, endtime=t + 3600)[0]
            d = tr.data.astype(float)
            frac = np.mean(np.abs(d) >= 0.98 * CLIP)
            row.append(
                "CLIP" if frac > 0.01 else ("DEAD" if d.std() < 1 else f"{d.std():.0e}")
            )
        out.append(f"{y}.{j:03d}:" + "/".join(row))
    return sta + " " + " ".join(out)


if __name__ == "__main__":
    with ProcessPoolExecutor(14) as pool:
        for line in pool.map(probe, STATIONS):
            print(line, flush=True)
