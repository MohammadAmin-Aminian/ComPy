"""V1: GlobSed v3 sediment thickness (Straume et al. 2019) at the archive stations.

An independent map (seismic reflection and refraction compilations), not used by
the inversion. Bilinear interpolation of the 5-arc-minute grid; NaN outside it.

Run: <python with netCDF4> -I v01_globsed.py <GlobSed-v3.nc> <out.csv>
(e.g. ~/Desktop/BRUIT-FM/ComPy-pressure-atlas/pressure_atlas/cache/venv/bin/python)
"""

import csv
import sys
from pathlib import Path

import numpy as np

PLAN = Path.home() / "Desktop/ComPy/12_inversion_all_stations/tables/station_plan.csv"


def read_grid(path):
    """Longitude, latitude and thickness [m] of the grid (netCDF4 needed: the file is HDF5)."""
    import netCDF4

    with netCDF4.Dataset(path) as nc:
        names = list(nc.variables)
        lon = nc.variables[
            [n for n in names if n.lower() in ("lon", "x", "longitude")][0]
        ][:]
        lat = nc.variables[
            [n for n in names if n.lower() in ("lat", "y", "latitude")][0]
        ][:]
        z = nc.variables[
            [n for n in names if n.lower() in ("z", "sedthick", "thickness")][0]
        ][:]
    z = np.ma.filled(np.ma.asarray(z, float), np.nan)
    return np.asarray(lon, float), np.asarray(lat, float), z


def sample(lon, lat, z, x, y):
    x = (x - lon[0]) % 360 + lon[0]
    i = np.clip(np.searchsorted(lon, x) - 1, 0, len(lon) - 2)
    j = np.clip(np.searchsorted(lat, y) - 1, 0, len(lat) - 2)
    tx = (x - lon[i]) / (lon[i + 1] - lon[i])
    ty = (y - lat[j]) / (lat[j + 1] - lat[j])
    return (
        (1 - tx) * (1 - ty) * z[j, i]
        + tx * (1 - ty) * z[j, i + 1]
        + (1 - tx) * ty * z[j + 1, i]
        + tx * ty * z[j + 1, i + 1]
    )


def main(grid, out):
    lon, lat, z = read_grid(grid)
    print(
        "grid",
        z.shape,
        "lon",
        lon[0],
        lon[-1],
        "lat",
        lat[0],
        lat[-1],
        "range",
        np.nanmin(z),
        np.nanmax(z),
    )
    rows = []
    with open(PLAN) as fh:
        for r in csv.DictReader(fh):
            v = sample(lon, lat, z, float(r["longitude"]), float(r["latitude"]))
            rows.append(
                dict(
                    sid=r["sid"],
                    latitude=r["latitude"],
                    longitude=r["longitude"],
                    globsed_m=round(float(v), 1) if np.isfinite(v) else "",
                )
            )
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(
        len(rows),
        "stations;",
        sum(1 for r in rows if r["globsed_m"] != ""),
        "inside the grid",
    )


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
