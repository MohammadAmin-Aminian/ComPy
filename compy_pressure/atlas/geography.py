"""Distances of seafloor gauges to the coastline and to the shelf break.

Both distances are great-circle distances to the nearest point of a densified
Natural Earth line (1:10 m): the coastline (islands included) and the boundary
of the area deeper than 200 m (``bathymetry_K_200``), taken as the shelf break.
The shelf-break distance is signed: negative for a gauge shallower than 200 m
(landward of the break), positive seaward. ``isobath_1000`` gives the same
distance to the 1000 m isobath (``bathymetry_J_1000``), the foot of most slopes.

Natural Earth lines are generalised for a 1:10 million map; their position error
is a few kilometres, which sets the floor of every distance here.
"""

from __future__ import annotations

import numpy as np

EARTH_RADIUS_KM = 6371.0088
NE_LAYERS = {"coast": ("physical", "coastline"), "shelf_break": ("physical", "bathymetry_K_200"),
             "isobath_1000": ("physical", "bathymetry_J_1000")}  # fmt: skip


def _unit(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.column_stack(
        [np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
    )


def _densify(coords, step_km):
    """Points along a lon/lat polyline, at most ``step_km`` apart (great circle)."""
    c = np.asarray(coords, float)[:, :2]
    if len(c) < 2:
        return c
    u = _unit(c[:, 1], c[:, 0])
    out = [c[:1]]
    for i in range(len(c) - 1):
        ang = np.arccos(np.clip(np.dot(u[i], u[i + 1]), -1, 1))
        n = max(1, int(np.ceil(ang * EARTH_RADIUS_KM / step_km)))
        if n == 1:
            out.append(c[i + 1 : i + 2])
            continue
        t = np.linspace(0, 1, n + 1)[1:]
        s = np.sin(ang)
        if s < 1e-12:
            out.append(np.repeat(c[i + 1 : i + 2], n, axis=0))
            continue
        p = (
            np.sin((1 - t) * ang)[:, None] * u[i] + np.sin(t * ang)[:, None] * u[i + 1]
        ) / s
        out.append(
            np.column_stack(
                [
                    np.degrees(np.arctan2(p[:, 1], p[:, 0])),
                    np.degrees(np.arcsin(np.clip(p[:, 2], -1, 1))),
                ]
            )
        )
    return np.vstack(out)


def line_points(layer, step_km=1.0):
    """Densified (lon, lat) points of a Natural Earth 1:10 m layer (``NE_LAYERS`` key)."""
    import cartopy.io.shapereader as shp
    from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon

    category, name = NE_LAYERS[layer]
    path = shp.natural_earth(resolution="10m", category=category, name=name)
    pts = []
    for geom in shp.Reader(path).geometries():
        parts = []
        if isinstance(geom, (Polygon, MultiPolygon)):
            for poly in getattr(geom, "geoms", [geom]):
                parts.append(poly.exterior.coords)
                parts += [r.coords for r in poly.interiors]
        elif isinstance(geom, (LineString, MultiLineString)):
            parts += [g.coords for g in getattr(geom, "geoms", [geom])]
        is_polygon = isinstance(geom, (Polygon, MultiPolygon))
        for p in parts:
            q = _densify(np.asarray(p), step_km)
            if is_polygon:
                # Natural Earth polygons are cut along the antimeridian and closed along the poles; those
                # edges are not isobaths and must not count as the shelf break.
                q = q[(np.abs(q[:, 0]) < 179.999) & (np.abs(q[:, 1]) < 89.999)]
            if len(q):
                pts.append(q)
    return np.vstack(pts)


class NearestLine:
    """Great-circle distance [km] to the nearest point of a set of (lon, lat) points."""

    def __init__(self, lonlat):
        from scipy.spatial import cKDTree

        self.lonlat = np.asarray(lonlat, float)
        self.tree = cKDTree(_unit(self.lonlat[:, 1], self.lonlat[:, 0]))

    def distance_km(self, lat, lon):
        d, i = self.tree.query(_unit(np.atleast_1d(lat), np.atleast_1d(lon)))
        return 2 * EARTH_RADIUS_KM * np.arcsin(np.clip(d / 2, 0, 1)), self.lonlat[i]


def gauge_distances(lat, lon, depth_m, step_km=1.0):
    """Distances of gauges [km]: coast, signed shelf break, signed 1000 m isobath.

    Args:
        lat, lon: gauge positions [deg].
        depth_m: water depth [m], positive down; sets the sign of the isobath distances.
    Returns:
        dict of arrays ``coast_km``, ``shelf_break_km``, ``isobath_1000_km``.
    """
    lat, lon, depth = (np.asarray(x, float) for x in (lat, lon, depth_m))
    out = {}
    for layer, iso in (
        ("coast", None),
        ("shelf_break", 200.0),
        ("isobath_1000", 1000.0),
    ):
        d, _ = NearestLine(line_points(layer, step_km)).distance_km(lat, lon)
        if iso is not None:
            d = np.where(depth < iso, -d, d)
        out[f"{layer}_km"] = d
    return out
