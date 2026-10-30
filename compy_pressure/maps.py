"""Map axes with continents and coastlines for station and network plots.

Create axes with ``subplot_kw={"projection": PLATE}``, call ``add_land`` and
plot stations with ``transform=PLATE``. Natural Earth outlines are downloaded
by Cartopy on first use and cached.
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature

PLATE = ccrs.PlateCarree()
LAND_COLOR = "#E6E1D6"
OCEAN_COLOR = "#F3F7F9"
COAST_COLOR = "#7D7D7D"
WORLD = (-180, 180, -80, 85)


def add_land(axis, extent=WORLD, *, scale="110m", labels=True):
    """Draw ocean, land, coastlines and a labelled graticule; return the gridliner.

    ``extent`` is (west, east, south, north) in degrees. Use ``scale="50m"``
    or ``"10m"`` for regional zooms, where the 110 m outlines are too coarse.
    """
    axis.set_extent(extent, crs=PLATE)
    axis.add_feature(cfeature.OCEAN.with_scale(scale), facecolor=OCEAN_COLOR, zorder=0)
    axis.add_feature(cfeature.LAND.with_scale(scale), facecolor=LAND_COLOR, zorder=1)
    axis.add_feature(
        cfeature.COASTLINE.with_scale(scale), edgecolor=COAST_COLOR, lw=0.5, zorder=2
    )
    grid = axis.gridlines(draw_labels=labels, lw=0.3, color="#999999", alpha=0.5)
    grid.top_labels = grid.right_labels = False
    return grid
