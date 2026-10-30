"""Identify PSD windows computed from constant (dead or zero-filled) data."""

import numpy as np


def dropout_windows(ppsd):
    """Return a boolean mask of saved PSD windows below the histogram floor.

    ObsPy stores zero power as about -3077 dB. Real pressure or seismic noise
    never reaches the -250 dB floor, so a window with such a bin came from
    constant data: on the archive these windows are below the floor in every
    bin. Saved PPSDs keep them; statistics and plot scaling leave them out.
    """
    values = np.asarray(ppsd.psd_values)
    if not values.size:
        return np.zeros(0, dtype=bool)
    return (values < ppsd.db_bin_edges[0]).any(axis=1)


def screened_values(ppsd):
    """PSD values without dropout windows, or all values if every one is."""
    values = np.asarray(ppsd.psd_values)
    dropouts = dropout_windows(ppsd)
    return values if dropouts.all() else values[~dropouts]
