"""Quick-look figures for checking processing windows."""

import matplotlib.pyplot as plt
import numpy as np


def plot_window_starts(stream_splitted):
    """Plot start-time distribution for any nonempty list of windows."""
    if not stream_splitted:
        raise ValueError("at least one window is required")
    starts = np.array([float(st[0].stats.starttime) for st in stream_splitted])
    indices = np.unique(np.linspace(0, len(starts) - 1, min(4, len(starts)), dtype=int))
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(starts, bins=max(1, len(starts)))
    ax.set_xticks(
        starts[indices],
        [str(stream_splitted[i][0].stats.starttime) for i in indices],
        rotation=45,
    )
    ax.set(xlabel="Window start", ylabel="Window count")
    fig.tight_layout()
    return fig
