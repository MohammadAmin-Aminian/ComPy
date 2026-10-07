#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 11:13:05 2023

@author: mohammadamin
"""

from compy_numerics import sliding_window, wavenumber, compliance_uncertainty, positive
from compy_streams import (
    split_stream,
    cut_stream_with_overlap,
    trim_streams_to_same_length,
)

__version__ = "2.1.0"

import numpy as np
import matplotlib.pyplot as plt
import tiskitpy as tiskit
from obspy.clients.fdsn import Client
import scipy
import obspy

nhnm = obspy.signal.spectral_estimation.get_nhnm()
nlnm = obspy.signal.spectral_estimation.get_nlnm()


def _water_depth(stream, depth):
    if depth is not None:
        return positive(depth, "water depth")
    if not stream:
        raise ValueError("stream must be nonempty")
    z = stream.select(component="Z")
    if len(z) != 1:
        raise ValueError("exactly one vertical trace is required")
    stats = z[0].stats
    inventory = Client("RESIF").get_stations(
        network=stats.network,
        station=stats.station,
        location=stats.location,
        channel=stats.channel,
        starttime=stats.starttime,
        endtime=stats.endtime,
        level="channel",
    )
    return positive(
        -inventory.get_coordinates(z[0].id, stats.starttime)["elevation"], "water depth"
    )


def _plot_compliance(result, stream):
    curves, coherence, _, fc, f, scatter, _, pressure, vertical = result
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    fig.suptitle(f"{stream[0].stats.network}.{stream[0].stats.station}")
    for ax, values, label in zip(
        axes.flat[:3],
        (vertical, pressure, coherence),
        ("Displacement PSD [m²/Hz]", "Pressure PSD [Pa²/Hz]", "Amplitude coherence"),
    ):
        ax.semilogx(f[1:], np.asarray(values)[:, 1:].T, alpha=0.2)
        ax.semilogx(f[1:], np.median(values, axis=0)[1:], color="black")
        ax.set(xlabel="Frequency [Hz]", ylabel=label)
        ax.grid(True)
    axes.flat[3].errorbar(
        fc, np.mean(curves, axis=0), yerr=scatter, label="Window scatter"
    )
    axes.flat[3].set(xlabel="Frequency [Hz]", ylabel="Normalized compliance [1/Pa]")
    axes.flat[3].grid(True)
    return fig


def Calculate_Compliance(
    stream,
    f_min_com=0.007,
    f_max_com=0.017,
    gain_factor=0.66,
    time_window=1,
    *,
    depth=None,
    nseg=4096,
    plot=True,
):
    """Compliance from displacement (*Z, m) and nominal pressure (*H, Pa).

    gain_factor multiplies nominal pressure. depth is positive water depth in
    metres; providing it avoids an online station metadata request. Return
    the legacy seven-tuple; uncertainty_theory uses amplitude coherence.
    """
    from compy_processing import estimate

    result = estimate(
        stream,
        _water_depth(stream, depth),
        gain_factor,
        time_window,
        f_min_com,
        f_max_com,
        False,
        nseg,
    )
    if plot:
        _plot_compliance(result, stream)
    return result[:7]


def Calculate_Compliance_beta(
    stream,
    f_min_com=0.007,
    f_max_com=0.02,
    gain_factor=0.66,
    time_window=2,
    *,
    depth=None,
    nseg=4096,
    plot=True,
):
    """Compliance with historical beta quality gates; returns six-tuple.

    Input must already have instrument responses removed: vertical
    displacement in m, nominal pressure in Pa. s returned here is standard
    deviation across selected windows, not standard error of the mean.
    """
    from compy_processing import estimate

    result = estimate(
        stream,
        _water_depth(stream, depth),
        gain_factor,
        time_window,
        f_min_com,
        f_max_com,
        True,
        nseg,
    )
    if plot:
        _plot_compliance(result, stream)
    return result[:6]


def start():
    print(" Plot Parametres Loaded")
    plt_params()


# %%


# %%


# %%
def gravitational_attraction(High_Dp, depth_s, f, pw=1025):
    """Surface-wave gravity correction, using pressure amplitude and SI units.

    The pressure/acceleration correction is evaluated without multiplying an
    overflowing cosh by an underflowing exponential. Surface height itself
    may be infinite outside the physically useful infragravity band.
    """
    pw = positive(pw, "water density")
    pressure = np.asarray(High_Dp, dtype=float)
    f = np.asarray(f, dtype=float)
    if np.any(~np.isfinite(pressure)) or np.any(pressure < 0):
        raise ValueError("pressure amplitude must be finite and nonnegative")
    kh = wavenumber(2 * np.pi * f, depth_s) * depth_s
    pa_ratio = 2 * np.pi * (6.6743e-11 / 9.8) * (1 + np.exp(-2 * kh)) / 2
    aw = pa_ratio * pressure
    with np.errstate(over="ignore"):
        hw = pressure * np.cosh(kh) / (pw * 9.8)
    return pa_ratio, aw, hw


# %%
def Comliance_uncertainty(compliance_function, coherence_function, number_of_window):
    """Legacy spelling of amplitude-coherence uncertainty (SI units)."""
    return compliance_uncertainty(
        compliance_function, coherence_function, number_of_window
    )


# %%
def Rotate(stream, time_window=1, *, plot=True):
    """Rotate and remove horizontal coherent noise in complete hour windows.

    Return (stream, azimuth_deg, angle_deg, variance_after_over_before).
    Processing failures raise with window context; input and seed metadata
    are preserved. A final incomplete window is omitted.
    """
    windows = split_stream(stream, positive(time_window, "time_window") * 3600)
    if not windows:
        raise ValueError("no complete rotation windows")
    azimuth, angles, variance_ratio, cleaned = [], [], [], []
    for i, window in enumerate(windows):
        try:
            before = np.var(window.select(component="Z")[0].data)
            if before <= 0:
                raise ValueError("vertical trace has zero variance")
            rotator = tiskit.CleanRotator(window, remove_eqs=False)
            rotated = rotator.apply(window, set_dtype=False)
            cleaner = tiskit.DataCleaner(
                rotated, ["*1", "*2"], n_to_reject=0, remove_eqs=False
            )
            rotated = cleaner.apply(rotated, set_dtype=False)
            after = np.var(rotated.select(component="Z")[0].data)
        except Exception as exc:
            raise RuntimeError(
                f"rotation failed in window {i} starting {window[0].stats.starttime}"
            ) from exc
        azimuth.append(rotator.azimuth % 360)
        angles.append(rotator.angle)
        variance_ratio.append(after / before)
        cleaned.extend(rotated)
    output = obspy.Stream(cleaned)
    output.merge(method=0)  # no interpolation of gaps or conflicting samples
    if plot:
        fig, axes = plt.subplots(2, 1, figsize=(12, 7), constrained_layout=True)
        for ax, values, label in zip(
            axes, (azimuth, angles), ("Azimuth [deg]", "Tilt [deg]")
        ):
            points = ax.scatter(
                np.arange(len(windows)) * time_window, values, c=variance_ratio
            )
            ax.set(xlabel="Hours since first complete window", ylabel=label)
            fig.colorbar(points, ax=ax, label="Variance after/before")
    return output, np.asarray(azimuth), np.asarray(angles), np.asarray(variance_ratio)


def Rotate_angles(stream, time_window=1):
    """Estimate angles and fractional variance reduction without applying cleaning."""
    windows = split_stream(stream, positive(time_window, "time_window") * 3600)
    if not windows:
        raise ValueError("no complete rotation windows")
    result = []
    for i, window in enumerate(windows):
        try:
            window.detrend("simple")
            rotator = tiskit.CleanRotator(window, remove_eqs=False)
        except Exception as exc:
            raise RuntimeError(f"angle estimation failed in window {i}") from exc
        result.append(
            (rotator.azimuth % 360, rotator.angle, rotator.variance_reduction)
        )
    return tuple(np.asarray(result).T)

    # return(azimuth,angle)


# %%

# return(azimuth,angle)


# %%

# %%
from obspy import read, Stream


# %%


# %%
def rms(arr):
    """Root-mean-square amplitude with floating-point and complex support."""
    from ffplot import rms as amplitude

    return amplitude(arr)


# %%
def split_and_save_stream(stream, interval_minutes, output_dir):
    """Save consecutive copied windows, including the final partial window."""
    from pathlib import Path
    from compy_streams import validate_stream, sample_count

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


# %%
def optimizer(
    Com,
    Czp,
    stream,
    f,
    alpha=0.95,
    beta=0.70,
    zeta=0.4,
    f_min_com=0.007,
    f_max_com=0.018,
):
    coherence_mask = (f >= f_min_com) & (f <= f_max_com)
    coherence_mask2 = (f >= 0.025) & (f <= 0.03)
    High_Czp = []
    High_Com = []
    High_Com_Stream = []
    for i in range(0, len(Czp)):
        if np.mean(Czp[i][coherence_mask]) > alpha:
            if np.mean(Czp[i][coherence_mask2]) < beta:
                if np.min(Czp[i][coherence_mask]) > zeta:
                    # Com1 = (k * Czp[i]* (np.sqrt(Dz[i]))) / (np.sqrt(Dp[i]) / gain_factor)

                    High_Czp.append(Czp[i])

                    High_Com.append(Com[i])

                    High_Com_Stream.append(stream[i])

                    print(i)

    # plt.rcParams.update({'font.size': 25})
    # plt.figure(dpi=300,figsize=(12,8))
    # for i in range(0,len(Czp)):
    #     plt.semilogx(f,Czp[i],linewidth = 0.5,color='r')
    # for i in range(0,len(High_Czp)):
    #     plt.semilogx(f,High_Czp[i],linewidth = 0.5,color='g')
    # plt.semilogx(f,np.median(High_Czp,axis=0),linewidth = 2,color='b',label='Median of optimized')
    # plt.xlabel('Frequency [Hz]')
    # plt.ylabel('Coherence')
    # plt.grid(True)
    # plt.xlim([0.001,1])
    # plt.vlines(x = f_min_com, ymin=0, ymax=1,color='black',linestyles="dashed",label="High Coherence Band")
    # plt.vlines(x = f_max_com, ymin=0, ymax=1,color='black',linestyles="dashed")

    # plt.legend(loc='upper right',fontsize=17)

    return (High_Com, High_Czp, High_Com_Stream)


# %%
def optimizer_rms(Com, Czp, stream, f, a1, a2, percentage=10, alpha=0.9, beta=0.5):

    High_Czp = []
    High_Com = []
    High_Com_Stream = []

    rms_data = rms(np.median(Com, axis=0)[a1:a2])
    for i in range(0, len(Com)):
        if (
            rms_data * (1 - (percentage) / 100)
            < rms(Com[i][a1:a2])
            < rms_data * (1 + (percentage) / 100)
        ):
            if 1 > np.median(Czp[i][a1:a2]) > alpha:
                if 1 > np.mean(Czp[i][a1:a2]) > beta:
                    High_Czp.append(Czp[i])

                    High_Com.append(Com[i])

                    High_Com_Stream.append(stream[i])

                    print(i)

    # plt.rcParams.update({'font.size': 25})
    # plt.figure(dpi=300,figsize=(12,8))
    # for i in range(0,len(Czp)):
    #     plt.semilogx(f,Czp[i],linewidth = 0.5,color='r')
    # for i in range(0,len(High_Czp)):
    #     plt.semilogx(f,High_Czp[i],linewidth = 0.5,color='g')
    # plt.semilogx(f,np.median(High_Czp,axis=0),linewidth = 2,color='b',label='Median of optimized')
    # plt.xlabel('Frequency [Hz]')
    # plt.ylabel('Coherence')
    # plt.grid(True)
    # plt.xlim([0.001,1])
    # plt.vlines(x = f_min_com, ymin=0, ymax=1,color='black',linestyles="dashed",label="High Coherence Band")
    # plt.vlines(x = f_max_com, ymin=0, ymax=1,color='black',linestyles="dashed")

    # plt.legend(loc='upper right',fontsize=17)

    return (High_Com, High_Czp, High_Com_Stream)


# %%


# %%
def overlap_checker(stream_splitted):
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


# %%
def plt_params():
    plt.rcParams["font.size"] = 40
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["font.family"] = ["STIXGeneral"]
    plt.rcParams["font.weight"] = "normal"
    plt.rcParams["mathtext.default"] = "regular"

    plt.rcParams["axes.grid"] = True
    plt.rcParams["axes.grid.which"] = "major"
    plt.rcParams["axes.labelcolor"] = "black"
    plt.rcParams["axes.labelpad"] = 4.0
    plt.rcParams["axes.labelsize"] = "medium"
    plt.rcParams["axes.labelweight"] = "normal"
    plt.rcParams["axes.linewidth"] = 0.6
    plt.rcParams["axes.titlecolor"] = "black"
    plt.rcParams["axes.titlepad"] = 5.0
    plt.rcParams["axes.titlesize"] = "large"
    plt.rcParams["axes.titleweight"] = "normal"
    plt.rcParams["axes.xmargin"] = 0.05
    plt.rcParams["axes.ymargin"] = 0.05

    plt.rcParams["grid.alpha"] = 1
    plt.rcParams["grid.color"] = "#b0b0b0"
    plt.rcParams["grid.linestyle"] = "--"
    plt.rcParams["grid.linewidth"] = 0.6

    plt.rcParams["xtick.color"] = "black"
    plt.rcParams["xtick.direction"] = "out"
    plt.rcParams["xtick.labelsize"] = 35
    plt.rcParams["xtick.major.pad"] = 2.0
    plt.rcParams["xtick.minor.pad"] = 2.0
    plt.rcParams["xtick.minor.visible"] = True
    plt.rcParams["ytick.color"] = "black"
    plt.rcParams["ytick.direction"] = "out"
    plt.rcParams["ytick.labelsize"] = 35
    plt.rcParams["ytick.major.pad"] = 2.0
    plt.rcParams["ytick.minor.pad"] = 2.0
    plt.rcParams["ytick.minor.visible"] = True
    plt.rcParams["xtick.minor.size"] = 10  # Length of minor ticks on x-axis
    plt.rcParams["ytick.minor.size"] = 10  # Length of minor ticks on y-axis
    plt.rcParams["ytick.minor.width"] = 1  # Width of minor ticks on y-axis

    # Setting parameters for major ticks
    plt.rcParams["xtick.major.size"] = 10  # Length of major ticks on x-axis
    plt.rcParams["ytick.major.size"] = 10  # Length of major ticks on y-axis
    plt.rcParams["xtick.major.width"] = 2  # Width of major ticks on x-axis
    plt.rcParams["ytick.major.width"] = 2  # Width of major ticks on y-axis
    plt.rcParams["xtick.minor.width"] = 1  # Width of minor ticks on x-axis

    plt.rcParams["legend.borderaxespad"] = 0.0
    plt.rcParams["legend.borderpad"] = 0.5
    plt.rcParams["legend.columnspacing"] = 1.5
    plt.rcParams["legend.edgecolor"] = "gray"
    plt.rcParams["legend.facecolor"] = "white"
    plt.rcParams["legend.fancybox"] = False
    plt.rcParams["legend.fontsize"] = "small"
    plt.rcParams["legend.framealpha"] = 0.8
    plt.rcParams["legend.handleheight"] = 1.0
    plt.rcParams["legend.handlelength"] = 2.0
    plt.rcParams["legend.handletextpad"] = 0.5
