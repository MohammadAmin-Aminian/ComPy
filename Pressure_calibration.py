#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 16:33:13 2023

@author: Mohammad-Amin Aminian

DFG Calibration

"""

from compy_numerics import plot_smooth, sliding_window, positive, residual
from obspy.signal.trigger import plot_trigger

# from obspy.signal.trigger import coincidence_trigger
# from obspy.signal.trigger import classic_sta_lta
import matplotlib.pyplot as plt
import scipy
from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
import numpy as np
import obspy
import tiskitpy as tiskit
from disba import PhaseDispersion


def calculate_spectral_ratio(
    stream,
    mag=7,
    coh_trsh=0.97,
    mean_trsh=0.97,
    f_min=0.02,
    f_max=0.06,
    plot_condition=False,
    *,
    inventory=None,
    event_spans=None,
    client=None,
    nseg=512,
    velocity_model=None,
):
    """Estimate pressure gain from response-corrected teleseismic P/a ratios.

    Input is raw counts. Supply inventory and event_spans to operate offline.
    Coherence thresholds here are magnitude-squared coherence, as in v1.
    Each event is processed independently; incomplete/low-coherence events
    are excluded and an empty accepted set raises ValueError.
    """
    from compy_numerics import integer
    from compy_processing import channel_pair, spectra
    from compy_streams import trim_streams_to_same_length

    channel_pair(stream)
    nseg = integer(nseg, "nseg", minimum=2)
    if not 0 < f_min < f_max < stream[0].stats.sampling_rate / 2:
        raise ValueError("calibration band must be positive and below Nyquist")
    if not 0 <= coh_trsh <= 1 or not 0 <= mean_trsh <= 1:
        raise ValueError("coherence thresholds must be in [0, 1]")
    z, p = channel_pair(stream)
    if inventory is None:
        client = Client("RESIF") if client is None else client
        inventory = client.get_stations(
            network=z.stats.network,
            station=z.stats.station,
            channel="*",
            location="*",
            starttime=z.stats.starttime,
            endtime=z.stats.endtime,
            level="response",
        )
    depth = positive(
        -inventory.get_coordinates(z.id, z.stats.starttime)["elevation"], "water depth"
    )
    if event_spans is None:
        event_spans = tiskit.TimeSpans.from_eqs(
            (z.stats.starttime, z.stats.endtime), minmag=mag, days_per_magnitude=0.5
        )
    start_times = list(event_spans.start_times)
    if not start_times:
        raise ValueError("no earthquake windows were found")
    ratios, coherences = [], []
    for start in start_times:
        event = (
            (stream.select(id=z.id) + stream.select(id=p.id))
            .slice(start, start + 7200)
            .copy()
        )
        if len(event) != 2 or min(len(tr) for tr in event) < nseg:
            continue
        event = trim_streams_to_same_length(event, (z.stats.channel, p.stats.channel))
        event.select(id=z.id).remove_response(inventory=inventory, output="ACC")
        event.select(id=p.id).remove_response(inventory=inventory, output="DEF")
        event.filter(
            "bandpass",
            freqmin=0.005,
            freqmax=min(0.1, event[0].stats.sampling_rate * 0.45),
        )
        event, _, _ = rayleigh_arrival(
            event, timelag=-5, window=20, plot_condition=False
        )
        if min(len(tr) for tr in event) < nseg:
            continue
        f, aa, pp, _, amplitude_coherence = spectra(event, nseg)
        band = (f >= f_min) & (f <= f_max)
        if not band.any():
            raise ValueError("calibration band contains no spectral bins")
        coherence = amplitude_coherence**2
        if (
            np.median(coherence[band]) <= coh_trsh
            or np.mean(coherence[band]) <= mean_trsh
        ):
            continue
        if np.any(aa[band] <= 0):
            continue
        ratio = np.zeros_like(f)
        np.divide(pp, aa, out=ratio, where=aa > 0)
        ratios.append(np.sqrt(ratio) / (1028 * depth))
        coherences.append(coherence)
    if not ratios:
        raise ValueError("no earthquake windows passed calibration quality thresholds")
    measured = np.median(ratios, axis=0)
    sound_speed = calculate_speed_of_sound_in_water(depth=depth)
    periods = np.sort(1 / f[band])
    theoretical_f, theoretical = theoretical_p_a_ratio(
        alpha=sound_speed,
        h=depth,
        t=periods,
        plot_condition=False,
        velocity_model=velocity_model,
    )
    order = np.argsort(theoretical_f)
    if len(theoretical_f) != len(periods) or not np.allclose(
        theoretical_f[order], f[band]
    ):
        raise ValueError(
            "dispersion solver did not return the full calibration frequency band"
        )
    target = theoretical[order]
    if not np.all(np.isfinite(target)) or np.any(target <= 0):
        raise ValueError(
            "calibration band crosses a resonance or invalid theoretical ratio"
        )
    gain = grid_search(measured[band], target)
    if plot_condition:
        fig, axes = plt.subplots(2, 1, figsize=(12, 9), constrained_layout=True)
        axes[0].semilogx(f[1:], np.asarray(coherences)[:, 1:].T, alpha=0.3)
        axes[0].set(ylabel="Magnitude-squared coherence", ylim=(0, 1))
        axes[1].semilogx(f[band], measured[band], label="Measured P/a / ρh")
        axes[1].semilogx(f[band], target, label="Theoretical")
        axes[1].semilogx(f[band], measured[band] * gain, label="Calibrated")
        axes[1].set(xlabel="Frequency [Hz]", ylabel="P/a / ρh")
        axes[1].legend()
    return gain


# %%
def rayleigh_arrival(stream, window=20, timelag=-2, plot_condition=False):

    max_index = np.argmax(
        np.abs(stream.select(channel="*Z")[0].data)
    )  # Find the maximum value in the trace

    fs = stream.select(channel="*Z")[0].stats.sampling_rate

    max_time = max_index / fs  # Calculate the time of the maximum value

    ray_arr = stream.select(channel="*Z")[0].stats.starttime + max_time + (timelag * 60)
    end_time = (
        stream.select(channel="*Z")[0].stats.starttime
        + max_time
        + ((timelag + window) * 60)
    )

    print(f"Fundamental Rayleigh arrival at {ray_arr} seconds.")

    st = stream.copy()

    st.trim(starttime=ray_arr, endtime=end_time)

    if plot_condition:
        plt.figure(dpi=300, figsize=(12, 6))
        plt.subplot(211)
        plt.plot(stream.select(channel="*Z")[0].data)
        plt.vlines(
            max_index + (timelag * 60) * fs,
            np.min(stream.select(channel="*Z")[0].data),
            np.max(stream.select(channel="*Z")[0].data),
            color="r",
        )
        plt.vlines(
            max_index + ((timelag + window) * 60 * fs),
            np.min(stream.select(channel="*Z")[0].data),
            np.max(stream.select(channel="*Z")[0].data),
            color="r",
        )

        plt.subplot(212)
        plt.plot(st.select(channel="*Z")[0].data)

    return (st, (max_time + (timelag * 60)), (max_time + ((timelag + window) * 60)))


# %%
def misfit(d, m, l=2, s=1):
    """Sum of absolute normalized residuals to power l (legacy calibration loss)."""
    l = positive(l, "norm power")
    return np.sum(np.abs(residual(d, m, s)) ** l)


# %%
# it can be better by writing the code for step size, you did before somewhere!!!
def grid_search(d, m):
    """Fit the positive gain multiplying measured d to match theoretical m.

    Analytic least squares replaces a 0.01-quantized search limited to [0,5),
    and avoids division by zero when the old search chose its zero endpoint.
    """
    d, m = np.asarray(d, float), np.asarray(m, float)
    residual(d, m)
    denominator = np.sum(d * d)
    if denominator == 0:
        raise ValueError("cannot calibrate an all-zero measured spectrum")
    return positive(np.sum(d * m) / denominator, "fitted gain")


# %%
def theoretical_p_a_ratio(
    alpha=1500, rho=1028, h=4760, plot_condition=True, t=None, velocity_model=None
):
    """Water-column P/a divided by rho*h, using angular frequency.

    Model units follow disba: km, km/s, km/s, g/cm³. The removable zero
    vertical-wavenumber singularity uses its analytic limit of one.
    """
    alpha, rho, h = (
        positive(alpha, "sound speed"),
        positive(rho, "water density"),
        positive(h, "water depth"),
    )
    f, phase_velocity = phase_dispersion(
        velocity_model=velocity_model, plot_condition=False, t=t
    )
    q = (
        2
        * np.pi
        * f
        * h
        * np.sqrt((alpha**-2 - (phase_velocity * 1000) ** -2).astype(complex))
    )
    ratio = np.ones_like(q)
    np.divide(np.tan(q), q, out=ratio, where=np.abs(q) > 1e-12)
    ratio = np.real_if_close(ratio).real
    if plot_condition:
        plt.figure(figsize=(10, 5))
        plt.semilogx(f, ratio)
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Theoretical P/a / ρh")
        plt.grid(True)
    return f, ratio


# %%
def phase_dispersion(velocity_model=None, plot_condition=True, t=None):
    """Fundamental Rayleigh dispersion: model [km, km/s, km/s, g/cm³]."""
    if velocity_model is None:
        velocity_model = np.array(
            [
                [0.26, 1.75, 0.34, 1.84],
                [0.33, 5.00, 2.70, 2.55],
                [0.33, 5.00, 2.70, 2.55],
                [0.60, 6.50, 3.70, 2.85],
                [0.60, 6.50, 3.70, 2.85],
                [0.60, 6.50, 3.70, 2.85],
                [1, 7.10, 4.05, 3.05],
                [1, 7.10, 4.05, 3.05],
                [1, 7.10, 4.05, 3.05],
                [1, 7.10, 4.05, 3.05],
                [1, 7.10, 4.05, 3.05],
                [2, 7.60, 4.35, 3.25],
                [2, 7.60, 4.35, 3.25],
                [2, 7.60, 4.35, 3.25],
                [2, 7.60, 4.35, 3.25],
                [2, 7.60, 4.35, 3.25],
                [10, 8.02, 4.39, 3.38],
                [10, 8.02, 4.39, 3.38],
                [10, 8.02, 4.39, 3.38],
                [10, 8.02, 4.39, 3.38],
                [10, 8.02, 4.39, 3.38],
                [20, 8.1, 4.4, 3.4],
                [20, 8.1, 4.4, 3.4],
                [20, 8.1, 4.4, 3.4],
                [20, 8.1, 4.4, 3.4],
                [20, 8.1, 4.4, 3.4],
                [25, 8.05, 4.43, 3.36],
                [25, 8.05, 4.43, 3.36],
                [50, 8.5, 4.6, 3.43],
                [50, 8.7, 4.7, 3.49],
                [50, 8.8, 4.75, 3.53],
                [100, 9.38, 5.07, 3.78],
                [100, 9.9, 5.37, 3.92],
                [200, 10.1, 5.53, 3.98],
                [400, 11.6, 6.45, 4.68],
            ]
        )

    velocity_model = np.asarray(velocity_model, float)
    if (
        velocity_model.ndim != 2
        or velocity_model.shape[1] != 4
        or len(velocity_model) < 2
    ):
        raise ValueError(
            "dispersion model must have at least two rows and four columns"
        )
    if (
        not np.all(np.isfinite(velocity_model))
        or np.any(velocity_model[:, 0] < 0)
        or np.any(velocity_model[:, 1:] <= 0)
    ):
        raise ValueError("invalid dispersion model values")
    if np.any(velocity_model[:, 1] <= velocity_model[:, 2]):
        raise ValueError("Vp must exceed Vs")
    t = np.logspace(0, 3, 1000) if t is None else np.asarray(t, float)
    if (
        t.ndim != 1
        or not len(t)
        or np.any(~np.isfinite(t))
        or np.any(t <= 0)
        or np.any(np.diff(t) <= 0)
    ):
        raise ValueError("periods must be finite, positive, and strictly increasing")
    curve = PhaseDispersion(*velocity_model.T)(t, mode=0, wave="rayleigh")
    if not len(curve.period):
        raise ValueError("dispersion solver returned no periods")
    if plot_condition:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
        depth = np.r_[0, np.cumsum(velocity_model[:, 0])]
        axes[0].step(
            np.r_[velocity_model[:, 2], velocity_model[-1, 2]],
            depth,
            where="post",
            label="Vs",
        )
        axes[0].step(
            np.r_[velocity_model[:, 1], velocity_model[-1, 1]],
            depth,
            where="post",
            label="Vp",
        )
        axes[0].invert_yaxis()
        axes[0].set(xlabel="Velocity [km/s]", ylabel="Depth [km]")
        axes[0].legend()
        axes[1].semilogx(1 / curve.period, curve.velocity)
        axes[1].set(xlabel="Frequency [Hz]", ylabel="Rayleigh phase velocity [km/s]")
    return 1 / curve.period, curve.velocity


# %%
def cut_signal_above_zero(signal):
    """
    Cuts a signal from where it has a value above zero and removes trailing zeros.

    Parameters:
    - signal (list): A list of values representing the signal.

    Returns:
    - cut_signal (list): The cut signal without trailing zeros.
    """
    start_index = 0

    for i, value in enumerate(signal):
        if value > 0:
            start_index = i
            break

    cut_signal = signal[start_index:]

    # Remove trailing zeros
    end_index = len(cut_signal)
    for i, value in enumerate(cut_signal[::-1]):
        if value != 0:
            end_index = len(cut_signal) - i
            break

    cut_signal = cut_signal[:end_index]
    return cut_signal


# %%
def calculate_speed_of_sound_in_water(temperature=4, salinity=35, depth=4760):
    """Calculate the speed of sound in water using Mackenzie's formula."""
    # Constants
    a1 = 1448.96
    a2 = 4.591
    a3 = -5.304 * 10 ** (-2)
    a4 = 2.374 * 10 ** (-4)
    a5 = 1.340
    a6 = 1.630 * 10 ** (-2)
    a7 = 1.675 * 10 ** (-7)
    a8 = -1.025 * 10 ** (-2)
    a9 = -7.139 * 10 ** (-13)

    t = temperature
    s = salinity
    d = depth

    # Calculate speed of sound
    speed_of_sound = (
        a1
        + a2 * t
        + a3 * t**2
        + a4 * t**3
        + a5 * (s - 35)
        + a6 * d
        + a7 * d**2
        + a8 * t * (s - 35)
        + a9 * t * d**3
    )
    return speed_of_sound


# %%
from math import sin, cos, sqrt, atan2, radians


def calculate_distance(lat1, lon1, lat2, lon2):
    # Convert coordinates to radians
    lat1_rad = radians(lat1)
    lon1_rad = radians(lon1)
    lat2_rad = radians(lat2)
    lon2_rad = radians(lon2)

    # Earth's radius in kilometers
    radius = 6371.0

    # Calculate the differences between the coordinates
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad

    # Haversine formula
    a = sin(dlat / 2) ** 2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    distance = radius * c

    return distance


# %%
# Calibrating via P-wave Arrival
# start =UTCDateTime("2012-10-12")
net = "YV"
sta = "RR38"


def pressure_calibration(stream, mag=7, i=1):
    """
    Pressure gauge calibration using p wave arrival

    """
    client = Client("RESIF")
    stream = stream.copy()
    invz = client.get_stations(
        network=stream[0].stats.network,
        station=stream[0].stats.station,
        channel="BHZ",
        location="*",
        level="response",
    )

    invp = client.get_stations(
        network=stream[0].stats.network,
        station=stream[0].stats.station,
        channel="BDH",
        location="*",
        level="response",
    )
    eq_spans = tiskit.TimeSpans.from_eqs(
        (
            stream.select(channel="*Z")[0].stats.starttime,
            stream.select(channel="*Z")[0].stats.endtime,
        ),
        minmag=mag,
        days_per_magnitude=0.5,
    )

    stream.trim(eq_spans.start_times[i], eq_spans.start_times[i] + 10 * 3600)

    stream.select(channel="*Z").remove_response(
        inventory=invz, output="ACC", plot=False
    )

    stream.select(channel="*H").remove_response(
        inventory=invp, output="DEF", plot=False
    )

    df = stream[0].stats.sampling_rate

    cft = obspy.signal.trigger.recursive_sta_lta(
        stream.select(channel="*Z")[0].data, int(5 * df), int(10 * df)
    )

    # trig = coincidence_trigger("recstalta", 7, 2, stream, 1, sta=int(5 * df), lta=int(10 * df))

    p_arrival_times = obspy.signal.trigger.trigger_onset(cft, 1, 0.2) / df
    plot_trigger(stream.select(channel="*Z")[0], cft, 1, 0.5)

    stream.trim(
        stream[0].stats.starttime + p_arrival_times[i - 1][0] - 2,
        stream[0].stats.starttime + p_arrival_times[i - 1][0] + 2,
    )
    stream.detrend()
    rho = 1028
    Acc = stream[0].data / (rho * -invz[0][0][0].elevation)

    plt.plot(stream[3].data), plt.plot(Acc)

    plt.plot(stream[3].data / -Acc), plt.plot(stream[3].data)

    scipy.stats.pearsonr(stream[0].data, stream[3].data)


# %%
def p_calibration(stream, gain_factor, rho=1025, mag=6):

    nseg = 2**11
    TP = 5

    server = "RESIF"
    client = Client(server)

    net = stream[0].stats.network
    sta = stream[0].stats.station

    invz = client.get_stations(
        network=net, station=sta, channel="BHZ", location="*", level="response"
    )

    invp = client.get_stations(
        network=net, station=sta, channel="BDH", location="*", level="response"
    )

    stream.filter("bandpass", freqmin=0.01, freqmax=0.05, corners=4, zerophase=True)

    st1 = stream.copy()
    st2 = stream.copy()
    st3 = stream.copy()

    stream22 = read()
    stream22.clear()

    stream2 = read()
    stream2.clear()

    st1.select(channel="*Z").remove_response(inventory=invz, output="ACC", plot=False)
    # st2.select(channel="*Z").remove_response(inventory=invz,
    #                                                           output="VEL", plot=False)
    # st3.select(channel="*Z").remove_response(inventory=invz,
    #                                                           output="DISP", plot=False)

    st1.select(channel="*H").remove_response(inventory=invp, output="DEF", plot=False)

    stream.sort()

    eq_spans = tiskit.TimeSpans.from_eqs(
        (
            stream.select(channel="*Z")[0].stats.starttime,
            stream.select(channel="*Z")[0].stats.endtime,
        ),
        minmag=mag,
        days_per_magnitude=0.5,
    )

    print(
        str(len(eq_spans))
        + " Earthquakes with magnitude greater than "
        + str(mag)
        + " Mw has been found"
    )

    for i in range(0, len(eq_spans)):
        stream1 = stream.copy()
        stream1 = st1.copy()
        stream22 = stream22 + stream1.trim(
            eq_spans.start_times[i], eq_spans.start_times[i] + 2 * 3600
        )

    Czp = np.zeros([len(eq_spans), int(nseg / 2 + 1)])

    # Calculation coherence

    print("Calculating Coherence and Spectral ratio")
    for i in range(0, len(eq_spans)):
        f, Czp[i] = scipy.signal.coherence(
            stream22[i * 4 : (i + 1) * 4].select(component="Z")[0].data,
            stream22[i * 4 : (i + 1) * 4].select(component="H")[0].data,
            fs=stream22[i * 4 : (i + 1) * 4][0].stats.sampling_rate,
            nperseg=nseg,
            noverlap=(nseg * 0.5),
            window=scipy.signal.windows.tukey(
                nseg,
                (TP * 60 * stream22[i * 4 : (i + 1) * 4][0].stats.sampling_rate) / nseg,
            ),
        )

    coherence_mask = (f >= 0.01) & (f <= 0.1)

    print(
        "Selecting Earthquakes with median of Coherence greater than "
        + str(0.95)
        + " and mean greater than "
        + str(0.8)
    )

    High_Czp = []
    coh_trsh = 0.95
    mean_trsh = 0.80
    for i in range(0, len(eq_spans)):
        if (
            np.median(Czp[i][coherence_mask]) > coh_trsh
            and np.mean(Czp[i][coherence_mask]) > mean_trsh
        ):
            High_Czp.append(Czp[i])
            # print(i)

    plt.rcParams.update({"font.size": 25})
    plt.figure(dpi=300, figsize=(16, 12))
    for i in range(0, len(Czp)):
        plt.semilogx(f, Czp[i], linewidth=1.5, color="r")

    for i in range(0, len(High_Czp)):
        plt.semilogx(f, High_Czp[i], linewidth=3, color="g")

    plt.title("Coherence Of Teleseismic events")
    plt.vlines(
        x=0.01,
        ymin=0,
        ymax=1,
        color="black",
        linestyles="dashed",
        label="Frequency limits ",
        linewidth=3,
    )
    plt.vlines(x=0.1, ymin=0, ymax=1, color="black", linestyles="dashed", linewidth=3)
    plt.ylabel("Coherence ")
    plt.xlabel("Frequency [Hz] ")
    plt.xlim([0.005, 0.2])
    plt.grid(True)
    plt.legend(loc="upper left")
    plt.tight_layout()

    plt.figure(dpi=300, figsize=(20, 12))
    plt.subplot(211)
    plt.plot(
        stream22[i * 4 : (i + 1) * 4].select(component="Z")[0].data[1500:2000],
        label="Acceleration Vertical",
        linewidth=3,
    )
    plt.plot(
        stream22[i * 4 : (i + 1) * 4].select(component="H")[0].data[1500:2000]
        / (rho * invz[0][0][0].elevation),
        linestyle="dotted",
        label="Pressure /ρh",
        linewidth=3,
    )
    plt.title("station " + sta)

    plt.legend(loc="upper right")
    plt.subplot(212)
    plt.plot(
        stream22[i * 4 : (i + 1) * 4].select(component="Z")[0].data[1500:2000],
        label="Acceleration Vertical",
        linewidth=3,
    )
    plt.plot(
        stream22[i * 4 : (i + 1) * 4].select(component="H")[0].data[1500:2000]
        * gain_factor
        / (rho * invz[0][0][0].elevation),
        linestyle="dotted",
        label="Pressure * gain factor /ρh",
        linewidth=3,
    )
    plt.title("Calibrated, Gain Factor = " + str(gain_factor))

    plt.legend(loc="upper right")
    plt.tight_layout()

    plt.figure(dpi=300, figsize=(12, 6))
    plt.plot(
        (st1.select(channel="*Z")[0].data / np.max(st1.select(channel="*Z")[0].data)),
        label="Acceleration Vertical",
    )
    # plt.plot((st1.select(channel="*H1")[0].data / np.max(st1.select(channel="*H1")[0].data)),label="Acceleration H1")
    # plt.plot((st1.select(channel="*H2")[0].data / np.max(st1.select(channel="*H2")[0].data)),label="Acceleration H2")
    plt.plot(
        (-st1.select(channel="*H")[0].data / np.max(st1.select(channel="*H")[0].data)),
        label="Pressure",
    )

    plt.legend(loc="upper left")

    plt.figure(dpi=300, figsize=(12, 6))
    plt.plot((st1.select(channel="*Z")[0].data), label="Acceleration Vertical")
    # plt.plot((st1.select(channel="*H1")[0].data / np.max(st1.select(channel="*H1")[0].data)),label="Acceleration H1")
    # plt.plot((st1.select(channel="*H2")[0].data / np.max(st1.select(channel="*H2")[0].data)),label="Acceleration H2")
    plt.plot(
        gain_factor
        * (st1.select(channel="*H")[0].data / (rho * invz[0][0][0].elevation)),
        label="Pressure",
    )

    plt.legend(loc="upper left")


# %%
import matplotlib as mpl
import scipy.signal


def _date_ticks(start, end, bins, *, spacing=7 * 86400):
    """Date ticks within the plotted bin range, including short recordings."""
    if bins < 1:
        raise ValueError("spectrogram must contain at least one time bin")
    duration = float(end - start)
    intervals = max(1, int(np.ceil(duration / spacing)))
    count = min(bins, intervals + 1)
    positions = np.linspace(0, bins - 1, count)
    times = np.linspace(0, max(0, duration), count)
    return positions, [str(start + offset)[:10] for offset in times]


def plot_spectrogram(raw_stream):
    """
    Plots the spectrogram of a raw stream and highlights specific frequency bands.

    Parameters:
    - raw_stream: The raw data stream containing seismic or other time series data.
    """
    # Define the number of segments for the spectrogram
    from compy_streams import validate_stream
    validate_stream(raw_stream)
    z, p = raw_stream.select(component="Z"), raw_stream.select(component="H")
    if len(z) != 1 or len(p) != 1:
        raise ValueError("exactly one vertical and pressure channel are required")
    nseg = min(2**14, len(z[0]), len(p[0]))
    if nseg < 2:
        raise ValueError("spectrogram needs at least two samples")

    # Compute the spectrogram
    f, t, Sp = scipy.signal.spectrogram(
        raw_stream.select(component="H")[0].data,
        fs=raw_stream[0].stats.sampling_rate,
        nperseg=nseg,
        noverlap=(nseg / 2),
        window="hann",
        scaling="density",
    )

    # Calculate the time difference between consecutive spectrogram points in hours
    # time_diff_hours = (t[1] - t[0]) / 3600

    tick_positions, dates = _date_ticks(
        raw_stream[0].stats.starttime, raw_stream[0].stats.endtime, len(t)
    )

    # Create a new time array for plotting
    t2 = np.arange(0, len(t))

    # Normalize pressure values for color mapping in the plot
    norm_p = mpl.colors.Normalize(vmin=0, vmax=60)

    # Define frequency bands for infra-gravity and microseismic analyses
    f_min_ig = 0.005
    f_max_ig = 0.02
    f1_ig = np.argmin(np.abs(f - f_min_ig))
    f2_ig = np.argmin(np.abs(f - f_max_ig))

    f_min_ms = 0.1
    f_max_ms = 0.5
    f1_ms = np.argmin(np.abs(f - f_min_ms))
    f2_ms = np.argmin(np.abs(f - f_max_ms))

    # Plot the spectrogram and frequency band analyses
    plt.figure(dpi=300, figsize=(25, 20))
    plt.suptitle(f"{raw_stream[0].stats.network}.{raw_stream[0].stats.station}")
    plt.subplot(121)

    # Spectrogram plot
    plt.pcolormesh(f, t2, 10 * np.log10(Sp.T), norm=norm_p)
    plt.vlines(f[f1_ig], 0, t2[-1], linewidth=3, linestyle="dashed", color="blue")
    plt.vlines(f[f2_ig], 0, t2[-1], linewidth=3, linestyle="dashed", color="blue")
    plt.vlines(f[f1_ms], 0, t2[-1], linewidth=3, linestyle="dashed", color="red")
    plt.vlines(f[f2_ms], 0, t2[-1], linewidth=3, linestyle="dashed", color="red")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Date")
    plt.xlim([0.001, 1])
    plt.yticks(tick_positions, dates, rotation=45)
    plt.colorbar(label="Pressure $(Pa^2/Hz)$[dB]", orientation="vertical")
    plt.xscale("log")

    # Frequency band analysis plots
    plt.subplot(122)
    plt.plot(
        10 * np.log10(np.mean(Sp[f1_ig:f2_ig], axis=0)),
        np.arange(0, len(Sp[0])),
        "b",
        label="Infra-Gravity",
    )
    plt.plot(
        10 * np.log10(np.mean(Sp[f1_ms:f2_ms], axis=0)),
        np.arange(0, len(Sp[0])),
        "red",
        label="Microsiesmic",
    )
    plt.ylim([0, len(t)])
    plt.xlabel("Pressure $(Pa^2/Hz)$[dB]")
    plt.yticks([])
    plt.xlim([0, 60])

    plt.legend(loc="upper left")
    plt.tight_layout()

    f, t, Sp = scipy.signal.spectrogram(
        raw_stream.select(component="Z")[0].data,
        fs=raw_stream[0].stats.sampling_rate,
        nperseg=nseg,
        noverlap=(nseg / 2),
        window="hann",
        scaling="density",
    )

    # Calculate the time difference between consecutive spectrogram points in hours
    # time_diff_hours = (t[1] - t[0]) / 3600

    tick_positions, dates = _date_ticks(
        raw_stream[0].stats.starttime, raw_stream[0].stats.endtime, len(t)
    )

    # Create a new time array for plotting
    t2 = np.arange(0, len(t))

    # Normalize pressure values for color mapping in the plot
    norm_p = mpl.colors.Normalize(vmin=0, vmax=60)

    # Define frequency bands for infra-gravity and microseismic analyses
    f_min_ig = 0.005
    f_max_ig = 0.02
    f1_ig = np.argmin(np.abs(f - f_min_ig))
    f2_ig = np.argmin(np.abs(f - f_max_ig))

    f_min_ms = 0.1
    f_max_ms = 0.5
    f1_ms = np.argmin(np.abs(f - f_min_ms))
    f2_ms = np.argmin(np.abs(f - f_max_ms))

    # Plot the spectrogram and frequency band analyses
    plt.figure(dpi=300, figsize=(25, 20))
    plt.suptitle(f"{raw_stream[0].stats.network}.{raw_stream[0].stats.station}")
    plt.subplot(121)

    # Spectrogram plot
    plt.pcolormesh(f, t2, 10 * np.log10(Sp.T))
    plt.vlines(f[f1_ig], 0, t2[-1], linewidth=3, linestyle="dashed", color="blue")
    plt.vlines(f[f2_ig], 0, t2[-1], linewidth=3, linestyle="dashed", color="blue")
    plt.vlines(f[f1_ms], 0, t2[-1], linewidth=3, linestyle="dashed", color="red")
    plt.vlines(f[f2_ms], 0, t2[-1], linewidth=3, linestyle="dashed", color="red")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Date")
    plt.xlim([0.001, 1])
    plt.yticks(tick_positions, dates, rotation=45)
    plt.colorbar(label="Acceleration $(m/s^2)^2$[dB]", orientation="vertical")
    plt.xscale("log")

    # Frequency band analysis plots
    plt.subplot(122)
    plt.plot(
        10 * np.log10(np.mean(Sp[f1_ig:f2_ig], axis=0)),
        np.arange(0, len(Sp[0])),
        "b",
        label="Infra-Gravity",
    )
    plt.plot(
        10 * np.log10(np.mean(Sp[f1_ms:f2_ms], axis=0)),
        np.arange(0, len(Sp[0])),
        "red",
        label="Microsiesmic",
    )
    plt.ylim([0, len(t)])
    plt.xlabel("Acceleration $(m/s^2)^2$[dB]")
    plt.yticks([])
    # plt.xlim([0, 60])

    plt.legend(loc="upper left")
    plt.tight_layout()


# %%


# %%
def coherogram_spectrogram_alpha(st, nseg=2**12, tw=1, f_min=0.005, f_max=0.02):

    Tresh_coh = 0.8
    Tresh_Dz = -170
    Tresh_Dp = -70

    plt.set_cmap("jet")
    norm_z = mpl.colors.Normalize(vmin=-190, vmax=-120)
    norm_p = mpl.colors.Normalize(vmin=0, vmax=60)
    norm_coh = mpl.colors.Normalize(vmin=0, vmax=1)

    TP = 5
    f, t, Sz = scipy.signal.spectrogram(
        st.select(component="Z")[0].data,
        fs=st[3].stats.sampling_rate,
        nperseg=nseg,
        noverlap=(nseg / 2),
        window="hann",
    )
    f, t, Sp = scipy.signal.spectrogram(
        st.select(component="H")[0].data,
        fs=st[0].stats.sampling_rate,
        nperseg=nseg,
        noverlap=(nseg / 2),
        window="hann",
    )

    ws = int(tw * 60 * 60 * st[0].stats.sampling_rate)
    Z, nd = sliding_window(st.select(component="Z")[0].data, ws=ws, hann=True)
    P, nd = sliding_window(st.select(component="H")[0].data, ws=ws, hann=True)

    f, Czp = scipy.signal.coherence(
        Z,
        P,
        fs=st[0].stats.sampling_rate,
        nperseg=nseg,
        noverlap=(nseg * 0.5),
        window=scipy.signal.windows.tukey(
            nseg, (TP * 60 * st[0].stats.sampling_rate) / nseg
        ),
    )

    f, Dzz = scipy.signal.welch(
        Z,
        fs=st[0].stats.sampling_rate,
        nperseg=nseg,
        noverlap=(nseg * 0.5),
        window=scipy.signal.windows.tukey(
            nseg, (TP * 60 * st[0].stats.sampling_rate) / nseg
        ),
    )

    f, Dpp = scipy.signal.welch(
        P,
        fs=st[0].stats.sampling_rate,
        nperseg=nseg,
        noverlap=(nseg * 0.5),
        window=scipy.signal.windows.tukey(
            nseg, (TP * 60 * st[0].stats.sampling_rate) / nseg
        ),
    )


    f1 = np.argmin(np.abs(f - f_min))
    f2 = np.argmin(np.abs(f - f_max))

    tick_positions, dates = _date_ticks(
        st[0].stats.starttime, st[0].stats.endtime, len(t)
    )

    t2 = np.arange(0, len(t))
    t1 = np.arange(0, len(t))

    # plt.yticks([])
    # plt.tight_layout()

    Dzz_smoothed = 10 * np.log10(
        (
            plot_smooth(
                np.median(Dzz[:, f1:f2] * (2 * np.pi * f[f1:f2]) ** 4, axis=1), 10, 1
            )
        )
    )
    Dpp_smoothed = 10 * np.log10(
        plot_smooth(np.median(Dpp[:, f1:f2], axis=1), 10, 1)
    )
    Czp_smoothed = plot_smooth(np.median(Czp[:, f1:f2], axis=1), 10, 1)

    good_windows = []
    bad_windows = []

    for i in range(0, len(Dzz)):
        if (
            10 * np.log10(np.median(Dzz[i, f1:f2] * (2 * np.pi * f[f1:f2]) ** 4))
            > Tresh_Dz
            and 10 * np.log10(np.median(Dpp[i, f1:f2])) > Tresh_Dp
            and np.median(Czp[i, f1:f2]) > Tresh_coh
        ):
            good_windows.append(i)
        else:
            bad_windows.append(i)

    tick_positions_2 = np.linspace(0, max(0, len(Dpp) - 1), len(dates))

    import matplotlib.gridspec as gridspec

    plt.figure(
        dpi=300, figsize=(40, 20)
    )  # Adjusted figsize to better suit the new layout

    gs = gridspec.GridSpec(
        1, 6, width_ratios=[3, 1, 3, 1, 3, 1]
    )  # Adjust ratios as needed
    plt.suptitle(st[0].stats.station)
    # Subplot 1 (twice the width)
    ax1 = plt.subplot(gs[0])
    # Your plotting commands for subplot 1...
    plt.pcolormesh(f, t2, 10 * np.log10(Sp.T), norm=norm_p)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Time [Date]")
    plt.xscale("log")
    plt.xlim(0.001, 1)
    plt.title("Spectrogram BDH", y=1.025)
    plt.vlines(f[f1], 0, t2[-1], linewidth=5, linestyle="dashed", color="black")
    plt.vlines(f[f2], 0, t2[-1], linewidth=5, linestyle="dashed", color="black")

    plt.yticks(tick_positions, dates, rotation=60)

    # plt.colorbar()
    plt.colorbar(label="Pressure $(Pa^2/Hz)$[dB]", orientation="vertical")
    plt.tight_layout()
    # Subplot 2 (half the width of subplot 1)
    ax2 = plt.subplot(gs[1])
    # Your plotting commands for subplot 2...
    f1 = np.argmin(np.abs(f - f_min))
    f2 = np.argmin(np.abs(f - f_max))
    plt.plot(
        10 * np.log10(np.median(Dpp[:, f1:f2], axis=1)), np.arange(0, len(Czp)), "b"
    )
    plt.plot(
        10
        * np.log10(plot_smooth(np.median(Dpp[:, f1:f2], axis=1), 10, 1)),
        np.arange(0, len(Czp)),
        "black",
        linewidth=5,
    )

    plt.xlabel("$Pa^2/Hz$ [dB]")
    plt.grid(True)
    plt.ylim([0, len(Czp)])
    plt.tight_layout()
    # plt.xlim([20,40])
    plt.yticks([])
    # Subplot 3 (twice the width, next to subplot 2)
    ax3 = plt.subplot(gs[2])
    # Your plotting commands for subplot 3...
    plt.pcolormesh(f, t2, 10 * np.log10(Sz.T * (2 * np.pi * f) ** 4))
    plt.xlabel("Frequency [Hz]")
    plt.xscale("log")
    plt.xlim(0.001, 1)
    plt.title("Spectrogram BHZ", y=1.025)
    plt.yticks([])
    plt.colorbar(label="Acceleration$((m/s^2)^2/Hz)$ [dB] ", orientation="vertical")
    plt.vlines(f[f1], 0, t2[-1], linewidth=5, linestyle="dashed", color="black")
    plt.vlines(f[f2], 0, t2[-1], linewidth=5, linestyle="dashed", color="black")

    # Subplot 4 (half the width of subplot 3)
    ax4 = plt.subplot(gs[3])
    # Your plotting commands for subplot 4...
    plt.plot(
        10 * np.log10((np.median(Dzz[:, f1:f2] * (2 * np.pi * f[f1:f2]) ** 4, axis=1))),
        np.arange(0, len(Czp)),
        "b",
    )
    plt.plot(
        10
        * np.log10(
            (
                plot_smooth(
                    np.median(Dzz[:, f1:f2] * (2 * np.pi * f[f1:f2]) ** 4, axis=1),
                    10,
                    1,
                )
            )
        ),
        np.arange(0, len(Czp)),
        "black",
        linewidth=5,
    )
    plt.xlabel("$(m/s^2)^2/Hz$ [dB]")
    plt.grid(True)
    plt.ylim([0, len(Czp)])
    plt.yticks([])
    plt.tight_layout()
    # plt.xlim([-180,-140])

    # Subplot 5 (twice the width, next to subplot 4)
    ax5 = plt.subplot(gs[4])
    # Your plotting commands for subplot 5...
    plt.pcolormesh(f, np.arange(0, len(Czp)), Czp, norm=norm_coh)
    plt.xscale("log")
    plt.xlim(0.001, 1)
    plt.xlabel("Frequency [Hz]")
    plt.xscale("log")
    plt.title("Coherogram", y=1.025)
    cbar = plt.colorbar(label="Coherence", orientation="vertical")
    plt.yticks([])
    plt.tight_layout()
    cbar.set_ticks([round(tick, 1) for tick in cbar.get_ticks()])
    plt.vlines(f[f1], 0, nd, linewidth=5, linestyle="dashed", color="black")
    plt.vlines(f[f2], 0, nd, linewidth=5, linestyle="dashed", color="black")
    plt.ylim([0, nd])
    # Subplot 6 (half the width of subplot 5)
    ax6 = plt.subplot(gs[5])
    # Your plotting commands for subplot 6...
    f1 = np.argmin(np.abs(f - f_min))
    f2 = np.argmin(np.abs(f - f_max))
    plt.plot(np.median(Czp[:, f1:f2], axis=1), np.arange(0, len(Czp)), "b")
    plt.plot(
        plot_smooth(np.median(Czp[:, f1:f2], axis=1), 10, 1),
        np.arange(0, len(Czp)),
        "black",
        linewidth=5,
    )
    plt.vlines(
        x=Tresh_coh,
        ymin=0,
        ymax=len(Czp),
        linestyles="dashed",
        color="r",
        label="0.80 Threshold",
        linewidth=5,
    )
    for ii in range(0, len(good_windows)):
        plt.axhspan(good_windows[ii], good_windows[ii] + 1, color="green", alpha=0.5)
    for ii in range(0, len(bad_windows)):
        plt.axhspan(bad_windows[ii], bad_windows[ii] + 1, color="lightcoral", alpha=0.5)
    plt.xlabel("Coherency")
    plt.grid(True)
    plt.ylim([0, len(Czp)])
    # plt.xlim([0.5,1])
    plt.yticks([])
    plt.subplots_adjust(
        left=0.05, right=0.95, top=0.93, bottom=0.07, wspace=0.02, hspace=0.2
    )
    plt.tight_layout(pad=1.0, w_pad=0, h_pad=2.0)


# %%
from scipy.signal import stft


def plot_stft(stream, nperseg=2**14):
    fs = stream[0].stats.sampling_rate

    """
  Performs a Short-Time Fourier Transform (STFT) on the provided signal and plots the results.

  Parameters:
  - signal: The input signal (a NumPy array).
  - fs: Sampling frequency of the signal (default is 1000Hz).
  - nperseg: Length of each segment (default is 256).
  - noverlap: Number of points to overlap between segments (default is None, which defaults to nperseg // 2).
  - nfft: Number of points in the FFT (default is None, which defaults to nperseg).
  - figsize: Tuple indicating the size of the plot (default is (10, 6)).
  """
    # Calculate the STFT
    f, t, Zxx = stft(stream[0], fs, nperseg=nperseg)
    norm_p = mpl.colors.Normalize(vmin=0, vmax=1)

    # Plotting
    plt.figure(dpi=300, figsize=(20, 15))
    plt.pcolormesh(t, f, np.abs(Zxx), norm=norm_p)
    plt.colorbar(label="Magnitude")
    plt.ylabel("Frequency [Hz]")
    plt.xlabel("Time [sec]")
    plt.yscale("log")
    plt.title("STFT Magnitude")
    plt.ylim(0.001, 1)
    plt.show()


# %%
def phase_frequency(stream, n=2**14):

    fft_result_z = scipy.fft.fft(stream.select(channel="*Z")[0].data, n=n)
    fft_result_p = scipy.fft.fft(stream.select(channel="*H")[0].data, n=n)
    fft_result_h1 = scipy.fft.fft(stream.select(channel="*1")[0].data, n=n)
    fft_result_h2 = scipy.fft.fft(stream.select(channel="*2")[0].data, n=n)

    # Phase Calculation: Compute the phase for each frequency component
    phase_z = np.angle(fft_result_z, deg=True)
    phase_p = np.angle(fft_result_p, deg=True)
    phase_h1 = np.angle(fft_result_h1, deg=True)
    phase_h2 = np.angle(fft_result_h2, deg=True)

    # Frequency axis (for plotting)
    sampling_rate = stream[0].stats.sampling_rate  # Example: 1000 samples per second

    freq = np.linspace(0, sampling_rate / 2, n // 2)

    # Plot the phase spectrum
    plt.figure(dpi=300, figsize=(35, 25))
    # First subplot
    plt.suptitle(
        "Phase Spectrum--"
        + str(stream[0].stats.station)
        + "--"
        + str(stream[0].stats.starttime)[0:10]
        + "--"
        + str(stream[0].stats.endtime)[0:10]
    )
    plt.subplot(2, 2, 1)  # (rows, columns, subplot number)
    plt.plot(
        freq, phase_z[: n // 2], "o", color="black", markersize=3
    )  # Smaller markers
    plt.xlim(0.001, 1)
    plt.title("Vertical")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Phase [°]")
    plt.xscale("log")
    plt.grid(True)

    # Second subplot
    plt.subplot(2, 2, 2)
    plt.plot(
        freq, phase_h1[: n // 2], "o", color="black", markersize=3
    )  # Smaller markers
    plt.xlim(0.001, 1)
    plt.title("Horizontal 1")
    plt.xscale("log")
    plt.grid(True)

    # Third subplot
    plt.subplot(2, 2, 3)
    plt.plot(
        freq, phase_p[: n // 2], "o", color="black", markersize=3
    )  # Smaller markers
    plt.xlim(0.001, 1)
    plt.title("Pressure")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Phase [°]")
    plt.xscale("log")
    plt.grid(True)

    # Fourth subplot
    plt.subplot(2, 2, 4)
    plt.plot(
        freq, phase_h2[: n // 2], "o", color="black", markersize=3
    )  # Smaller markers
    plt.xlim(0.001, 1)
    plt.title("Horizontal 2")
    plt.xlabel("Frequency [Hz]")
    plt.xscale("log")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
