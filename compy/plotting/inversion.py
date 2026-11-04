"""Figures for Metropolis inversion results: posteriors, misfits and histograms.

The layouts reproduce the figures of Aminian et al. (2025); heat-map rows
use the fixed depth sampling of those figures (``LEGACY_HEATMAP_ROW_M``).
"""

import matplotlib.pyplot as plt
import numpy as np

from compy._deprecation import renamed
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from compy.inversion.dram import serpentinization
from compy.inversion.models import (
    reference_models,
    station_starting_model,
    station_starting_model_mean,
)

#: Depth of one row of the legacy posterior heat maps [m]: the 1-m Vs profiles
#: of the first 10 km are averaged over 100 samples (docs/INVERSION.md, F1).
LEGACY_HEATMAP_ROW_M = 100.0


#: Depth of one element of the 100-element reference arrays of
#: ``start_model_plot``/``start_model_plot_mean``/``reference_models`` [m]
#: (the CRUST1.0 Moho at element 42 lies at about 7 km).
LEGACY_REFERENCE_ROW_M = 10000.0 / 60.0


_LEGACY_REFERENCE_ROWS = np.arange(100) * LEGACY_REFERENCE_ROW_M / LEGACY_HEATMAP_ROW_M


def _legacy_rows(depth_m):
    """Heat-map row index of a depth [m] in the legacy figures."""
    return depth_m / LEGACY_HEATMAP_ROW_M


def _legacy_depth_ticks(max_depth_m=14000.0, step_m=2000.0):
    """Row positions and depth labels [m, negative down] for legacy heat maps.

    Before 2.4 the ticks put row 12 at 2000 m (167-m rows) while each row is
    100 m, so features were labelled 5/3 too deep (docs/INVERSION.md, F1).
    """
    depths = np.arange(0.0, max_depth_m + step_m / 2, step_m)
    return list(depths / LEGACY_HEATMAP_ROW_M), [int(-d) for d in depths]


def plot_inversion_summary(
    starting_model,
    vs,
    mis_fit,
    ncompl,
    Data,
    likelihood_data,
    freq,
    sta,
    iteration,
    s,
    sigma_v,
    sigma_h,
    n_layer,
    alpha,
    burnin=50000,
    mis_fit_trsh=1,
):
    """Summary of one Metropolis run: data fit, misfit and likelihood traces, and the Vs profiles as an image."""

    depth = np.arange(0, -int(np.sum(starting_model[0:-1, 0, 0])), -1)

    # for i in range(0, vs.shape[0]):
    Vs_final = np.zeros(vs[0].shape)
    N = 0
    vs_filtered = []

    for i in range(burnin, iteration):
        if mis_fit[0][i] < mis_fit_trsh:
            # if starting_model[layer_number][3][i] < 50000:
            Vs_final = Vs_final + vs[i]
            N = N + 1
            vs_filtered.append(vs[i])
    vs_filtered = np.array(vs_filtered)

    Vs_final = Vs_final / N

    plt.figure(dpi=300, figsize=(15, 25))
    plt.plot(Vs_final, depth, color="black", label="Final Result", linewidth=3)

    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m Below Seafloor]")
    plt.ylim([-12000, 0])
    plt.grid(True)
    plt.legend(loc="lower left")
    plt.tight_layout()

    dd = np.zeros([len(vs[0]), 1000])

    for i in range(0, len(vs[0])):
        dd[i] = np.histogram(vs_filtered[:, i], bins=5000, range=([0, 5000]))[0]
        dd[i] = dd[i] / np.max(dd[i])
    plt.rcParams.update({"font.size": 40})
    nn = (
        int((iteration - burnin) / 1000) + 1
    )  # I want to see just 1000 points of the data
    plt.figure(dpi=300, figsize=(35, 25))
    plt.suptitle(
        "Sigma V = "
        + str(sigma_v)
        + ", Sigma h = "
        + str(sigma_h)
        + ", Alpha = "
        + str(alpha)
        + ",N of Layer ="
        + str(n_layer)
    )
    plt.subplot(121)

    for i in range(burnin, ncompl.shape[0], nn):
        if mis_fit[0][i] < mis_fit_trsh:
            plt.plot(freq, ncompl[i], color="green", linewidth=0.25)

        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Normalized Compliance")

    plt.errorbar(
        freq,
        Data,
        yerr=s,
        ecolor=("black"),
        color="black",
        fmt="none",
        linewidth=5,
        label="Measured Compliance",
        capsize=15,
    )

    # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
    #             label='Median of Brun-in')
    plt.plot(
        freq,
        ncompl[1],
        color="black",
        label="Start Compliance",
        linewidth=5,
        linestyle="dashed",
    )
    plt.ylim([10e-13, 10e-11])

    plt.grid(True)
    plt.legend(loc="upper left", fontsize=30)

    plt.subplot(122)
    plt.title(sta)
    plt.imshow(dd, aspect="auto")
    plt.colorbar()

    # Add labels, legend, and colorbar as in your code
    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m Below Seafloor]")
    plt.ylim([12000, 0])
    plt.grid(True)
    plt.legend(loc="lower left")
    plt.tight_layout()

    plt.figure(dpi=300, figsize=(16, 24))
    plt.subplot(211)
    plt.plot(
        likelihood_data[0, 1 : iteration - 1], color="green", label="Data Likelihood"
    )
    plt.xscale("log")
    plt.vlines(
        x=burnin, ymin=0, ymax=1, color="r", label="Burn-in Region", linestyles="dashed"
    )
    plt.legend(loc="upper left")

    plt.title("Likelihood")
    plt.xlabel("Iteration")
    plt.ylabel("Liklihood")

    plt.grid(True)

    plt.subplot(212)
    plt.plot(mis_fit[0, 1 : iteration - 1])
    plt.xscale("log")

    plt.vlines(
        x=burnin,
        ymin=0,
        ymax=np.max(mis_fit[0, 1 : iteration - 1]),
        color="r",
        label="Burn-in Region",
        linestyles="dashed",
        linewidth=3,
    )

    plt.title("Data Misfit")
    plt.xlabel("Iteration")
    plt.ylabel("$\\chi^2$")
    plt.grid(True)
    plt.minorticks_on()

    plt.tight_layout()

    plt.legend(loc="lower left")
    plt.tight_layout()


def plot_inversion(
    starting_model,
    vs,
    mis_fit,
    ncompl,
    Data,
    likelihood_data,
    freq,
    sta,
    iteration,
    s,
    sigma_v,
    sigma_h,
    n_layer,
    alpha,
    burnin=50000,
    mis_fit_trsh=1,
):
    """Summary of one Metropolis run: data fit, misfit and likelihood traces with burn-in, and Vs profiles."""
    linewidth_models = 0.05
    linewidth_compliance = 0.5
    depth = np.arange(0, -int(np.sum(starting_model[0:-1, 0, 0])), -1)

    # for i in range(0, vs.shape[0]):
    Vs_final = np.zeros(vs[0].shape)
    N = 0
    for i in range(burnin, iteration):
        if mis_fit[0][i] < mis_fit_trsh:
            # if starting_model[layer_number][3][i] < 50000:
            Vs_final = Vs_final + vs[i]
            N = N + 1

    Vs_final = Vs_final / N

    plt.figure(dpi=300, figsize=(15, 25))
    plt.title(sta)

    plt.plot(Vs_final, depth, color="black", label="Final Result", linewidth=3)

    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m Below Seafloor]")
    plt.ylim([-12000, 0])
    plt.grid(True)
    plt.legend(loc="lower left")
    plt.tight_layout()

    # for i in range(burnin, vs.shape[0], nn):
    #     if mis_fit[0][i] < mis_fit_trsh:

    plt.rcParams.update({"font.size": 40})
    nn = int((iteration - burnin) / 1000)  # I want to see just 1000 points of the data
    plt.figure(dpi=300, figsize=(35, 25))
    plt.suptitle(
        "Sigma V = "
        + str(sigma_v)
        + ", Sigma h = "
        + str(sigma_h)
        + ", Alpha = "
        + str(alpha)
        + ",N of Layer ="
        + str(n_layer)
    )
    plt.subplot(121)

    for i in range(burnin, ncompl.shape[0], nn):
        if mis_fit[0][i] < mis_fit_trsh:
            plt.plot(freq, ncompl[i], color="green", linewidth=linewidth_compliance)

        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Normalized Compliance")

    plt.errorbar(
        freq,
        Data,
        yerr=s,
        ecolor=("black"),
        color="black",
        fmt="none",
        linewidth=5,
        label="Measured Compliance",
        capsize=15,
    )

    # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
    #             label='Median of Brun-in')
    plt.plot(
        freq,
        ncompl[1],
        color="black",
        label="Start Compliance",
        linewidth=5,
        linestyle="dashed",
    )

    plt.grid(True)
    plt.legend(loc="upper left", fontsize=30)

    plt.subplot(122)
    plt.title("YV." + sta)
    for i in range(burnin, vs.shape[0]):
        if mis_fit[0][i] < mis_fit_trsh:
            plt.plot(vs[i], depth, color="green", linewidth=linewidth_models)
    plt.plot(
        vs[0],
        depth,
        color="black",
        label="Start Model",
        linewidth=5,
        linestyle="dashed",
    )

    # Add labels, legend, and colorbar as in your code
    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m Below Seafloor]")
    plt.ylim([-12000, 0])
    plt.grid(True)
    plt.legend(loc="lower left")
    plt.tight_layout()

    # for i in range(burnin, vs.shape[0], nn):
    #     if mis_fit[0][i] < mis_fit_trsh:

    plt.figure(dpi=300, figsize=(16, 24))

    plt.subplot(211)
    plt.suptitle("YV." + sta)

    plt.plot(
        likelihood_data[0, 1 : iteration - 1], color="green", label="Data Likelihood"
    )
    plt.xscale("log")
    plt.vlines(
        x=burnin, ymin=0, ymax=1, color="r", label="Burn-in Region", linestyles="dashed"
    )
    plt.legend(loc="upper left")

    plt.title("Likelihood")
    plt.xlabel("Iteration")
    plt.ylabel("Liklihood")

    plt.grid(True)

    plt.subplot(212)
    plt.plot(mis_fit[0, 1 : iteration - 1])
    plt.xscale("log")

    plt.vlines(
        x=burnin,
        ymin=0,
        ymax=np.max(mis_fit[0, 1 : iteration - 1]),
        color="r",
        label="Burn-in Region",
        linestyles="dashed",
        linewidth=3,
    )

    plt.title("Data Misfit")
    plt.xlabel("Iteration")
    plt.ylabel("$\\chi^2$")
    plt.grid(True)
    plt.minorticks_on()

    plt.tight_layout()

    plt.legend(loc="lower left")
    plt.tight_layout()


def plot_inversion_with_start(
    starting_model,
    vs,
    vs0,
    mis_fit,
    ncompl,
    Data,
    likelihood_data,
    freq,
    sta,
    iteration,
    s,
    sigma_v,
    sigma_h,
    n_layer,
    alpha,
    burnin=50000,
    mis_fit_trsh=1,
):
    """As ``plot_inversion``, also drawing the starting profile ``vs0``."""
    linewidth_models = 0.05
    linewidth_compliance = 0.5
    depth = np.arange(0, -int(np.sum(starting_model[0:-1, 0, 0])), -1)

    # for i in range(0, vs.shape[0]):
    Vs_final = np.zeros(vs[0].shape)
    N = 0
    for i in range(burnin, iteration):
        if mis_fit[0][i] < mis_fit_trsh:
            # if starting_model[layer_number][3][i] < 50000:
            Vs_final = Vs_final + vs[i]
            N = N + 1

    Vs_final = Vs_final / N

    plt.figure(dpi=300, figsize=(15, 25))
    plt.title(sta)

    plt.plot(Vs_final, depth, color="black", label="Final Result", linewidth=3)

    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m Below Seafloor]")
    plt.ylim([-12000, 0])
    plt.xlim([0, 5000])
    plt.grid(True)
    plt.legend(loc="lower left")
    plt.tight_layout()

    # for i in range(burnin, vs.shape[0], nn):
    #     if mis_fit[0][i] < mis_fit_trsh:

    plt.rcParams.update({"font.size": 40})

    nn = int((iteration - burnin) / 1000)  # I want to see just 1000 points of the data
    plt.figure(dpi=300, figsize=(35, 25))
    plt.suptitle(
        "Sigma V = "
        + str(sigma_v)
        + ", Sigma h = "
        + str(sigma_h)
        + ", Alpha = "
        + str(alpha)
        + ",N of Layer ="
        + str(n_layer)
    )
    plt.subplot(121)

    for i in range(burnin, ncompl.shape[0], nn):
        if mis_fit[0][i] < mis_fit_trsh:
            if np.mean(vs[i]) < np.mean(vs0[:, 0 : len(vs[0])]):
                plt.plot(freq, ncompl[i], color="red", linewidth=linewidth_compliance)
            else:
                plt.plot(freq, ncompl[i], color="green", linewidth=linewidth_compliance)

        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Normalized Compliance")

    plt.errorbar(
        freq,
        Data,
        yerr=s,
        ecolor=("black"),
        color="black",
        fmt="none",
        linewidth=5,
        label="Measured Compliance",
        capsize=15,
    )

    # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
    #             label='Median of Brun-in')
    plt.plot(
        freq,
        ncompl[1],
        color="black",
        label="Start Compliance",
        linewidth=5,
        linestyle="dashed",
    )

    plt.grid(True)
    plt.legend(loc="upper left", fontsize=30)

    plt.subplot(122)
    plt.title("YV." + sta)
    for i in range(burnin, vs.shape[0]):
        if mis_fit[0][i] < mis_fit_trsh:
            if np.mean(vs[i]) < np.mean(vs0[:, 0 : len(vs[0])]):
                plt.plot(vs[i], depth, color="red", linewidth=linewidth_models)
            else:
                plt.plot(vs[i], depth, color="green", linewidth=linewidth_models)

    plt.plot(
        vs[0],
        depth,
        color="black",
        label="Start Model",
        linewidth=5,
        linestyle="dashed",
    )

    # Add labels, legend, and colorbar as in your code
    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m Below Seafloor]")
    plt.ylim([-12000, 0])
    plt.grid(True)
    plt.legend(loc="lower left")
    plt.tight_layout()

    # for i in range(burnin, vs.shape[0], nn):
    #     if mis_fit[0][i] < mis_fit_trsh:

    plt.figure(dpi=300, figsize=(16, 24))

    plt.subplot(211)
    plt.suptitle("YV." + sta)

    plt.plot(
        likelihood_data[0, 1 : iteration - 1], color="green", label="Data Likelihood"
    )
    plt.xscale("log")
    plt.vlines(
        x=burnin, ymin=0, ymax=1, color="r", label="Burn-in Region", linestyles="dashed"
    )
    plt.legend(loc="upper left")

    plt.title("Likelihood")
    plt.xlabel("Iteration")
    plt.ylabel("Liklihood")

    plt.grid(True)

    plt.subplot(212)
    plt.plot(mis_fit[0, 1 : iteration - 1])
    plt.xscale("log")

    plt.vlines(
        x=burnin,
        ymin=0,
        ymax=np.max(mis_fit[0, 1 : iteration - 1]),
        color="r",
        label="Burn-in Region",
        linestyles="dashed",
        linewidth=3,
    )

    plt.title("Data Misfit")
    plt.xlabel("Iteration")
    plt.ylabel("$\\chi^2$")
    plt.grid(True)
    plt.minorticks_on()

    plt.tight_layout()

    plt.legend(loc="lower left")
    plt.tight_layout()


def plot_inversion_density(
    vs, vs0, mis_fit, Data, s, freq, sta, burnin, ncompl, iteration, mis_fit_trsh=8
):
    """Posterior Vs probability versus depth for one station, with the measured and predicted compliance."""
    bins = 100
    jj = 0
    vs_good = []

    for ii in range(0, len(mis_fit[0])):
        if mis_fit[0][ii] < mis_fit_trsh:
            # if vs[ii][9000][0] < vs[ii][6000][0]:

            vs_good.append(vs[ii])
            jj = jj + 1

    vs_good = np.array(vs_good)

    linewidth_compliance = 0.5
    nn = int((iteration - burnin) / 1000)

    # for i in range(0,len(vs[0])):

    # # Setting custom y-axis ticks and labels

    vs_good = vs_good[:, 0:10000, 0]

    downsample_factor = vs_good.shape[1] // 100

    # Reshape the array to prepare for averaging
    # New shape will be (20000, 100, downsample_factor, 1)
    vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

    # Take the mean along the downsample_factor dimension
    vs_downsampled = vs_reshaped.mean(axis=2)

    c = np.zeros([bins, bins])
    b = plt.hist(vs_reshaped[:, 10, 0], bins=bins, range=([0, 5000]))[1]

    for i in range(0, len(vs_downsampled[0])):
        c[i] = plt.hist(
            vs_downsampled[:, i], bins=bins, range=([0, 5000]), density=True, log=False
        )[0]
        c[i] = c[i] / np.max(c[i])

    from matplotlib.colors import LinearSegmentedColormap

    # Define the colors for the colormap (from white to red to black)
    colors = ["white", "red", "black"]

    # Define the transition points for the colors
    n_bins = [0, 0.4, 1]  # You can adjust these thresholds based on your data

    # Create the custom colormap
    custom_colormap = LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(n_bins, colors))
    )

    plt.figure(dpi=300, figsize=(30, 20))
    plt.suptitle("YV." + str(sta))
    plt.subplot(121)
    for i in range(burnin, ncompl.shape[0], nn):
        if mis_fit[0][i] < mis_fit_trsh:
            plt.plot(freq, ncompl[i], color="green", linewidth=linewidth_compliance)

    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Normalized Compliance")

    plt.errorbar(
        freq,
        Data,
        yerr=s,
        ecolor=("black"),
        color="black",
        fmt="none",
        linewidth=5,
        label="Measured Compliance",
        capsize=15,
    )
    plt.legend(loc="upper left", fontsize=30)

    plt.subplot(122)
    selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
    selected_labels = b[selected_indices]
    selected_labels = selected_labels.astype(int)

    ff = plt.imshow(
        c, aspect="auto", cmap=custom_colormap, norm=plt.Normalize(vmin=0, vmax=1)
    )

    # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
    #             label='Median of Brun-in')

    plt.grid(True)

    # Setting custom y-axis ticks and labels
    plt.xticks(ticks=selected_indices, labels=selected_labels)
    cbar = plt.colorbar(ff)

    selected_indices_y, selected_labels_y = _legacy_depth_ticks()
    plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
    plt.xlabel("Shear Velocity [m/s]")
    plt.ylabel("Depth [m]")
    cbar.set_label("Probability")
    plt.legend(loc="lower left", fontsize=35)
    plt.ylim(_legacy_rows(10000), 0)
    plt.tight_layout()
    plt.show()


def plot_inversion_density_all(Inversion_container):
    """Posterior Vs probability versus depth for every station in ``Inversion_container``."""

    plt.figure(dpi=300, figsize=(40, 25))

    for ii in range(0, len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        sta = Inversion_container[ii]["Station"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        bins = 100
        jj = 0
        vs_good = []

        for i in range(0, len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                # if vs[ii][9000][0] < vs[ii][6000][0]:

                vs_good.append(vs[i])
                jj = jj + 1

        vs_good = np.array(vs_good)

        vs_good = vs_good[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100

        # Reshape the array to prepare for averaging
        # New shape will be (20000, 100, downsample_factor, 1)
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

        # Take the mean along the downsample_factor dimension
        vs_downsampled = vs_reshaped.mean(axis=2)

        c = np.zeros([bins, bins])
        b = np.histogram(vs_reshaped[:, 10, 0], bins=bins, range=(0, 5000))[1]

        for i in range(0, len(vs_downsampled[0])):
            c[i] = np.histogram(
                vs_downsampled[:, i], bins=bins, range=(0, 5000), density=True
            )[0]
            c[i] = c[i] / np.max(c[i])

        from matplotlib.colors import LinearSegmentedColormap

        # Define the colors for the colormap (from white to red to black)
        colors = ["white", "red", "black"]
        viridis = plt.cm.get_cmap("viridis", 256)

        # Sample the colormap
        colors = [
            viridis(0.0),  # dark blue/lilac at the start of the colormap
            viridis(0.5),  # greenish in the middle of the colormap
            viridis(1.0),
        ]  # bright yellow at the end of the colormap

        # Define the transition points for the colors
        n_bins = [0, 0.4, 1]  # You can adjust these thresholds based on your data

        # Create the custom colormap
        custom_colormap = LinearSegmentedColormap.from_list(
            "custom_cmap", list(zip(n_bins, colors))
        )

        plt.subplot(2, 4, int(ii + 1))
        plt.title(str("YV.") + sta)
        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        ff = plt.imshow(
            c, aspect="auto", cmap=custom_colormap, norm=plt.Normalize(vmin=0, vmax=1)
        )
        plt.plot(
            np.median(vs_downsampled, axis=0) / 50,
            depth,
            color="red",
            linewidth=8,
            linestyle="dashed",
            label="Median",
        )

        # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
        #             label='Median of Brun-in')

        plt.grid(True)

        # Setting custom y-axis ticks and labels
        plt.xticks(ticks=selected_indices, labels=selected_labels)
        cbar = plt.colorbar(ff)

        selected_indices_y, selected_labels_y = _legacy_depth_ticks()
        plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
        if ii == 0 or ii == 4:
            plt.ylabel("Depth [m]")
        if ii == 4 or ii == 5 or ii == 6 or ii == 7:
            plt.xlabel("Shear Velocity [m/s]")
        if ii == 3 or ii == 7:
            cbar.set_label("Probability")
        if ii == 7:
            plt.legend(loc="lower left", fontsize=30)

        if ii == 0:
            plt.text(
                0.5, 0.95, "a)", fontsize=40, fontweight="bold", va="top", color="white"
            )
        if ii == 1:
            plt.text(
                0.05,
                0.95,
                "b)",
                fontsize=40,
                fontweight="bold",
                va="top",
                color="white",
            )
        if ii == 2:
            plt.text(
                0.001,
                0.95,
                "c)",
                fontsize=40,
                fontweight="bold",
                va="top",
                color="white",
            )
        if ii == 3:
            plt.text(
                0.4, 0.95, "d)", fontsize=40, fontweight="bold", va="top", color="white"
            )
        if ii == 4:
            plt.text(
                0.9, 0.95, "e)", fontsize=40, fontweight="bold", va="top", color="white"
            )
        if ii == 5:
            plt.text(
                0.1, 0.95, "f)", fontsize=40, fontweight="bold", va="top", color="white"
            )
        if ii == 6:
            plt.text(
                0.1, 0.95, "g)", fontsize=40, fontweight="bold", va="top", color="white"
            )
        if ii == 7:
            plt.text(
                0.1, 0.95, "h)", fontsize=40, fontweight="bold", va="top", color="white"
            )

        plt.ylim(_legacy_rows(8000), 0)
    plt.tight_layout()


def plot_publication_figure(Inversion_container, *, image_dir):
    """Publication layout; image_dir contains the six named bathymetry images."""
    from pathlib import Path

    image_dir = Path(image_dir)
    import matplotlib.image as mpimg  # Import the image module
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    # Create a new figure with specified dimensions
    fig = plt.figure(figsize=(35, 30), constrained_layout=True)

    # Define the grid layout with customized row heights for the third column
    gs = gridspec.GridSpec(3, 3, width_ratios=[1, 1, 2], height_ratios=[1, 1, 1])

    img1 = mpimg.imread(image_dir / "RR36.png")
    img2 = mpimg.imread(image_dir / "RR38.png")
    img3 = mpimg.imread(image_dir / "RR40.png")
    img4 = mpimg.imread(image_dir / "RR50.png")
    img5 = mpimg.imread(image_dir / "RR52.png")
    img6 = mpimg.imread(image_dir / "Rift_valley.png")

    # First plot
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title("RR36")
    ax1.imshow(img1, extent=[0, 100, 0, 100], aspect="auto")
    ax1.set_xlabel("Distance East/West [km]")
    ax1.set_ylabel("Distance North/South [km]")
    plt.text(
        0.05,
        0.95,
        "a)",
        transform=ax1.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )

    # Second plot
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_title("RR38")
    ax2.imshow(img2, extent=[0, 100, 0, 100], aspect="auto")
    ax2.set_xlabel("Distance East/West [km]")
    ax2.set_ylabel("Distance North/South [km]")
    plt.text(
        0.05,
        0.95,
        "b)",
        transform=ax2.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )

    # Third plot
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.set_title("RR40")
    ax3.imshow(img3, extent=[0, 100, 0, 100], aspect="auto")
    ax3.set_xlabel("Distance East/West [km]")
    ax3.set_ylabel("Distance North/South [km]")
    plt.text(
        0.05,
        0.95,
        "c)",
        transform=ax3.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )

    # Fourth plot
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.set_title("RR50")
    ax4.imshow(img4, extent=[0, 100, 0, 100], aspect="auto")
    ax4.set_xlabel("Distance East/West [km]")
    ax4.set_ylabel("Distance North/South [km]")
    plt.text(
        0.05,
        0.95,
        "d)",
        transform=ax4.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )

    # Fifth plot
    ax5 = fig.add_subplot(gs[2, 0])
    ax5.set_title("RR52")
    ax5.imshow(img5, extent=[0, 100, 0, 100], aspect="auto")
    ax5.set_xlabel("Distance East/West [km]")
    ax5.set_ylabel("Distance North/South [km]")
    plt.text(
        0.05,
        0.95,
        "e)",
        transform=ax5.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )

    # Sixth plot
    ax6 = fig.add_subplot(gs[2, 1])
    ax6.imshow(img6, extent=[-25, 25, -5, 5], aspect="auto")

    ax6.set_xlabel("Distance [km]")
    ax6.set_ylabel("Depth [km]")
    plt.text(
        0.05,
        0.95,
        "f)",
        transform=ax6.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )
    ax6.grid(False)

    # Seventh plot
    ax7 = fig.add_subplot(gs[0:3, 2])  # Spanning the top two rows of the third column

    for ii in range(0, len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        sta = Inversion_container[ii]["Station"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        bins = 100
        jj = 0
        vs_good = []

        for i in range(0, len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                # if vs[ii][9000][0] < vs[ii][6000][0]:

                vs_good.append(vs[i])
                jj = jj + 1

        vs_good = np.array(vs_good)

        vs_good = vs_good[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100

        # Reshape the array to prepare for averaging
        # New shape will be (20000, 100, downsample_factor, 1)
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

        # Take the mean along the downsample_factor dimension
        vs_downsampled = vs_reshaped.mean(axis=2)

        c = np.zeros([bins, bins])
        b = np.histogram(vs_reshaped[:, 10, 0], bins=bins, range=(0, 5000))[1]

        for i in range(0, len(vs_downsampled[0])):
            c[i] = np.histogram(
                vs_downsampled[:, i], bins=bins, range=(0, 5000), density=True
            )[0]
            c[i] = c[i] / np.max(c[i])

        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        plt.vlines(90, 0, 100, linestyles="dashed", colors="grey", linewidth=3)
        plt.text(
            95,
            25,
            "Fresh peridotite",
            fontsize=35,
            color="grey",
            ha="center",
            rotation=90,
        )

        plt.plot(
            np.median(vs_downsampled, axis=0) / 50,
            depth,
            linewidth=5,
            linestyle="solid",
            label=sta,
        )

        plt.grid(True)

        # Setting custom y-axis ticks and labels
        plt.xticks(ticks=selected_indices, labels=selected_labels)

        selected_indices_y, selected_labels_y = _legacy_depth_ticks()
        plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
        plt.legend(loc="upper right")
        plt.ylim(_legacy_rows(8000), 0)

        ax7.set_xlabel("Shear Velocity [m/s]")
        ax7.set_ylabel("Depth [m]")
        plt.text(
            0.025,
            0.98,
            "g)",
            transform=ax7.transAxes,
            fontsize=40,
            fontweight="bold",
            va="top",
        )

    # Eighth plot as inset
    ax8 = inset_axes(ax7, width="40%", height="40%", loc="lower left", borderpad=3.5)

    for ii in range(0, len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        sta = Inversion_container[ii]["Station"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        bins = 100
        jj = 0
        vs_good = []

        for i in range(0, len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                # if vs[ii][9000][0] < vs[ii][6000][0]:

                vs_good.append(vs[i])
                jj = jj + 1

        vs_good = np.array(vs_good)

        vs_good = vs_good[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100

        # Reshape the array to prepare for averaging
        # New shape will be (20000, 100, downsample_factor, 1)
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

        # Take the mean along the downsample_factor dimension
        vs_downsampled = vs_reshaped.mean(axis=2)

        c = np.zeros([bins, bins])
        b = np.histogram(vs_reshaped[:, 10, 0], bins=bins, range=(0, 5000))[1]

        for i in range(0, len(vs_downsampled[0])):
            c[i] = np.histogram(
                vs_downsampled[:, i], bins=bins, range=(0, 5000), density=True
            )[0]
            c[i] = c[i] / np.max(c[i])
        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        # Fig. 10c relation (4.2 km/s: 0 %, 2.2 km/s: 100 %); was (9000 - 2 Vs)/30 (F2).
        serpentinization_percentage = serpentinization(
            np.median(vs_downsampled, axis=0)
        )
        plt.plot(
            serpentinization_percentage,
            depth,
            linewidth=5,
            linestyle="solid",
            label=sta,
        )

        plt.grid(True)
        plt.xlim([0, 100])
        selected_indices_y, selected_labels_y = _legacy_depth_ticks()
        plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
        plt.ylim(_legacy_rows(8000), 0)
        plt.legend(loc="lower right")

        ax8.set_xlabel("Serpentinization (%)")
    plt.text(
        0.05,
        0.95,
        "h)",
        transform=ax8.transAxes,
        fontsize=40,
        fontweight="bold",
        va="top",
    )

    plt.tight_layout()
    plt.show()


def plot_inversion_all(Inversion_container):
    """Median Vs profiles of every station with the SWIR reference profiles."""
    reference_model = reference_models() / 50

    plt.figure(dpi=300, figsize=(15, 25))

    for ii in range(0, len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        sta = Inversion_container[ii]["Station"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        bins = 100
        jj = 0
        vs_good = []

        for i in range(0, len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                # if vs[ii][9000][0] < vs[ii][6000][0]:

                vs_good.append(vs[i])
                jj = jj + 1

        vs_good = np.array(vs_good)

        vs_good = vs_good[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100

        # Reshape the array to prepare for averaging
        # New shape will be (20000, 100, downsample_factor, 1)
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

        # Take the mean along the downsample_factor dimension
        vs_downsampled = vs_reshaped.mean(axis=2)

        c = np.zeros([bins, bins])
        b = np.histogram(vs_reshaped[:, 10, 0], bins=bins, range=(0, 5000))[1]

        for i in range(0, len(vs_downsampled[0])):
            c[i] = np.histogram(
                vs_downsampled[:, i], bins=bins, range=(0, 5000), density=True
            )[0]
            c[i] = c[i] / np.max(c[i])

        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        plt.plot(np.median(vs_downsampled, axis=0) / 50, depth, linewidth=5, label=sta)

        # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
        #             label='Median of Brun-in')

        plt.grid(True)

        # Setting custom y-axis ticks and labels
    plt.xticks(ticks=selected_indices, labels=selected_labels)

    selected_indices_y, selected_labels_y = _legacy_depth_ticks()
    plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
    plt.ylabel("Depth [m]")
    plt.xlabel("Shear Velocity [m/s]")
    plt.plot(
        reference_model[0],
        _LEGACY_REFERENCE_ROWS,
        color="red",
        linewidth=5,
        linestyle="dashed",
        label="SWIR 64°30'E (EW)",
    )
    plt.plot(
        reference_model[1],
        _LEGACY_REFERENCE_ROWS,
        color="purple",
        linewidth=5,
        linestyle="dotted",
        label="SWIR 64°30'E (NS) ",
    )
    plt.plot(
        reference_model[2],
        _LEGACY_REFERENCE_ROWS,
        color="green",
        linewidth=5,
        linestyle="dotted",
        label="SWIR 65-64°E ",
    )
    plt.legend(loc="lower left", fontsize=30)
    plt.ylim(_legacy_rows(8000), 0)
    plt.tight_layout()


def plot_inversion_density_mean_all(Inversion_container):
    """Median Vs profiles of every station against CRUST1.0, Atlantis Bank and SWIR references."""

    plt.figure(dpi=300, figsize=(40, 25))

    for ii in range(0, len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        sta = Inversion_container[ii]["Station"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        start_model = station_starting_model(sta) / 50
        reference_model = reference_models() / 50

        bins = 100
        jj = 0
        vs_good = []

        for i in range(0, len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                # if vs[ii][9000][0] < vs[ii][6000][0]:

                vs_good.append(vs[i])
                jj = jj + 1

        vs_good = np.array(vs_good)

        vs_good = vs_good[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100

        # Reshape the array to prepare for averaging
        # New shape will be (20000, 100, downsample_factor, 1)
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

        # Take the mean along the downsample_factor dimension
        vs_downsampled = vs_reshaped.mean(axis=2)

        c = np.zeros([bins, bins])
        b = np.histogram(vs_reshaped[:, 10, 0], bins=bins, range=(0, 5000))[1]

        for i in range(0, len(vs_downsampled[0])):
            c[i] = np.histogram(
                vs_downsampled[:, i], bins=bins, range=(0, 5000), density=True
            )[0]
            c[i] = c[i] / np.max(c[i])

        # Define the colors for the colormap (from white to red to black)

        # Define the transition points for the colors

        # Create the custom colormap

        plt.subplot(2, 4, int(ii + 1))
        plt.title(str("YV.") + sta)
        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        plt.plot(
            start_model,
            _LEGACY_REFERENCE_ROWS,
            color="blue",
            linewidth=5,
            linestyle="dashdot",
            label="Crust 1 ",
        )
        if sta == "RR40":
            plt.plot(
                reference_model[0],
                _LEGACY_REFERENCE_ROWS,
                color="red",
                linewidth=5,
                linestyle="dashed",
                label="SWIR 64°30'E (EW)",
            )
            plt.plot(
                reference_model[1],
                _LEGACY_REFERENCE_ROWS,
                color="purple",
                linewidth=5,
                linestyle="dotted",
                label="SWIR 64°30'E (NS) ",
            )
            plt.plot(
                reference_model[2],
                _LEGACY_REFERENCE_ROWS,
                color="green",
                linewidth=5,
                linestyle="dotted",
                label="SWIR 65-64°E ",
            )

        if sta == "RR36" or sta == "RR38":
            plt.plot(
                reference_model[4],
                _LEGACY_REFERENCE_ROWS,
                color="red",
                linewidth=5,
                linestyle="dashed",
                label="Atlantis Bank",
            )

        plt.plot(
            np.median(vs_downsampled, axis=0) / 50,
            depth,
            color="black",
            linewidth=6,
            linestyle="solid",
            label="Median (This Study)",
        )

        # plt.plot(freq, np.median(ncompl[burnin:iteration],axis=0), color='blue',
        #             label='Median of Brun-in')

        plt.grid(True)

        # Setting custom y-axis ticks and labels
        plt.xticks(ticks=selected_indices, labels=selected_labels)

        selected_indices_y, selected_labels_y = _legacy_depth_ticks()
        plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
        if ii == 0 or ii == 4:
            plt.ylabel("Depth [m]")
        if ii == 4 or ii == 5 or ii == 6 or ii == 7:
            plt.xlabel("Shear Velocity [m/s]")
        # if ii == 3 or ii == 7:
        # if ii == 7:
        plt.legend(loc="lower left", fontsize=30)
        plt.ylim(_legacy_rows(8000), 0)
    plt.tight_layout()


def plot_inversion_serpentinization(Inversion_container):
    """Vs and degree of serpentinization versus depth for every station, with uncertainty bands."""

    plt.figure(dpi=300, figsize=(45, 25))

    for ii in range(len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        sta = Inversion_container[ii]["Station"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        vs_good = []
        for i in range(len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                vs_good.append(vs[i])

        vs_good = np.array(vs_good)[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)
        vs_downsampled = vs_reshaped.mean(axis=2)

        ax_main = plt.subplot(2, len(Inversion_container), ii + 1)
        plt.title(str("YV.") + sta)
        plt.vlines(90, 0, 100, linestyles="dashed", colors="grey", linewidth=3)
        plt.text(
            95,
            25,
            "Fresh peridotite",
            fontsize=35,
            color="grey",
            ha="center",
            rotation=90,
        )

        depth = np.arange(0, 100, 1)
        plt.plot(
            np.median(vs_downsampled, axis=0) / 50,
            depth,
            color="black",
            linewidth=5,
            linestyle="solid",
            label="Median",
        )
        medians = np.median(vs_downsampled, axis=0) / 50
        std_dev = np.std(vs_downsampled, axis=0) / 50
        plt.fill_betweenx(
            depth, medians - std_dev, medians + std_dev, color="gray", alpha=0.5
        )
        plt.grid(True)
        plt.ylim(100, 0)
        plt.xlabel("Shear Velocity [m/s]")
        plt.ylabel("Depth [m]")

        # Adding inset plot at a new position
        ax_inset = inset_axes(
            ax_main, width="30%", height="50%", loc="lower left", borderpad=1.5
        )
        # Fig. 10c relation (4.2 km/s: 0 %, 2.2 km/s: 100 %); was (9000 - 2 Vs)/30 (F2).
        serpentinization_percentage = serpentinization(
            np.median(vs_downsampled, axis=0)
        )
        ax_inset.plot(
            serpentinization_percentage,
            depth,
            color="blue",
            linewidth=5,
            linestyle="solid",
        )
        ax_inset.fill_betweenx(
            depth,
            serpentinization_percentage - std_dev,
            serpentinization_percentage + std_dev,
            color="blue",
            alpha=0.5,
        )
        ax_inset.set_xlabel("Serpentinization [%]")
        ax_inset.set_ylabel("Depth [m]")
        ax_inset.set_ylim(100, 0)
        ax_inset.grid(True)

    plt.tight_layout()
    plt.show()


# Assuming Inversion_container is predefined somewhere in your code with the necessary data
def plot_inversion_serpentinization_panels(Inversion_container):
    """Variant of ``plot_inversion_serpentinization`` with a different panel layout."""

    plt.figure(dpi=300, figsize=(45, 25))

    for ii in range(0, len(Inversion_container)):
        vs = Inversion_container[ii]["Shear Velocity"]
        mis_fit = Inversion_container[ii]["Misfit Fucntion"]
        sta = Inversion_container[ii]["Station"]
        mis_fit_trsh = Inversion_container[ii]["mis_fit_trsh"]

        bins = 100
        jj = 0
        vs_good = []

        for i in range(0, len(mis_fit[0])):
            if mis_fit[0][i] < mis_fit_trsh:
                # if vs[ii][9000][0] < vs[ii][6000][0]:

                vs_good.append(vs[i])
                jj = jj + 1

        vs_good = np.array(vs_good)

        vs_good = vs_good[:, 0:10000, 0]

        downsample_factor = vs_good.shape[1] // 100

        # Reshape the array to prepare for averaging
        # New shape will be (20000, 100, downsample_factor, 1)
        vs_reshaped = vs_good.reshape(vs_good.shape[0], 100, downsample_factor)

        # Take the mean along the downsample_factor dimension
        vs_downsampled = vs_reshaped.mean(axis=2)

        c = np.zeros([bins, bins])
        b = np.histogram(vs_reshaped[:, 10, 0], bins=bins, range=(0, 5000))[1]

        for i in range(0, len(vs_downsampled[0])):
            c[i] = np.histogram(
                vs_downsampled[:, i], bins=bins, range=(0, 5000), density=True
            )[0]
            c[i] = c[i] / np.max(c[i])

        # Define the colors for the colormap (from white to red to black)

        # Define the transition points for the colors

        # Create the custom colormap
        plt.subplot(2, len(Inversion_container), int(ii + 1))
        plt.title(str("YV.") + sta)
        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        plt.vlines(90, 0, 100, linestyles="dashed", colors="grey", linewidth=3)
        plt.text(
            95,
            25,
            "Fresh peridotite",
            fontsize=35,
            color="grey",
            ha="center",
            rotation=90,
        )

        plt.plot(
            np.median(vs_downsampled, axis=0) / 50,
            depth,
            color="black",
            linewidth=5,
            linestyle="solid",
            label="Median",
        )
        medians = np.median(vs_downsampled, axis=0) / 50
        std_dev = np.std(vs_downsampled, axis=0) / 50

        plt.fill_betweenx(
            depth, medians - std_dev, medians + std_dev, color="gray", alpha=0.5
        )
        plt.grid(True)

        # Setting custom y-axis ticks and labels
        plt.xticks(ticks=selected_indices, labels=selected_labels)

        selected_indices_y, selected_labels_y = _legacy_depth_ticks()
        plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
        if ii == 0 or ii == 4:
            plt.ylabel("Depth [m]")
        if ii in range(len(Inversion_container)):
            plt.xlabel("Shear Velocity [m/s]")
        # if ii == 3 or ii == 7:
        # if ii == 7:
        plt.legend(loc="lower left", fontsize=30)
        plt.ylim(_legacy_rows(8000), 0)

        plt.subplot(2, len(Inversion_container), int(ii + len(Inversion_container) + 1))
        selected_indices = np.linspace(0, len(b) - 1, 6, dtype=int)
        selected_labels = b[selected_indices]
        selected_labels = selected_labels.astype(int)

        depth = np.arange(0, 100, 1)

        # Fig. 10c relation (4.2 km/s: 0 %, 2.2 km/s: 100 %); was (9000 - 2 Vs)/30 (F2).
        serpentinization_percentage = serpentinization(
            np.median(vs_downsampled, axis=0)
        )
        plt.plot(
            serpentinization_percentage,
            depth,
            color="blue",
            linewidth=5,
            linestyle="solid",
        )

        medians = np.median(vs_downsampled, axis=0) / 50
        std_dev = np.std(vs_downsampled, axis=0) / 50

        plt.fill_betweenx(
            depth,
            serpentinization_percentage - std_dev,
            serpentinization_percentage + std_dev,
            color="blue",
            alpha=0.5,
        )

        plt.grid(True)
        plt.xlim([0, 100])
        selected_indices_y, selected_labels_y = _legacy_depth_ticks()
        plt.yticks(ticks=selected_indices_y, labels=selected_labels_y)
        if ii == 0 or ii == len(Inversion_container):
            plt.ylabel("Depth [m]")
        if ii in range(len(Inversion_container)):
            plt.xlabel("Serpentinization [%]")
        # if ii == 3 or ii == 7:
        # if ii == 7:
        plt.ylim(_legacy_rows(8000), 0)

    plt.tight_layout()


def plot_hist(starting_model, burnin, mis_fit, mis_fit_trsh):
    """Histograms of each layer's Vs and thickness for models below the misfit threshold after burn-in."""

    starting_model_opt = []
    for i in range(burnin, len(starting_model[0][0][:]) - 1):
        if mis_fit[0][i] < mis_fit_trsh:
            starting_model_opt.append(starting_model[:, :, i])
    starting_model_opt = np.array(starting_model_opt)
    plt.rcParams.update({"font.size": 40})

    plt.figure(dpi=300, figsize=(35, 70))
    plt.suptitle("Shear Velocity Distribution of Layers", y=0.99)
    for i in range(0, (len(starting_model_opt[0]) - 2)):
        plt.subplot((len(starting_model_opt[0]) - 1), 4, i + 1)

        plt.hist(
            starting_model_opt[:, i, 3],
            bins=20,
            density=True,
            histtype="barstacked",
            facecolor="k",
            alpha=1,
            stacked=True,
            label="Layer " + str(i + 1),
        )
        plt.xlim([0, 4500])
        plt.legend(loc="upper left")
        plt.tight_layout()

    plt.figure(dpi=300, figsize=(35, 70))
    plt.suptitle("Thickness Distribution of Layers", y=0.99)
    for i in range(0, (len(starting_model_opt[0]) - 2)):
        plt.subplot((len(starting_model_opt[0]) - 1), 4, i + 1)

        plt.hist(
            starting_model_opt[:, i, 0],
            bins=20,
            histtype="barstacked",
            facecolor="k",
            alpha=1,
            stacked=True,
            label="Layer " + str(i + 1),
        )
        plt.xlim([0, 2000])
        plt.legend(loc="upper left")
        plt.tight_layout()


def plot_hist2d(starting_model, burnin, mis_fit, mis_fit_trsh):
    """2-D Vs-thickness histograms per layer for models below the misfit threshold after burn-in."""

    starting_model_opt = []
    for i in range(burnin, len(starting_model[0][0][:]) - 1):
        if mis_fit[0][i] < mis_fit_trsh:
            starting_model_opt.append(starting_model[:, :, i])
    starting_model_opt = np.array(starting_model_opt)

    plt.rcParams.update({"font.size": 40})
    plt.figure(dpi=300, figsize=(50, 100))
    plt.suptitle("Shear Velocity and Thickness Distribution of Layers", y=0.99)
    for i in range(0, (len(starting_model_opt[0]) - 2)):
        plt.subplot((len(starting_model_opt[0]) - 1), 4, i + 1)

        plt.hist2d(
            starting_model_opt[:, i, 0],
            starting_model_opt[:, i, 3],
            density=True,
            bins=10,
            cmap="jet",
        )
        # bins = 50, density = True, histtype = "barstacked",facecolor='k', alpha=1,stacked = True,label='Layer ' + str(i+1))
        plt.xlabel("Thickness[m]")
        plt.ylabel("Vs [m/s]")
        plt.colorbar()
        plt.tight_layout()


def plot_autocorrelation(starting_model, N):
    """Autocorrelation of each layer's thickness over ``N`` lags and the effective sample size (Gelman et al. 2013)."""

    cc = np.zeros([len(starting_model), len(starting_model[0, 3, :])])

    for i in range(0, len(cc) - 1):
        s = starting_model[i, 3, :]
        dd = np.correlate(s - np.mean(s), s - np.mean(s), "full") / np.sum(
            (s - np.mean(s)) ** 2
        )
        # Auto-correlations.
        cc[i] = dd[N - 1 :]

    # Estimate of the effective sample size (Gelman et al., 2013).
    Neff = np.zeros([len(starting_model)])
    for i in range(0, len(cc)):
        for j in range(N - 1):
            if cc[i][j] + cc[i][j + 1] > 0.0:
                Neff[i] += cc[i][j]

    Neff = N / (1.0 + 2.0 * Neff)
    # for i in range(0,len(cc)):

    # Plot autocorrelation function.
    # for i in range(0,len(cc)-1):

    plt.figure(dpi=300, figsize=(30, 40))
    plt.rcParams.update({"font.size": 30})
    for i in range(0, len(cc) - 1):
        plt.subplot(len(cc) - 1, 3, i + 1)
        plt.plot(cc[i][0:N], "k", linewidth=3, label="Neff = : %f" % Neff[i])
        plt.xlabel("Iteration", labelpad=15)
        plt.xlim([0, N])
        plt.title("Layer " + str(i + 1))
        plt.legend(loc="upper right", fontsize=15)
        plt.grid()
        plt.tight_layout()

    plt.figure(dpi=300, figsize=(30, 40))
    plt.rcParams.update({"font.size": 30})
    plt.suptitle("Vs")
    for i in range(0, len(cc) - 1):
        plt.subplot(len(cc) - 1, 3, i + 1)
        plt.plot(starting_model[i, 3, :], "k", linewidth=2)
        plt.xlabel("Iteration", labelpad=15)
        plt.xlim([0, N])
        plt.title("Layer " + str(i + 1))
        plt.grid()
        plt.tight_layout()

    plt.figure(dpi=300, figsize=(30, 40))
    plt.rcParams.update({"font.size": 30})
    plt.suptitle("Thickness")
    for i in range(0, len(cc) - 1):
        plt.subplot(len(cc) - 1, 3, i + 1)
        plt.plot(starting_model[i, 0, :], "k", linewidth=2)
        plt.xlabel("Iteration", labelpad=15)
        plt.xlim([0, N])
        plt.title("Layer " + str(i + 1))
        plt.grid()
        plt.tight_layout()

    return Neff


__getattr__ = renamed(
    __name__,
    {
        "plot_inversion_v2": "plot_inversion_summary",
        "plot_inversion_beta": "plot_inversion_with_start",
        "plot_inversion_serpentinization1": "plot_inversion_serpentinization_panels",
    },
)
