#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file name: blazar_plots.py
Program purpose: Code for various plots

Parameters are always listed in the following order:
[delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R]
All parameters are the logarithm of the true value except for delta, n1, and n2
---------------------------------------------------
delta       doppler factor                  linear
K           particle density [cm^-3]        log
n1          n_1 (first index)           linear
n2          n_2 (second index)          linear
gamma_min   low-energy cutoff               log
gamma_max   high-energy cutoff              log
gamma_break energy break                    log
B           magnetic field strength [G]     log
R           blob radius (cm)                log
---------------------------------------------------

Note that data (observed) is expected in linear scale and model data is expected in log scale unless log param is set to False.

"""
import glob
import os
import pickle
import random

import corner
import numpy as np
import tqdm
import subprocess
from matplotlib import pyplot as plt
from scipy import interpolate
from scipy import stats

from bjet_mcmc import blazar_model
from bjet_mcmc.blazar_properties import *
from bjet_mcmc import blazar_utils

__all__ = [
    "cooling_time_Thomson",
    "corner_plot",
    "get_min_max_per_point",
    "get_params_1sigma_ranges",
    "N_e_BknPowLaw",
    "plot_1sigma",
    "plot_1sigma_plots",
    "plot_chain",
    "plot_chi_squared",
    "plot_cooling_times",
    "plot_data",
    "plot_likelihood_profiles",
    "plot_model",
    "plot_model_and_data",
    "plot_particle_spectrum",
    "residual_plot",
    "scale_to_values",
    "sig_T",
    "m_e",
    "c",
]

cmap = plt.get_cmap("tab10")
image_type = "svg"
plt.ioff()


# plots of SED ---------------------------------------------------------------------------------------------------------
def plot_model(
    model,
    title=None,
    no_title=False,
    line=True,
    points=True,
    point_style=".",
    line_style="-",
    line_alpha=1.0,
    file_name=RESULTS_FOLDER + "/model." + image_type,
    clear_plot=True,
    save=False,
    show=True,
    log=True,
):
    """
    Plots frequency against energy flux using matplotlib
    Args:
        results: an index-able type with at least two elements
            The first two elements are 1D numpy arrays of floats. Only the first
            2 elements are used. results[0] should be logv values and results[1]
            logvFv values.
        save (optional): bool
            Whether model should be saved as image; default is False
        show (optional): bool
            Specifies if plot is shown; default is true
        file_name (optional):str
            Where model should be saved; default is "<RESULTS_FOLDER>/model.<image_type>"
        line (optional): bool
            Specifies if the line between points should be shown; default is True
        points (optional): bool
            Specifies if the points should be shown (or just the line); default is True
        point_style (optional): str that is a valid matplotlib point style
            Specifies the point style; default is '.'
        line_style (optional): str that is a valid matplotlib line style
            Specifies the line style; default is '-'
        clear_plot (optional): bool
            Whether the plot should be cleared before new data is plotted; default
            is True

    After function call:
        If show is true, the plot is shown.
        If save is true, the plot is saved as file.
    """
    if log:
        x = np.power(10, model[0])
        y = np.power(10, model[1])
    else:
        x = model[0]
        y = model[1]

    if clear_plot:
        plt.figure("SED Plot")
    plt.xlabel(r"$\nu$")
    plt.ylabel(r"$\nu F_{\nu}$")
    if title is None and not no_title:
        title = "SED frequency/energy flux plot"
    if title is not None:
        plt.title(title)
    if points:
        plt.plot(x, y, point_style)
    if line:
        plt.plot(x, y, line_style, alpha=line_alpha)
    if save:
        plt.savefig(file_name)
    if show:
        plt.show()


def plot_data(
    data_file,
    title=None,
    no_title=False,
    adjust_scale=True,
    lower_adjust_multiplier=None,
    upper_adjust_multiplier=None,
    file_name=RESULTS_FOLDER + "/data." + image_type,
    clear_plot=True,
    save=False,
    show=True,
):
    """
    Create a plot of the data
    Args:
        v_data: 1D np array of floats
            Data
        vFv_data: 1D np array of floats
            Data
        title: str
            Title for plot
        file_name (optional): str
            Relative path to where the image should be saved; default is
            "<RESULTS_FOLDER>>/data.<image_type>"
        save (optional): bool
            If the plot should be saved; default is False
        show (optional): bool
            Specifies if plot is shown; default is True
        adjust_scale (optional): bool
            True = the plot should use the scale of the data
            False = scaled to the model
            default is True
        lower_adjust_multiplier (optional): float
            How far below the data the plot should be scaled; default is 1.1
        upper_adjust_multiplier (optional): float
            How far above the data the plot should be scaled; default is 1.1
        clear_plot (optional): bool
            Whether the plot should be cleared before new data is plotted; default
            is True

    After function call:
        If show is true, the plot is shown.
        If save is true, the plot is saved as file.
    """
    filled_markers = (
        "o",
        "^",
        "<",
        ">",
        "8",
        "s",
        "p",
        "*",
        "h",
        "H",
        "D",
        "d",
        "P",
        "X",
    )
    marker_size = 3

    if clear_plot:
        plt.figure("Plot of data")
    fig, ax = plt.subplots()
    plt.xscale("log")
    plt.xlabel(r"Frequency $\nu$ (Hz)")
    ax.set_ylabel(r"Energy Flux $\nu F_{\nu}$ (erg cm$^{-2}$ s$^{-1})$")
    ax.set_yscale("log")
    if title is None and not no_title:
        title = "Data"
    if title is not None:
        plt.title(title)

    data = blazar_utils.read_data(data_file, instrument=True)
    v_data, vFv_data, err_data, instrument_data, nubin_data = data
    # setting limits
    uplims = [False] * len(v_data)
    lolims = [False] * len(v_data)
    for i in range(len(err_data[1])):
        if err_data[0][i] == 0:
            uplims[i] = True
            err_data[0][i] = vFv_data[i] / 4
        if err_data[0][i] == -1:
            lolims[i] = True
            err_data[1][i] = vFv_data[i] / 4

    list_intruments = [instrument_data[0]]
    v_data_inst = [v_data[0]]
    vFv_data_inst = [vFv_data[0]]
    err_data_inst_down = [err_data[0][0]]
    err_data_inst_up = [err_data[1][0]]
    nubin_data_inst_low = [nubin_data[0][0]]
    nubin_data_inst_high = [nubin_data[1][0]]
    uplims_inst = [uplims[0]]
    lolims_inst = [lolims[0]]
    tmp = 0
    tmp2 = 0

    for i in range(1, len(instrument_data)):
        # cycle colors & markers
        if len(list_intruments) - tmp >= cmap.N:
            tmp += cmap.N
        color_index = len(list_intruments) - tmp
        if len(list_intruments) - tmp2 >= len(filled_markers):
            tmp2 += len(filled_markers)
        marker_index = len(list_intruments) - tmp2

        if instrument_data[i] != list_intruments[-1]:
            ax.errorbar(
                v_data_inst,
                vFv_data_inst,
                xerr=(nubin_data_inst_low, nubin_data_inst_high),
                yerr=(err_data_inst_down, err_data_inst_up),
                uplims=uplims_inst,
                lolims=lolims_inst,
                fmt=filled_markers[marker_index - 1],
                label=str(list_intruments[-1]),
                elinewidth=1,
                markersize=marker_size,
                color=cmap(color_index),
            )
            list_intruments.append(instrument_data[i])
            v_data_inst = [v_data[i]]
            vFv_data_inst = [vFv_data[i]]
            err_data_inst_down = [err_data[0][i]]
            err_data_inst_up = [err_data[1][i]]
            nubin_data_inst_low = [nubin_data[0][i]]
            nubin_data_inst_high = [nubin_data[1][i]]
            uplims_inst = [uplims[i]]
            lolims_inst = [lolims[i]]
        else:
            v_data_inst.append(v_data[i])
            vFv_data_inst.append(vFv_data[i])
            err_data_inst_down.append(err_data[0][i])
            err_data_inst_up.append(err_data[1][i])
            nubin_data_inst_low.append(nubin_data[0][i])
            nubin_data_inst_high.append(nubin_data[1][i])
            uplims_inst.append(uplims[i])
            lolims_inst.append(lolims[i])
        if i == len(instrument_data) - 1:
            ax.errorbar(
                v_data_inst,
                vFv_data_inst,
                xerr=(nubin_data_inst_low, nubin_data_inst_high),
                yerr=(err_data_inst_down, err_data_inst_up),
                uplims=uplims_inst,
                lolims=lolims_inst,
                fmt=filled_markers[marker_index - 1],
                label=str(list_intruments[-1]),
                elinewidth=1,
                markersize=marker_size,
                color=cmap(color_index),
            )
    ax.legend(loc="best", ncol=2)

    secax = ax.secondary_xaxis(
        "top", functions=(blazar_utils.v_to_e, blazar_utils.e_to_v)
    )
    secax.set_xlabel("Energy (eV)")

    if adjust_scale:
        new_min, new_max = scale_to_values(
            vFv_data,
            lower_adjust_multiplier=lower_adjust_multiplier,
            upper_adjust_multiplier=upper_adjust_multiplier,
        )
        plt.ylim(new_min, new_max)

    if save:
        plt.savefig(BASE_PATH + file_name)
    if show:
        plt.show()


def plot_model_and_data(
    model,
    data_file,
    flat_samples,
    indices_within_1sigma,
    redshift,
    eic,
    title=None,
    no_title=False,
    adjust_scale=True,
    lower_adjust_multiplier=None,
    upper_adjust_multiplier=None,
    file_name=RESULTS_FOLDER + "/model_and_data." + image_type,
    save=False,
    show=True,
    log=True,
    fixed_params=None,
    verbose=False,
):
    """
    Create a plot with data and model data
    Args:
        model: tuple of 1D np arrays of floats
            The data for the model; only the first 2 elems are used, which are
            logv and logvFv
        title: str
            Title for plot
        file_name (optional): str
            Relative path to where the image should be saved; default is
            "<RESULTS_FOLDER>/model_and_data.<image_type>"
        save (optional): bool
            If the plot should be saved; default is False
        show (optional): bool
            Specifies if plot is shown; default is True
        adjust_scale (optional): bool
            True = the plot should use the scale of the data
            False = scaled to the model
            default is True
        lower_adjust_multiplier (optional): float
            How far below the data the plot should be scaled; default is 1.1
        upper_adjust_multiplier (optional): float
            How far above the data the plot should be scaled; default is 1.1

    After function call:
        If show is true, the plot is shown.
        If save is true, the plot is saved as file.
    """

    params = {  #'backend': 'ps',
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
    }
    plt.rcParams.update(params)

    # if clear_plot:
    #     plt.figure("Model and Data")

    plot_data(
        data_file,
        title=title,
        no_title=no_title,
        adjust_scale=adjust_scale,
        lower_adjust_multiplier=lower_adjust_multiplier,
        upper_adjust_multiplier=upper_adjust_multiplier,
        clear_plot=False,
        save=False,
        show=False,
    )

    # I would rather not having this to read the configs
    # what if not standard config. It should only read from the output
    # configs = blazar_utils.read_configs()
    # eic = configs["eic"]
    # redshift = configs["redshift"]

    name_stem = "plot"
    min_per_param, max_per_param = get_params_1sigma_ranges(
        flat_samples, indices_within_1sigma, eic=eic, fixed_params=fixed_params
    )
    to_plot = np.concatenate((np.array(min_per_param), np.array(max_per_param)))

    command_params_1, command_params_2 = blazar_model.command_line_sub_strings(
        name_stem=name_stem, redshift=redshift, prev_files=False, eic=eic
    )
    command_params_2[3] = "99"  # number of points for the SED
    lowest_per_point, highest_per_point = get_min_max_per_point(
        model[0],
        to_plot,
        name_stem=name_stem,
        redshift=redshift,
        command_params_1=command_params_1,
        command_params_2=command_params_2,
        eic=eic,
        fixed_params=fixed_params,
    )
    x, y = (model[2], model[3])
    plt.plot(x, y, label="Best model")
    plt.fill_between(
        x,
        np.power(10, lowest_per_point),
        np.power(10, highest_per_point),
        alpha=0.5,
        label=r"Within 1$\sigma$",
    )
    plt.legend(loc="best", ncol=3)
    plt.xlim(1e8, 1e28)
    plt.tight_layout()
    if save:
        plt.savefig(BASE_PATH + file_name)
    if show:
        plt.show()


# MCMC Plots -----------------------------------------------------------------------------------------------------------
def corner_plot(
    values,
    param_min_vals,
    param_max_vals,
    best_params,
    sigma_below_params,
    sigma_above_params,
    title=None,
    no_title=False,
    param_names=None,
    file_name=RESULTS_FOLDER + "/corner." + image_type,
    save=False,
    show=True,
    dpi=300,
    eic=False,
    fixed_params=None,
):
    """
    Create a corner plot:
    Args:
        values: 2D np array of arrays w/ NUM_DIM columns
            Flat samples--a list of sets of parameters
        param_min_vals: list or np array of NUM_DIM floats
            Minimum value for each param
        param_max_vals: list or np array of NUM_DIM floats
            Maximum value for each param
        best_params: list or np array of NUM_DIM floats
            Parameter values for the model with the best chi squared value
        sigma_below_params: list or np array of NUM_DIM floats
            Value for cutoff of lowest param value within 1 sigma for each param
        sigma_above_params: list or np array of NUM_DIM floats
            Value for cutoff of highest param value within 1 sigma for each param
        param_names (optional): list of strings
            Names of the parameters (param_names for the plot). These should be
            formatted params for math text. Default is None; then they will be
            set to FORMATTED_PARAM_NAMES.
        save (optional): bool
            If the plot should be saved; default is False
        show (optional): bool
            Specifies if plot is shown; default is True
        file_name (optional): str
            Relative path to where the image should be saved; default is
            "corner.<image_type>"
    After function call:
        mpl plot is shown, if save is true, plot is saved to a file
    """
    if param_names is None:
        param_names = modelProperties(
            eic, fixed_params=fixed_params
        ).FORMATTED_PARAM_NAMES

    min_maxes = []
    for i in range(modelProperties(eic, fixed_params=fixed_params).NUM_DIM):
        min_maxes.append([param_min_vals[i], param_max_vals[i]])
    fig = corner.corner(
        values,
        labels=param_names,
        range=min_maxes,
        use_math_text=True,
        plot_datapoints=False,
        fill_contours=True,
        label_kwargs={"fontsize": 15},
        labelpad=0.28,
        max_n_ticks=4,
    )

    fig.subplots_adjust(top=0.97, bottom=0.08, left=0.07, right=0.88)

    # extract axes
    dims = modelProperties(eic, fixed_params=fixed_params).NUM_DIM
    axes = np.array(fig.axes).reshape((dims, dims))
    # loop over the diagonal, add lines on histogram
    for i in range(dims):
        ax = axes[i, i]
        ax.tick_params(axis="both", labelsize=15)
        ax.axvline(best_params[i], color="r")
        ax.axvline(sigma_below_params[i], color="b", ls="--")
        ax.axvline(sigma_above_params[i], color="b", ls="--")
        best = "{:.3f}".format(best_params[i])
        below = "{:.3f}".format(sigma_below_params[i] - best_params[i])
        above = "{:.3f}".format(sigma_above_params[i] - best_params[i])
        text = param_names[i] + r" = $" + best + "^{" + above + "}_{" + below + "}$"
        ax.set_title(text, loc="left", fontsize=15)

    # add crossing lines for best point
    for i in range(1, dims):
        for j in range(i):
            ax = axes[i, j]
            ax.tick_params(axis="both", labelsize=15)
            ax.axhline(best_params[i], color="red", ls="dotted")
            ax.axvline(best_params[j], color="red", ls="dotted")
            ax.plot(best_params[j], best_params[i], ".", color="red", markersize=2)
    if not eic:
        fig.set_size_inches(11, 10)
    else:
        fig.set_size_inches(17, 17)
    if save:
        plt.savefig(BASE_PATH + file_name, dpi=dpi)
    if show:
        plt.show()
    return fig


def plot_chain(
    chain,
    param_names=None,
    file_name="chain." + image_type,
    save=False,
    show=True,
    eic=False,
):
    """
    Args:
        chain (numpy.ndarray): The chain of samples to plot, usually a 3-dimensional array.
        param_names (list[str], optional): List of parameter names corresponding to the dimensions of the chain. Defaults to None.
        file_name (str, optional): The name of the file to save the plot as. Defaults to "chain.svg".
        save (bool, optional): Whether to save the plot as an image file. Defaults to False.
        show (bool, optional): Whether to display the plot. Defaults to True.
        eic (bool, optional): Whether to use EIC (electrochemical impedance spectroscopy) properties. Defaults to False.

    """
    if param_names is None:
        param_names = modelProperties(eic).FORMATTED_PARAM_NAMES
    fig, axes = plt.subplots(modelProperties(eic).NUM_DIM, sharex="all")
    for i in range(modelProperties(eic).NUM_DIM):
        ax = axes[i]
        ax.plot(chain[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(chain))
        ax.set_ylabel(param_names[i])
    if save:
        plt.savefig(BASE_PATH + file_name)
    if show:
        plt.show()


def plot_chi_squared(
    values,
    discard_number,
    use_log_probs=True,
    plot_type="med",
    title=None,
    no_title=False,
    fmt="",
    file_name="med_chi_squared_plot." + image_type,
    save=False,
    show=True,
    clear_plot=True,
):
    """
    Args:
        values (np.ndarray): The values of the chi-squared statistic.
        discard_number (int): The number of steps to discard from the values.
        use_log_probs (bool): Whether to use logarithmic values for the chi-squared statistic. Defaults to True.
        plot_type (str): The type of plot to generate. Must be one of "med", "best", or "all". Defaults to "med".
        title (str): The title of the plot. If not provided, a default title will be used depending on the plot_type.
        no_title (bool): Whether to show a title on the plot. Defaults to False.
        fmt (str): The format string for the plot lines. Defaults to an empty string.
        file_name (str): The name of the file to save the plot as. Defaults to "med_chi_squared_plot.svg", where image_type should be replaced with the desired image type.
        save (bool): Whether to save the plot as an image file. Defaults to False.
        show (bool): Whether to display the plot. Defaults to True.
        clear_plot (bool): Whether to clear the plot before generating the new plot. Defaults to True.
    """
    if plot_type not in ["med", "best", "all"]:
        raise ValueError(plot_type + " is not a valid plot type")

    if clear_plot:
        plt.figure("Plot_of_chi^2")
    if use_log_probs:
        chi_sq = -2.0 * values
    else:
        chi_sq = 1.0 * values

    if plot_type == "med":
        chi_sq = np.median(chi_sq, axis=1)
        if title is None and not no_title:
            title = r"Median $\chi^2$ by Step"
        plt.plot(
            np.arange(discard_number, discard_number + len(chi_sq)), chi_sq[:], fmt
        )
    elif plot_type == "best":
        chi_sq = np.min(chi_sq, axis=1)
        if title is None and not no_title:
            title = r"Best $\chi^2$ by Step"
        plt.plot(
            np.arange(discard_number, discard_number + len(chi_sq)), chi_sq[:], fmt
        )
    else:
        if title is None and not no_title:
            title = r"$\chi^2$ by Step"
    if title is not None:
        plt.title(title)
    if plot_type == "all":
        plt.plot(
            np.arange(discard_number, discard_number + len(chi_sq)),
            chi_sq[:, :],
            fmt,
            alpha=0.3,
        )
        plt.semilogy()

    plt.xlim(discard_number, discard_number + len(chi_sq))

    plt.xlabel("step")
    plt.ylabel(r"$\chi^2$")

    if save:
        plt.savefig(FOLDER_PATH + file_name)
    if show:
        plt.show()
    if clear_plot and show == False:
        plt.close("Plot_of_chi^2")


def plot_1sigma(
    v_data,
    vFv_data,
    err_data,
    indices_within_1sigma,
    flat_samples,
    min_chi_squared_index,
    both=False,
    extreme=True,
    title=None,
    no_title=False,
    folder=None,
    file=None,
    save=False,
    show=True,
    serialize=False,
    lower_adjust_multiplier=None,
    upper_adjust_multiplier=1.02,
    max_num_lines_to_graph=1000,
    dims=None,
    eic=False,
    theta=None,
    redshift=None,
    min_freq=None,
    max_freq=None,
    executable=None,
    data_folder=None,
    name_stem=None,
    command_params_full=None,
    command_params_1=None,
    command_params_2=None,
    torus_temp=None,
    torus_luminosity=None,
    torus_frac=None,
    verbose=False,
):
    """
    Args:
        v_data : The observed values of the independent variable.
        vFv_data : The observed values of the dependent variable.
        err_data : The error values associated with the observed dependent variable.
        indices_within_1sigma : The indices of the models within 1 sigma range.
        flat_samples : The flat samples obtained from the MCMC algorithm.
        min_chi_squared_index : The index of the model with minimum chi-squared value.
        both : Boolean value indicating whether both extreme and random models should be plotted. Default is False.
        extreme : Boolean value indicating whether extreme models should be plotted. Default is True.
        title : The title of the plot. If None, a default title will be used.
        no_title : Boolean value indicating whether to hide the title. Default is False.
        folder : The folder path where the plot will be saved. If None, the plot will not be saved.
        file : The file name to save the plot. If None, a default file name will be used.
        save : Boolean value indicating whether to save the plot. Default is False.
        show : Boolean value indicating whether to display the plot. Default is True.
        serialize : Boolean value indicating whether to serialize the plot using pickle. Default is False.
        lower_adjust_multiplier : The lower adjustment multiplier for scaling the y-axis. If None, the minimum y-value of the observed dependent variable will be used.
        upper_adjust_multiplier : The upper adjustment multiplier for scaling the y-axis. Default is 1.02.
        max_num_lines_to_graph : The maximum number of lines to graph for random models. Default is 1000.
        dims : The number of dimensions in the model.
        eic : Boolean value indicating whether to include the EIC parameter. Default is False.
        theta : The theta parameter for the model.
        redshift : The redshift parameter for the model.
        min_freq : The minimum frequency parameter for the model.
        max_freq : The maximum frequency parameter for the model.
        executable : The executable path for the model.
        data_folder : The data folder path for the model.
        name_stem : The name stem for the model.
        command_params_full : The command parameters for the full model.
        command_params_1 : The command parameters for the first model.
        command_params_2 : The command parameters for the second model.
        torus_temp : The torus temperature parameter for the model.
        torus_luminosity : The torus luminosity parameter for the model.
        torus_frac : The torus fraction parameter for the model.
        verbose : Boolean value indicating whether to display verbose output during model creation. Default is False.

    Notes:
        TODO what happens when folder is None??
        The parameters within 1 sigma that have the biggest and smallest values for
        each parameter are found, resulting in 2 arrays of dimension NUM_DIM * NUM_DIM.
        Models are created from these, and for each frequency value,
        the minimum and the maximum are found.
        The graph is made by filling in the space between the minimum and maximum
        for each frequency value.
        The best model and the actual data with error bars are plotted on top of this.
    """
    if both:
        descriptor = ""
    elif extreme:
        descriptor = " extreme"
    else:
        descriptor = " random"

    if title is None and not no_title:
        title = (
            "MCMC results with the range from" + descriptor + " models within 1 sigma"
        )

    fig = plt.figure(title)
    ax = fig.add_subplot()

    if name_stem is None:
        delete_after = True
        name_stem = "for_plotting" + str(random.getrandbits(20))
    else:
        delete_after = False
    best_model = blazar_model.make_model(
        flat_samples[min_chi_squared_index],
        name_stem,
        theta=theta,
        redshift=redshift,
        min_freq=min_freq,
        max_freq=max_freq,
        torus_temp=torus_temp,
        torus_luminosity=torus_luminosity,
        torus_frac=torus_frac,
        data_folder=data_folder,
        executable=executable,
        command_params_full=command_params_full,
        command_params_1=command_params_1,
        command_params_2=command_params_2,
        prev_files=False,
        use_param_file=False,
        verbose=verbose,
        eic=eic,
    )

    # going to get the values to fill between
    logv = best_model[0].copy()
    if dims is None:
        dims = modelProperties(eic).NUM_DIM
    to_plot = np.empty((0, dims))

    if extreme or both:
        min_per_param, max_per_param = get_params_1sigma_ranges(
            flat_samples, indices_within_1sigma, eic=eic
        )
        to_plot = np.concatenate((np.array(min_per_param), np.array(max_per_param)))
    if not extreme or both:
        samples = np.unique(flat_samples[indices_within_1sigma], axis=0)
        num_lines_to_graph = min(len(samples), max_num_lines_to_graph)
        if num_lines_to_graph != len(samples):
            indices = np.random.choice(
                np.arange(0, len(samples)), size=num_lines_to_graph, replace=False
            )
        else:
            indices = np.arange(0, len(samples))
        to_plot = np.concatenate(
            (to_plot, samples[indices])
        )  # if both, both sets will be plotted

    lowest_per_point, highest_per_point = get_min_max_per_point(
        logv,
        to_plot,
        name_stem=name_stem,
        theta=theta,
        redshift=redshift,
        min_freq=min_freq,
        max_freq=max_freq,
        torus_temp=torus_temp,
        torus_luminosity=torus_luminosity,
        torus_frac=torus_frac,
        data_folder=data_folder,
        executable=executable,
        command_params_full=command_params_full,
        command_params_1=command_params_1,
        command_params_2=command_params_2,
        verbose=False,
        eic=eic,
    )

    ax.fill_between(
        logv, lowest_per_point, highest_per_point, alpha=0.5, label=r"Within 1 sigma"
    )

    plt.plot(
        best_model[0], best_model[1], "-", color="g", linewidth=2, label="Best model"
    )

    plt.errorbar(
        v_data, vFv_data, yerr=(err_data), fmt=".", color="b", ecolor="k", label="Data"
    )

    plt.xlabel(r"log $\nu$")
    plt.ylabel(r"log $\nu F_{\nu}$")
    new_min, new_max = scale_to_values(
        vFv_data,
        lower_adjust_multiplier=lower_adjust_multiplier,
        upper_adjust_multiplier=upper_adjust_multiplier,
    )
    plt.ylim(new_min, new_max)
    plt.legend(loc="best")
    if title is not None:
        plt.title(title)

    if delete_after:
        for f in glob.glob(BASE_PATH + DATA_FOLDER + "/" + name_stem + "_*"):
            os.remove(f)

    if save or serialize:
        if folder is None:
            folder = ""
        if folder[-1] != "/":
            folder = folder + "/"
        if descriptor != "" and descriptor is not None:
            file_name = BASE_PATH + folder + "plot_with_" + descriptor[1:] + "_params."
        else:
            file_name = BASE_PATH + folder + "plot_with_error."
        if save:
            if file is None:
                if descriptor != " extreme":
                    plt.savefig(file_name + "png")  # too big for svg
                else:
                    plt.savefig(file_name + image_type)
            else:
                plt.savefig(BASE_PATH + file)
        if serialize:
            with open(file_name + "pickle", "wb") as f:
                pickle.dump(fig, f)
    if show:
        plt.show()


def plot_1sigma_plots(
    v_data,
    vFv_data,
    err_data,
    indices_within_1sigma,
    flat_samples,
    min_chi_squared_index,
    both=False,
    extreme=True,
    title=None,
    no_title=False,
    folder=None,
    file=None,
    save=False,
    show=True,
    serialize=False,
    lower_adjust_multiplier=None,
    upper_adjust_multiplier=1.02,
    max_num_lines_to_graph=1000,
    dims=None,
    eic=False,
    return_models=False,
    theta=None,
    redshift=None,
    min_freq=None,
    max_freq=None,
    executable=None,
    data_folder=None,
    command_params_full=None,
    command_params_1=None,
    command_params_2=None,
    name_stem=None,
    torus_temp=None,
    torus_luminosity=None,
    torus_frac=None,
    verbose=False,
):
    """
    Args:
        v_data : The observed velocities data.
        vFv_data : The observed flux data.
        err_data : The error in the observed flux data.
        indices_within_1sigma : The indices of models within 1 sigma.
        flat_samples : The samples from the MCMC chains.
        min_chi_squared_index : The index of the minimum chi-squared value in the samples.
        both : A flag indicating whether to plot both extreme and random models within 1 sigma. Default is False.
        extreme : A flag indicating whether to plot extreme models within 1 sigma. Default is True.
        title : The title of the plot. If not provided, a default title will be used.
        no_title : A flag indicating whether to display a title in the plot. Default is False.
        folder : The folder to save the plot. If not provided, the plot will not be saved.
        file : The file name of the saved plot. If not provided, a default name will be used.
        save : A flag indicating whether to save the plot. Default is False.
        show : A flag indicating whether to display the plot. Default is True.
        serialize : A flag indicating whether to serialize the plot. Default is False.
        lower_adjust_multiplier : The lower adjust multiplier for scaling the y-axis. If not provided, a default value will be used.
        upper_adjust_multiplier : The upper adjust multiplier for scaling the y-axis. Default is 1.02.
        max_num_lines_to_graph : The maximum number of lines to graph. Default is 1000.
        dims : The dimensions of the model. If not provided, the dimensions will be determined based on the model properties.
        eic : A flag indicating whether the model is an EIC model. Default is False.
        return_models : A flag indicating whether to return the generated models. Default is False.
        theta : The theta parameter of the model.
        redshift : The redshift parameter of the model.
        min_freq : The minimum frequency parameter of the model.
        max_freq : The maximum frequency parameter of the model.
        executable : The path to the executable file for the model.
        data_folder : The folder containing the data files.
        command_params_full : The full command parameters of the model.
        command_params_1 : The first set of command parameters of the model.
        command_params_2 : The second set of command parameters of the model.
        name_stem : The name stem of the generated models. If not provided, a random name stem will be used.
        torus_temp : The temperature of the torus component of the model.
        torus_luminosity : The luminosity of the torus component of the model.
        torus_frac : The fraction of the torus component of the model.
        verbose : A flag indicating whether to display verbose output. Default is False.

    Returns:
        The generated models if return_models is True, otherwise None.

    Notes:
        TODO what happens when folder is None??
        The parameters within 1 sigma that have the biggest and smallest values for
        each parameter are found, resulting in 2 arrays of dimension NUM_DIM * NUM_DIM.
        Models are created from these, and for each frequency value,
        the minimum and the maximum are found.
        The graph is made by filling in the space between the minimum and maximum
        for each frequency value.
        The best model and the actual data with error bars are plotted on top of this.
    """
    if both:
        descriptor = ""
    elif extreme:
        descriptor = " extreme"
    else:
        descriptor = " random"

    if title is None and not no_title:
        title = "MCMC results with" + descriptor + " models within 1 sigma"

    fig = plt.figure(title)
    ax = fig.add_subplot()

    if name_stem is None:
        delete_after = True
        name_stem = "for_plotting" + str(random.getrandbits(20))
    else:
        delete_after = False
    best_model = blazar_model.make_model(
        flat_samples[min_chi_squared_index],
        name_stem,
        theta=theta,
        redshift=redshift,
        min_freq=min_freq,
        max_freq=max_freq,
        torus_temp=torus_temp,
        torus_luminosity=torus_luminosity,
        torus_frac=torus_frac,
        data_folder=data_folder,
        executable=executable,
        command_params_full=command_params_full,
        command_params_1=command_params_1,
        command_params_2=command_params_2,
        prev_files=False,
        use_param_file=False,
        verbose=verbose,
        eic=eic,
    )

    # going to get the values to fill between
    logv = best_model[0].copy()
    if dims is None:
        dims = modelProperties(eic).NUM_DIM
    to_plot = np.empty((0, dims))

    if extreme or both:
        min_per_param, max_per_param = get_params_1sigma_ranges(
            flat_samples, indices_within_1sigma, eic=eic
        )
        to_plot = np.concatenate((np.array(min_per_param), np.array(max_per_param)))
    if not extreme or both:
        samples = np.unique(flat_samples[indices_within_1sigma], axis=0)
        num_lines_to_graph = min(len(samples), max_num_lines_to_graph)
        if num_lines_to_graph != len(samples):
            indices = np.random.choice(
                np.arange(0, len(samples)), size=num_lines_to_graph, replace=False
            )
        else:
            indices = np.arange(0, len(samples))
        to_plot = np.concatenate(
            (to_plot, samples[indices])
        )  # if both, both sets will be plotted

    name_stem = "plotting" + str(random.getrandbits(20))
    models = []
    for i in tqdm.trange(0, len(to_plot)):
        params = to_plot[i]
        model = blazar_model.make_model(
            params,
            name_stem,
            theta=theta,
            redshift=redshift,
            min_freq=min_freq,
            max_freq=max_freq,
            torus_temp=torus_temp,
            torus_luminosity=torus_luminosity,
            torus_frac=torus_frac,
            data_folder=data_folder,
            executable=executable,
            command_params_full=command_params_full,
            command_params_1=command_params_1,
            command_params_2=command_params_2,
            prev_files=False,
            use_param_file=False,
            verbose=verbose,
            eic=eic,
        )
        models.append(model)
        plot_model(
            model,
            points=False,
            line_style="r-",
            line_alpha=0.5,
            clear_plot=False,
            show=show,
        )
        if return_models:
            models.append(model)

    plt.xscale("log")
    plt.yscale("log")
    plt.plot(
        np.power(10, best_model[0]),
        np.power(10, best_model[1]),
        "-",
        color="g",
        linewidth=2,
        label="Best model",
    )

    plt.errorbar(
        v_data,
        vFv_data,
        yerr=(err_data, err_data),
        fmt=".",
        color="b",
        ecolor="k",
        label="Data",
    )

    plt.xlabel(r"$\nu$")
    plt.ylabel(r"$\nu F_{\nu}$")
    new_min, new_max = scale_to_values(
        vFv_data,
        lower_adjust_multiplier=lower_adjust_multiplier,
        upper_adjust_multiplier=upper_adjust_multiplier,
    )
    plt.ylim(new_min, new_max)
    plt.legend(loc="best")
    if title is not None:
        plt.title(title)

    if delete_after:
        for f in glob.glob(BASE_PATH + DATA_FOLDER + "/" + name_stem + "_*"):
            os.remove(f)

    if folder is not None and folder[-1] != "/":
        folder = folder + "/"
    if folder is None:
        folder = ""
    if descriptor != "" and descriptor is not None:
        file_name = (
            BASE_PATH + folder + "plot_with_lines_" + descriptor[1:] + "_params."
        )
    else:
        file_name = BASE_PATH + folder + "plot_with_error_lines."
    if save:
        if file is None:
            if descriptor != " extreme":
                plt.savefig(file_name + "png")  # too big for svg
            else:
                plt.savefig(file_name + image_type)
        else:
            plt.savefig(BASE_PATH + file)
    if serialize:
        with open(file_name + "pickle", "wb") as f:
            pickle.dump(fig, f)
    if show:
        plt.show()

    return models if return_models else None


# utils ----------------------------------------------------------------------------------------------------------------
def scale_to_values(values, upper_adjust_multiplier=None, lower_adjust_multiplier=None):
    """
    Scales the given values to new minimum and maximum values based on the upper and lower adjust multipliers.

    Args:
        values (array-like): The values to be scaled.
        upper_adjust_multiplier (float, optional): The multiplier to adjust the upper limit of the scaled values. By default, it is set to 5.
        lower_adjust_multiplier (float, optional): The multiplier to adjust the lower limit of the scaled values. By default, it is set to 5.

    Returns:
        tuple: A tuple containing the new minimum and maximum scaled values.
    """
    if upper_adjust_multiplier is None:
        upper_adjust_multiplier = 5
    if lower_adjust_multiplier is None:
        lower_adjust_multiplier = 5
    data_min = np.min(values)
    data_max = np.max(values)
    if data_min < 0:
        new_min = lower_adjust_multiplier * data_min
    else:
        new_min = 1 / lower_adjust_multiplier * data_min
    if data_max < 0:
        new_max = 1 / upper_adjust_multiplier * data_max
    else:
        new_max = upper_adjust_multiplier * data_max
    return new_min, new_max


def get_params_1sigma_ranges(
    flat_samples, indices_within_1sigma, eic=False, fixed_params=None
):
    """
    Finds the array of parameters that has the minimum and maximum value for
    each of the parameters.
    For example, with samples [[0, 1, 2], [5, 2, 1], [4, 6, 1], [3, 2, 0]], the
    minima would be [[0, 1, 2], [0, 1, 2], [3, 2, 0]] and the maxima would be
    [[5, 2, 1], [4, 6, 1], [0, 1, 2]]

    Arguments:
    flat_samples:
    indices_within_1sigma:

    Returns:
    """
    # set the minima and maxima to the first set of params for all vals
    dims = modelProperties(eic, fixed_params=fixed_params).NUM_DIM
    # print(np.shape(indices_within_1sigma))
    # print(eic)
    minima = [flat_samples[indices_within_1sigma[0]].copy() for _ in range(dims)]
    maxima = [flat_samples[indices_within_1sigma[0]].copy() for _ in range(dims)]
    for index in indices_within_1sigma:
        params = flat_samples[index]
        # print(len(params), len(minima))
        for i in range(dims):
            if params[i] < minima[i][i]:
                minima[i] = params.copy()
            if params[i] > maxima[i][i]:
                maxima[i] = params.copy()
    # print(minima, maxima)
    # print(len(minima))
    return minima, maxima


def get_min_max_per_point(
    v_vals,
    model_params_list,
    name_stem=None,
    theta=None,
    redshift=None,
    min_freq=None,
    max_freq=None,
    torus_temp=None,
    torus_luminosity=None,
    torus_frac=None,
    data_folder=None,
    executable=None,
    command_params_full=None,
    command_params_1=None,
    command_params_2=None,
    verbose=False,
    eic=False,
    fixed_params=None,
):
    """
    Args:
        v_vals (array-like): The array of values at which the model will be computed.
        model_params_list (list): The list of model parameter arrays.
        name_stem (str, optional): The name stem for the output files. If not provided, a random stem will be generated. Defaults to None.
        theta (float, optional): The viewing angle. Defaults to None.
        redshift (float, optional): The redshift of the source. Defaults to None.
        min_freq (float, optional): The minimum frequency for the plot. Defaults to None.
        max_freq (float, optional): The maximum frequency for the plot. Defaults to None.
        torus_temp (float, optional): The temperature of the torus. Defaults to None.
        torus_luminosity (float, optional): The luminosity of the torus. Defaults to None.
        torus_frac (float, optional): The fraction of the torus emission. Defaults to None.
        data_folder (str, optional): The folder where the output files will be saved. Defaults to None.
        executable (str, optional): The path to the executable file for computing the model. Defaults to None.
        command_params_full (str, optional): The command parameters for the full model computation. Defaults to None.
        command_params_1 (str, optional): The command parameters for the first model computation. Defaults to None.
        command_params_2 (str, optional): The command parameters for the second model computation. Defaults to None.
        verbose (bool, optional): Whether to print verbose output. Defaults to False.
        eic (bool, optional): Whether to use effective index of refraction correction. Defaults to False.
        fixed_params (array-like, optional): The array of fixed parameters. Defaults to None.

    Returns:
        tuple: A tuple containing the lowest and highest interpolated values for each point.
    """
    num_points = len(v_vals)
    lowest_per_point = np.array([np.inf for _ in range(num_points)])
    highest_per_point = np.array([-np.inf for _ in range(num_points)])
    print("min_max_per_point")
    if name_stem is None:
        delete_after = True
        name_stem = "for_plotting" + str(random.getrandbits(20))
    else:
        delete_after = False

    for i in tqdm.trange(0, len(model_params_list)):
        params = model_params_list[i]
        # all frozen parameters need to be reimplemented for the model computation
        if fixed_params:
            for i in range(len(fixed_params)):
                if fixed_params[i] != -np.inf:
                    params = np.insert(params, i, fixed_params[i])

        model = blazar_model.make_model(
            params,
            name_stem,
            theta=theta,
            redshift=redshift,
            min_freq=min_freq,
            max_freq=max_freq,
            torus_temp=torus_temp,
            torus_luminosity=torus_luminosity,
            torus_frac=torus_frac,
            data_folder=data_folder,
            executable=executable,
            command_params_full=command_params_full,
            command_params_1=command_params_1,
            command_params_2=command_params_2,
            prev_files=False,
            use_param_file=False,
            verbose=verbose,
            eic=eic,
        )
        f = interpolate.interp1d(model[0], model[1], fill_value="extrapolate")
        vFv_interpolated = f(v_vals)
        for j in range(num_points):
            if vFv_interpolated[j] < lowest_per_point[j]:
                lowest_per_point[j] = vFv_interpolated[j]
            if vFv_interpolated[j] > highest_per_point[j]:
                highest_per_point[j] = vFv_interpolated[j]

    if delete_after:
        for f in glob.glob(BASE_PATH + DATA_FOLDER + "/" + name_stem + "_*"):
            os.remove(f)
    return lowest_per_point, highest_per_point


def residual_plot(data, best_model, lowest_per_point, highest_per_point):
    """
    Args:
        data (tuple): A tuple containing data for the residual plot. It should have the following elements:
            - data[0] (array-like): Frequencies (Hz) for the plot.
            - data[1] (array-like): Residual values for the plot.
            - data[2] (array-like): Errors for the residual values.

        best_model (tuple): A tuple containing data for the best model. It should have the following elements:
            - best_model[2] (array-like): Frequencies (Hz) for the best model.
            - best_model[3] (array-like): Predicted values for the best model.

        lowest_per_point (array-like): Per-point values for the lowest range of the shaded area.

        highest_per_point (array-like): Per-point values for the highest range of the shaded area.
    """
    f_best = interpolate.interp1d(best_model[2], best_model[3])
    f_low = interpolate.interp1d(best_model[2], np.power(10, lowest_per_point))
    f_high = interpolate.interp1d(best_model[2], np.power(10, highest_per_point))
    plt.fill_between(
        data[0],
        f_low(data[0]) - data[1],
        f_high(data[0]) - data[1],
        alpha=0.5,
        zorder=1,
    )
    plt.errorbar(
        data[0], f_best(data[0]) - data[1], yerr=(data[2], data[2]), fmt=".", zorder=2
    )
    plt.axhline(0, zorder=0, color="k", alpha=0.5)
    plt.xlabel(r"Frequency $\nu$ (Hz)")
    plt.ylabel(r"Residual Energy Flux $\nu F_{\nu}$ (erg cm$^{-2}$ s$^{-1})$")
    plt.xscale("log")


def N_e_BknPowLaw(gamma, K, n_1, n_2, gamma_break):
    """
    Args:
        gamma (float): The value of gamma parameter.
        K (int): The value of K parameter.
        n_1 (int): The value of n_1 parameter.
        n_2 (int): The value of n_2 parameter.
        gamma_break (float): The value of gamma_break parameter.

    Returns:
        float: The calculated value of N_e using the broken power law formula.

    """
    K = 10**K
    gamma_break = 10**gamma_break
    K2 = K * gamma_break ** (n_2 - n_1)
    diff = gamma_break - gamma
    sign = diff > 0
    return K * gamma ** (-n_1) * sign + K2 * gamma ** (-n_2) * ~sign


def plot_particle_spectrum(
    best_params,
    min_1sigma_params,
    max_1sigma_params,
    fixed_params,
    file_name=BASE_PATH + RESULTS_FOLDER + "/particle_spectrum.svg",
    save=False,
    show=True,
):
    """
    Plot the broken power-law particle spectrum with 1sigma contour based on Bjet_MCMC outputs
    Note the errors of gamma_min and gamma_max are not included in the contour

    Parameters
    ----------
    best_params : 1D np array of floats
        Array of best parameters found by the MCMC fit.
        For simple SSC without fixed parameters, the order is:
            [delta, K, n_1, n_2, gamma_min, gamma_max, gamma_break, B, R]
        Additional description:
            ---------------------------------------------------
            delta       doppler factor                  linear
            K           particle density [cm^-3]        log
            n1          n_1 (first index)           linear
            n2          n_2 (second index)          linear
            gamma_min   low-energy cutoff               log
            gamma_max   high-energy cutoff              log
            gamma_break energy break                    log
            B           magnetic field strength [G]     log
            R           blob radius [cm]                log
    min_1sigma_params : 1D np array of floats
        Array of 1 sigma lower boundary of free parameters found by the MCMC fit.
        The order should match the one of best_params
    max_1sigma_params : 1D np array
        Array of 1 sigma upper boundary of free parameters found by the MCMC fit.
        The order should match the one of best_params
    fixed_params : list of floats
        List of user-fixed parameters. If no fixed parameters in simple SSC model:
            fixed_params = [-inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf]
    file_name : str, optional
        Absolute path and name for the saved plot.
        The default is BASE_PATH + RESULTS_FOLDER + "/particle_spectrum.svg".
    save : boolean, optional
        True for saving. The default is False.
    show : boolean, optional
        True to display. The default is True.

    Returns
    -------
    The particle spectrum plot.
    """

    # crate a list of particle spectrum parameters (free and fixed)
    params = np.zeros(6)
    min_params = np.zeros(6)
    max_params = np.zeros(6)

    count = 0
    for i in range(7):
        if i == 0 and fixed_params[i] != -np.inf:
            # case of fixed delta
            count = 1
        elif fixed_params[i] != -np.inf:
            params[i - 1] = fixed_params[i]
            min_params[i - 1] = fixed_params[i]
            max_params[i - 1] = fixed_params[i]
            count += 1
        else:
            params[i - 1] = best_params[i - count]
            min_params[i - 1] = min_1sigma_params[i - count]
            max_params[i - 1] = max_1sigma_params[i - count]

    K, n_1, n_2, gamma_min, gamma_max, gamma_break = params
    params_error_down = params - min_params
    params_error_up = max_params - params

    gamma = np.logspace(gamma_min, gamma_max, 500)
    plt.figure("Particle Spectrum")
    plt.plot(
        gamma, N_e_BknPowLaw(gamma, K, n_1, n_2, gamma_break), label="Particle Spectrum"
    )
    plt.loglog()

    if count != 6:
        # count = 6 means that all particle spectrum parameters are fixed, there is no error bars

        # consider asymmetric errors in parameters (two-sided multivariate normal method)
        cov_matrix_down = np.diag(params_error_down) ** 2
        cov_matrix_up = np.diag(params_error_up) ** 2
        rng = np.random.RandomState(seed=30)
        parameter_samples_down = rng.multivariate_normal(params, cov_matrix_down, 5000)
        parameter_samples_up = rng.multivariate_normal(params, cov_matrix_up, 5000)
        parameter_samples = []
        for i in range(len(parameter_samples_down)):
            diff = params - parameter_samples_down[i]
            sign = diff > 0
            # be sure that parameter_samples_up are aways more than the mean
            # by default from our method parameter_samples_down will always be below the mean
            parameter_samples_up[i] = params + np.abs(parameter_samples_up[i] - params)
            # when a projection is below the mean, use parameter_samples_down, when above use parameter_samples_up
            parameter_samples_temp = (
                parameter_samples_down[i] * sign + parameter_samples_up[i] * ~sign
            )
            # remove non-accepted parameter sets from bjet
            if (
                parameter_samples_temp[1] <= parameter_samples_temp[2]
                and parameter_samples_temp[3]
                <= parameter_samples_temp[5]
                <= parameter_samples_temp[4]
            ):
                parameter_samples.append(parameter_samples_temp)
        realizations = np.array(
            [
                N_e_BknPowLaw(gamma, pars[0], pars[1], pars[2], pars[5])
                for pars in parameter_samples
            ]
        )
        q = 100 * stats.norm.cdf(-1)  # 1 is the 1 sigma
        y_low = np.percentile(realizations, q, axis=0)
        q = 100 * stats.norm.cdf(1)  # 1 is the 1 sigma
        y_high = np.percentile(realizations, q, axis=0)

        plt.fill_between(
            gamma, y_low, y_high, edgecolor="None", alpha=0.5, label=r"Within 1$\sigma$"
        )
    plt.xlabel(r"$\gamma$", fontsize=13)
    plt.ylabel(r"Density [cm$^{-3}$]", fontsize=13)
    plt.legend()
    if save:
        plt.savefig(file_name)


sig_T = 6.652453 * 1.0e-25  # [cm^2]
m_e = 9.109558 * 1.0e-28  # [g]
c = 2.997924 * 1.0e10  # [cm / s]


def cooling_time_Thomson(gamma, U_B, U_syn, U_blr):
    """
    Args:
        gamma (float): Lorentz factor of the particle.
        U_B (float): Energy density of the magnetic field.
        U_syn (float): Energy density of the synchrotron radiation.
        U_blr (float): Energy density of the broad-line region radiation.

    Returns:
        float: The cooling time of the particle according to Thomson scattering.

    """
    # see e.g. Inoue & Takahara 1996
    return 3 * m_e * c / (4 * (U_B + U_syn + U_blr) * sig_T * gamma)


def plot_cooling_times(
    logfile,
    best_params,
    fixed_params,
    file_name=RESULTS_FOLDER + "/cooling_times.svg",
    save=False,
    show=True,
    eic=False,
    redshift=None,
):
    """
    Args:
        logfile (str): The name of the log file.
        best_params (tuple): A tuple containing the best parameters.
        fixed_params (tuple): A tuple containing the fixed parameters.
        file_name (str, optional): The name of the file to save the plot. Defaults to RESULTS_FOLDER + "/cooling_times.svg".
        save (bool, optional): Whether to save the plot. Defaults to False.
        show (bool, optional): Whether to show the plot. Defaults to True.
        eic (bool, optional): Whether to include the energy density of the BLR. Defaults to False.
        redshift (float, optional): The redshift of the object. Defaults to None.

    """
    # retrieve the energy densities from bjet.log
    grep = subprocess.run(
        ["grep", "(U_B)", BASE_PATH + logfile], capture_output=True, text=True
    )
    U_B = float(grep.stdout.split()[2])
    grep = subprocess.run(
        ["grep", "(U_syn)", BASE_PATH + logfile], capture_output=True, text=True
    )
    U_syn = float(grep.stdout.split()[2])
    U_blr = 0.0
    if eic:
        grep = subprocess.run(
            ["grep", "(U_blr)", BASE_PATH + logfile], capture_output=True, text=True
        )
        U_blr = float(grep.stdout.split()[2])

    gamma_min, gamma_max = best_params[4:6]
    gamma = np.logspace(gamma_min, gamma_max, 500)
    doppler = best_params[0]

    # plot cooling time in Thomson limit
    plt.figure("Observed cooling time (Thomson limit)")
    plt.plot(
        gamma,
        cooling_time_Thomson(gamma, U_B, U_syn, U_blr) * (1 + redshift) / (doppler),
    )
    plt.loglog()
    plt.xlabel(r"$\gamma$", fontsize=13)
    plt.ylabel(r"$\tau_{\mathrm{cool}}$ [s]", fontsize=13)
    if save:
        plt.savefig(BASE_PATH + file_name)


def plot_likelihood_profiles(
    flat_samples,
    flat_log_probs,
    best_params,
    min_1sigma_params,
    max_1sigma_params,
    save=False,
    show=True,
    fixed_params=None,
    eic=False,
    folder_path=BASE_PATH + RESULTS_FOLDER,
):
    """
    Args:
        flat_samples (numpy.ndarray): The flattened samples of the free parameters from the MCMC run.
        flat_log_probs (numpy.ndarray): The flattened log probabilities corresponding to the flat_samples.
        best_params (numpy.ndarray): The best-fit parameter values.
        min_1sigma_params (numpy.ndarray): The parameter values at the lower 1-sigma confidence level.
        max_1sigma_params (numpy.ndarray): The parameter values at the upper 1-sigma confidence level.
        save (bool, optional): Whether to save the plot. Default is False.
        show (bool, optional): Whether to display the plot. Default is True.
        fixed_params (dict, optional): A dictionary of fixed parameter values. Default is None.
        eic (bool, optional): Whether to include extra inelastic channels. Default is False.
        folder_path (str, optional): The folder path to save the plot. Default is BASE_PATH + RESULTS_FOLDER.

    """
    # plot posterior likelihood profiles for all free parameters.

    param_names = modelProperties(eic, fixed_params=fixed_params).PARAM_NAMES
    fromatted_param_names = modelProperties(
        eic, fixed_params=fixed_params
    ).FORMATTED_PARAM_NAMES

    for j in range(len(best_params)):
        # test to plot likelihood profiles
        param_array = flat_samples[:, j]
        # sort parameter with associated log_prob
        param_array, sorted_log_probs = zip(*sorted(zip(param_array, flat_log_probs)))
        nbins = 40
        Drange = param_array[-1] - param_array[0]
        binsize = Drange / (nbins - 1)
        param_binned = np.linspace(param_array[0], param_array[-1], nbins)
        bin_min = param_array[0]
        bin_max = bin_min + binsize
        prob_tmp = []
        prob_max = []
        ln_prob_max = []
        param_binned_mid = []
        for i in range(len(param_array)):
            if bin_min <= param_array[i] < bin_max:
                prob_tmp.append(sorted_log_probs[i])
            else:
                if len(prob_tmp) == 0:
                    prob_tmp.append(-np.inf)
                prob_max.append(np.exp(max(prob_tmp)))
                ln_prob_max.append(max(prob_tmp))
                param_binned_mid.append(bin_min + binsize / 2.0)
                prob_tmp = []
                bin_min += binsize
                bin_max += binsize

        figure_name = "likelihood_" + param_names[j]
        plt.figure(figure_name)
        plt.plot(param_binned_mid, prob_max, ds="steps-mid", color="0")
        plt.axvline(best_params[j], color="r")
        plt.axvline(min_1sigma_params[j], color="b", ls="--")
        plt.axvline(max_1sigma_params[j], color="b", ls="--")
        plt.xlabel(str(fromatted_param_names[j]), fontsize=13)
        plt.ylabel(r"Likelihood", fontsize=13)

        if save:
            plt.savefig(folder_path + "/" + figure_name + ".svg")
        if show:
            plt.show()