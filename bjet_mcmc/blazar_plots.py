#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program purpose: Code for various plots

Parameters are always listed in the following order:
[delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R]

All parameters are the logarithm of the true value except for delta, n1, and n2

---------------------------------------------------

===================  =============================  ======
Parameter            Description                    Scale
===================  =============================  ======
delta                doppler factor                 linear
K                    particle density [cm^-3]       log
n1                   n_1 (first index)              linear
n2                   n_2 (second index)             linear
gamma_min            low-energy cutoff              log
gamma_max            high-energy cutoff             log
gamma_break          energy break                   log
B                    magnetic field strength [G]    log
R                    blob radius (cm)               log
===================  =============================  ======

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
from scipy.optimize import curve_fit

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
    Plot Model

    :param model: The model data to plot. The first two elements are 1D numpy arrays of floats. Only the first 2 elements are used. results[0] should be logv values and results[1] logvFv values.
    :type model: tuple
    :param title: Optional title for the plot.
    :type title: str
    :param no_title: If True, no title will be displayed on the plot.
    :type no_title: bool
    :param line: If True, a line plot will be created.
    :type line: bool
    :param points: If True, scatter points will be plotted.
    :type points: bool
    :param point_style: Marker style for the scatter points.
    :type point_style: str
    :param line_style: Line style for the line plot.
    :type line_style: str
    :param line_alpha: Alpha value for the line plot.
    :type line_alpha: float
    :param file_name: File name to save the plot.
    :type file_name: str
    :param clear_plot: If True, clears the existing plot before creating a new one.
    :type clear_plot: bool
    :param save: If True, saves the plot to a file.
    :type save: bool
    :param show: If True, displays the plot.
    :type show: bool
    :param log: If True, applies logarithmic scaling to the model data before plotting.
    :type log: bool
    :return: None
    :rtype: None
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
    Plots data from a given file.

    After function call:
        If show is true, the plot is shown.
        If save is true, the plot is saved as file.

    :param data_file: Path to the data file.
    :type data_file: str
    :param title: Title of the plot. Default is None.
    :type title: str, optional
    :param no_title: Whether to include a title in the plot. Default is False.
    :type no_title: bool, optional
    :param adjust_scale: Whether to adjust the y-axis scale of the plot. Default is True.
    :type adjust_scale: bool, optional
    :param lower_adjust_multiplier: Lower multiplier for adjusting the y-axis scale. Default is None.
    :type lower_adjust_multiplier: float, optional
    :param upper_adjust_multiplier: Upper multiplier for adjusting the y-axis scale. Default is None.
    :type upper_adjust_multiplier: float, optional
    :param file_name: Name of the output file to save the plot. Default is "<RESULTS_FOLDER>/<data>.<IMAGE_TYPE>".
    :type file_name: str, optional
    :param clear_plot: Whether to clear the current plot before plotting. Default is True.
    :type clear_plot: bool, optional
    :param save: Whether to save the plot as an image file. Default is False.
    :type save: bool, optional
    :param show: Whether to display the plot. Default is True.
    :type show: bool, optional
    :return: None
    :rtype: None
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

    After function call:
        If show is true, the plot is shown.
        If save is true, the plot is saved as file.

    :param model: The data for the model; only the first 2 elems are used, which are logv and logvFv
    :type model: Any
    :param data_file: The data file to plot.
    :type data_file: str
    :param flat_samples: The flat samples.
    :type flat_samples: Any
    :param indices_within_1sigma: The indices within 1 sigma.
    :type indices_within_1sigma: Any
    :param redshift: The redshift value.
    :type redshift: Any
    :param eic: The eic value.
    :type eic: Any
    :param title: The title of the plot. (optional)
    :type title: str, default=None
    :param no_title: Whether to exclude the title from the plot. (optional)
    :type no_title: bool, default=False
    :param adjust_scale: Whether to adjust the scale of the plot. True = the plot should use the scale of the data. False = scaled to the model. (optional)
    :type adjust_scale: bool, default=True
    :param lower_adjust_multiplier: How far below the data the plot should be scaled. (optional)
    :type lower_adjust_multiplier: Any, default=1.1
    :param upper_adjust_multiplier: How far above the data the plot should be scaled. (optional)
    :type upper_adjust_multiplier: Any, default=1.1
    :param file_name:  Relative path to where the image should be saved. (optional)
    :type file_name: str, default=RESULTS_FOLDER + "/model_and_data." + image_type
    :param save: Whether to save the plot to a file. (optional)
    :type save: bool, default=False
    :param show: Whether to show the plot. (optional)
    :type show: bool, default=True
    :param log: Whether to use a logarithmic scale on the plot. (optional)
    :type log: bool, default=True
    :param fixed_params: The fixed parameters. (optional)
    :type fixed_params: Any, default=None
    :param verbose: Whether to display verbose output. (optional)
    :type verbose: bool, default=False
    :return: None
    :rtype: None
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
    This function generates a corner plot for given values and parameters.

    After function call:
        mpl plot is shown, if save is true, plot is saved to a file

    :param values: The values to plot on the corner plot. 2D np array of arrays w/ NUM_DIM columns. Flat samples--a list of sets of parameters
    :type values: numpy.ndarray

    :param param_min_vals: The minimum values for each parameter. Array of NUM_DIM floats.
    :type param_min_vals: list

    :param param_max_vals: The maximum values for each parameter. Array of NUM_DIM floats.
    :type param_max_vals: list

    :param best_params: Parameter values for the model with the best chi squared value. Array of NUM_DIM floats.
    :type best_params: list

    :param sigma_below_params: Value for cutoff of lowest param value within 1 sigma for each param.Array of NUM_DIM floats.
    :type sigma_below_params: list

    :param sigma_above_params: Value for cutoff of highest param value within 1 sigma for each param. Array of NUM_DIM floats.
    :type sigma_above_params: list

    :param title: The title of the corner plot. (optional)
    :type title: str

    :param no_title: If True, the corner plot will not have a title. (optional)
    :type no_title: bool

    :param param_names: Names of the parameters (param_names for the plot). These should be formatted params for math text. Default is None; then they will be set to FORMATTED_PARAM_NAMES.
    :type param_names: list

    :param file_name: Relative path to where the image should be saved; default is
            "corner.<image_type>"
    :type file_name: str

    :param save: If True, the plot will be saved. (optional, default: False)
    :type save: bool

    :param show: If True, the plot will be displayed. (optional, default: True)
    :type show: bool

    :param dpi: The DPI (dots per inch) for saving the plot. (optional, default: 300)
    :type dpi: int

    :param eic: If True, the plot will be generated for EIC data. (optional, default: False)
    :type eic: bool

    :param fixed_params: The fixed parameters. (optional)
    :type fixed_params: dict

    :return: The generated corner plot figure.
    :rtype: matplotlib.figure.Figure
    """
    if param_names is None:
        param_names = modelProperties(
            eic, fixed_params=fixed_params
        ).FORMATTED_PARAM_NAMES

    min_maxes = []
    for i in range(modelProperties(eic, fixed_params=fixed_params).NUM_DIM):
        #automatically rescale x axis if the parameter space is too narrow
        custom_low = best_params[i]-5*(best_params[i]-sigma_below_params[i])
        custom_high = best_params[i]+5*(sigma_above_params[i]-best_params[i])
        if param_min_vals[i] < custom_low:
            min_val = custom_low
        else:
            min_val = param_min_vals[i]
        if param_max_vals[i] > custom_high:
            max_val = custom_high
        else:
            max_val = param_max_vals[i]
        min_maxes.append([custom_low, custom_high])
        
    
        
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
    
    if len(param_names) > 6:
        fig.subplots_adjust(top=0.97, bottom=0.08, left=0.07, right=0.88)
    else:
        fig.subplots_adjust(top=0.96, bottom=0.11, left=0.1, right=0.9)

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
    Plot the given chain data.

    :param chain: The chain of samples to plot, usually a 3-dimensional array.
    :type chain: numpy.ndarray
    :param param_names: List of parameter names corresponding to the dimensions of the chain. Default is None.
    :type param_names: list, optional
    :param file_name: The name of the file to save the plot as. Default is "chain.svg".
    :type file_name: str, optional
    :param save: Whether to save the plot to a file. Default is False.
    :type save: bool, optional
    :param show: Whether to display the plot. Default is True.
    :type show: bool, optional
    :param eic: Whether the chain contains EIC data. Default is False.
    :type eic: bool, optional
    :return: None
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

# x values for the data
# y values for the data
# params to pass the fitting function
# names to express the parameters
# fit_type: optional string to use one of the built in fitting functions
# fit_func: optional callable to pass a custom fitting function
def create_curve_fit(x_data, y_data, params, param_names=[], 
                     fit_type="", fit_func=None, fit_func_tex=None):
    func = None
    func_tex = None

    # fitted y_vals, params, params_string
    ret = []

    # fitting functions should take, x, and then a list of parameters
    def exponential_decay(x, *p):
        return p[0] * np.exp(-x/p[1]) + p[2]

    if fit_func is not None:
        func = fit_func
    elif fit_type == "" or fit_type=="exponential_decay":
        func = exponential_decay
    else:
        ret = [np.zeros(len(x_data)), np.zeros(len(params)), ""]
        return ret

    if fit_func_tex is not None:
        func_tex = fit_func_tex
    elif fit_type=="exponential_decay":
        func_tex = r"$N_0e^{-x/\tau}+c$"
    else:
        func_tex = ""



    popt, pcov = curve_fit(func, x_data, y_data, p0=params)

    fitted_curve = exponential_decay(x_data, *popt)

    params_str=""
    for i in range(len(popt)):
        pname = ""
        if (i < len(param_names)):
            pname = param_names[i]
        else:
            pname = f"P{i}"
        params_str+=f"{pname}: {popt[i]:.3e}\n"

    ret = [fitted_curve, popt, params_str, func_tex]
    return ret




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
    This function plots the chi-squared values against the step number. It provides options for customizing the plot type, title, format, file name, saving the plot, showing the plot, and clearing the plot.

    :param values: An array of chi-squared values.
    :type values: numpy.ndarray
    :param discard_number: The number of initial steps to discard.
    :type discard_number: int
    :param use_log_probs: Whether to use log probabilities for calculating chi-squared values. Default is True.
    :type use_log_probs: bool
    :param plot_type: The type of plot to generate. Must be one of "med", "best", or "all". Default is 'med'.
    :type plot_type: str
    :param title: The title of the plot. If not provided, a default title will be used depending on the plot_type.
    :type title: str
    :param no_title: Determines whether to show the title on the plot. Default is False.
    :type no_title: bool
    :param fmt: The format string for the plot. Default is an empty string.
    :type fmt: str
    :param file_name: The name of the file to save the plot. Default is 'med_chi_squared_plot.<image_type>', where <image_type> is the file format specified in the code.
    :type file_name: str
    :param save: Determines whether to save the plot. Default is False.
    :type save: bool
    :param show: Determines whether to show the plot. Default is True.
    :type show: bool
    :param clear_plot: Determines whether to clear the plot before generating a new one. Default is True.
    :type clear_plot: bool
    :return: None
    :rtype: None
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
        x_data = np.arange(discard_number, discard_number+len(chi_sq))
        plt.plot(
            x_data, chi_sq[:], fmt
        )
        fit = create_curve_fit(x_data, chi_sq, 
                               params=[chi_sq[0], 1e4, 100],
                               param_names=[r"$N_0$", r"$\tau$", r"$c$"],
                               fit_type="exponential_decay",
                               )
        plt.plot(x_data, fit[0], 'r-')
        plt.figtext(0.89, 0.80, fit[2],
                    horizontalalignment="right", verticalalignment="top")
        plt.figtext(0.89, 0.855, fit[3], size="large",
                    horizontalalignment="right", verticalalignment="top")

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
    # TODO what happens when folder is None??
    """Plot the range from models within 1 sigma along with the best model and the data.

    .. note::
        - The parameters within 1 sigma that have the biggest and smallest values for each parameter are found, resulting in 2 arrays of dimension NUM_DIM * NUM_DIM.
        - Models are created from these, and for each frequency value, the minimum and the maximum are found.
        - The graph is made by filling in the space between the minimum and maximum for each frequency value.
        - The best model and the actual data with error bars are plotted on top of this.

    :param v_data: The observed values of the frequency.
    :type v_data:
    :param vFv_data: The observed values of the nuFnu.
    :type vFv_data:
    :param err_data: The error values of the data.
    :type err_data:
    :param indices_within_1sigma: The indices of the samples within 1 sigma.
    :type indices_within_1sigma:
    :param flat_samples: The flat samples obtained from the MCMC algorithm.
    :type flat_samples:
    :param min_chi_squared_index: The index of the model with minimum chi-squared value.
    :type min_chi_squared_index:
    :param both: Boolean value indicating whether both extreme and random models should be plotted. Default is False.
    :type both: bool, optional
    :param extreme: Boolean value indicating whether extreme models should be plotted. Default is True.
    :type extreme: bool, optional
    :param title: The title of the plot.
    :type title: str, optional
    :param no_title: If True, do not display a title on the plot.
    :type no_title: bool, optional
    :param folder: The folder in which to save the plot.
    :type folder: str, optional
    :param file: The name of the file to save the plot as.
    :type file: str, optional
    :param save: If True, save the plot to a file.
    :type save: bool, optional
    :param show: If True, display the plot.
    :type show: bool, optional
    :param serialize: Boolean value indicating whether to serialize the plot using pickle. Default is False.
    :type serialize: bool, optional
    :param lower_adjust_multiplier: The lower adjustment multiplier for scaling the y-axis. If None, the minimum y-value of the observed dependent variable will be used.
    :type lower_adjust_multiplier: float, optional
    :param upper_adjust_multiplier: The upper adjustment multiplier for scaling the y-axis. Default is 1.02.
    :type upper_adjust_multiplier: float, optional
    :param max_num_lines_to_graph: The maximum number of lines to graph for random models. Default is 1000.
    :type max_num_lines_to_graph: int, optional
    :param dims: The number of dimensions in the model.
    :type dims: int, optional
    :param eic: If True, include the EIC parameter.
    :type eic: bool, optional
    :param theta: The theta parameter for the model.
    :type theta: float, optional
    :param redshift: The redshift parameter for the model.
    :type redshift: float, optional
    :param min_freq: The minimum frequency for the model.
    :type min_freq: float, optional
    :param max_freq: The maximum frequency for the model.
    :type max_freq: float, optional
    :param executable: The path to the executable for the model.
    :type executable: str, optional
    :param data_folder: The path to the folder containing the data.
    :type data_folder: str, optional
    :param name_stem: The stem of the name of the plot files.
    :type name_stem: str, optional
    :param command_params_full: The command parameters for the model.
    :type command_params_full: str, optional
    :param command_params_1: The first set of command parameters for the model.
    :type command_params_1: str, optional
    :param command_params_2: The second set of command parameters for the model.
    :type command_params_2: str, optional
    :param torus_temp: The torus temperature parameter for the model.
    :type torus_temp: float, optional
    :param torus_luminosity: The torus luminosity parameter for the model.
    :type torus_luminosity: float, optional
    :param torus_frac: The torus fraction parameter for the model.
    :type torus_frac: float, optional
    :param verbose: If True, display verbose output.
    :type verbose: bool, optional
    :return: None
    :rtype: None
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
    # TODO what happens when folder is None??
    """
    .. note::
        - The parameters within 1 sigma that have the biggest and smallest values for each parameter are found, resulting in 2 arrays of dimension NUM_DIM * NUM_DIM.
        - Models are created from these, and for each frequency value, the minimum and the maximum are found.
        - The graph is made by filling in the space between the minimum and maximum for each frequency value.
        - The best model and the actual data with error bars are plotted on top of this.

    :param v_data: The observed values for frequency (ν) data.
    :type v_data: numpy.ndarray
    :param vFv_data: The observed values for νFν data.
    :type vFv_data: numpy.ndarray
    :param err_data: The error values for νFν data.
    :type err_data: numpy.ndarray
    :param indices_within_1sigma: The indices of the samples within 1 sigma.
    :type indices_within_1sigma: numpy.ndarray
    :param flat_samples: The flattened MCMC samples.
    :type flat_samples: numpy.ndarray
    :param min_chi_squared_index: The index of the minimum chi-squared value.
    :type min_chi_squared_index: int
    :param both: Whether to plot both sets of models (within 1 sigma and extreme). Default is False.
    :type both: bool
    :param extreme: Whether to plot extreme models within 1 sigma. Default is True.
    :type extreme: bool
    :param title: The title of the plot. Default is None.
    :type title: str
    :param no_title: Whether to display a title on the plot. Default is False.
    :type no_title: bool
    :param folder: The folder to save the plot in. Default is None.
    :type folder: str
    :param file: The filename to save the plot as. Default is None.
    :type file: str
    :param save: Whether to save the plot. Default is False.
    :type save: bool
    :param show: Whether to display the plot. Default is True.
    :type show: bool
    :param serialize: Whether to serialize the plot as a pickle file. Default is False.
    :type serialize: bool
    :param lower_adjust_multiplier: The lower adjustment multiplier for scaling the y-axis. Default is None.
    :type lower_adjust_multiplier: float
    :param upper_adjust_multiplier: The upper adjustment multiplier for scaling the y-axis. Default is 1.02.
    :type upper_adjust_multiplier: float
    :param max_num_lines_to_graph: The maximum number of lines to plot. Default is 1000.
    :type max_num_lines_to_graph: int
    :param dims: The dimensions of the model. Default is None.
    :type dims: int
    :param eic: Whether to use Extended Isothermal Cone (EIC) model. Default is False.
    :type eic: bool
    :param return_models: Whether to return the generated models. Default is False.
    :type return_models: bool
    :param theta: The theta value for the model. Default is None.
    :type theta: float
    :param redshift: The redshift value for the model. Default is None.
    :type redshift: float
    :param min_freq: The minimum frequency value for the model. Default is None.
    :type min_freq: float
    :param max_freq: The maximum frequency value for the model. Default is None.
    :type max_freq: float
    :param executable: The path to the executable file. Default is None.
    :type executable: str
    :param data_folder: The path to the data folder. Default is None.
    :type data_folder: str
    :param command_params_full: The full command parameters. Default is None.
    :type command_params_full: str
    :param command_params_1: The first command parameters. Default is None.
    :type command_params_1: str
    :param command_params_2: The second command parameters. Default is None.
    :type command_params_2: str
    :param name_stem: The name stem for the plots. Default is None.
    :type name_stem: str
    :param torus_temp: The temperature of the torus. Default is None.
    :type torus_temp: float
    :param torus_luminosity: The luminosity of the torus. Default is None.
    :type torus_luminosity: float
    :param torus_frac: The fractional size of the torus. Default is None.
    :type torus_frac: float
    :param verbose: Whether to display verbose output. Default is False.
    :type verbose: bool
    :return: The generated models, if return_models is True.
    :rtype: list
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


    :param values: The input values to be scaled.
    :type values: list, numpy.ndarray
    :param upper_adjust_multiplier: Optional. The multiplier to adjust the upper range.
        If not provided, the default value is 5.
    :type upper_adjust_multiplier: int or float
    :param lower_adjust_multiplier: Optional. The multiplier to adjust the lower range.
        If not provided, the default value is 5.
    :type lower_adjust_multiplier: int or float
    :return: The scaled minimum and maximum values.
    :rtype: tuple

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
    Returns the minimum and maximum parameter values within the specified range.

    Finds the array of parameters that has the minimum and maximum value for
    each of the parameters.

    For example, with samples [[0, 1, 2], [5, 2, 1], [4, 6, 1], [3, 2, 0]], the
    minima would be [[0, 1, 2], [0, 1, 2], [3, 2, 0]] and the maxima would be
    [[5, 2, 1], [4, 6, 1], [0, 1, 2]]

    :param flat_samples: A numpy array of flattened parameter samples.
    :type flat_samples: numpy.ndarray
    :param indices_within_1sigma: A list of indices representing samples within the desired range.
    :type indices_within_1sigma: list
    :param eic: A boolean value indicating whether the model properties should be adjusted for electron-ion collisions. Default is False.
    :type eic: bool
    :param fixed_params: A dictionary of fixed parameters. Default is None.
    :type fixed_params: dict
    :return: A tuple containing the minimum and maximum parameter values within the specified range.
    :rtype: tuple
    """
    # set the minima and maxima to the first set of params for all vals
    dims = modelProperties(eic, fixed_params=fixed_params).NUM_DIM
    # print(np.shape(indices_within_1sigma))
    # print(eic)
    minima = [flat_samples[indices_within_1sigma[0]].copy() for _ in range(dims)]
    maxima = [flat_samples[indices_within_1sigma[0]].copy() for _ in range(dims)]
    index_min = [0]*dims
    index_max = [0]*dims
    for index in indices_within_1sigma:
        params = flat_samples[index]
        # print(len(params), len(minima))
        for i in range(dims):
            if params[i] < minima[i][i]:
                minima[i] = params.copy()
                index_min[i] = index
            if params[i] > maxima[i][i]:
                maxima[i] = params.copy()
                index_max[i] = index
    
    #we have sometime multiple time the same index, this doesn't optimize the contours as it runs multiple time the exact same parameters
    #just select index-1 for the duplicate indexes
    # all_indexes = index_min+index_max
    # for i in range(dims):
    #     if index_min.count(index_min[i]) > 1 or index_min[i] in index_max:
    #         if index_min[i] < len()
    #         index_min[i] += 1
            
    
    
    #print(index_min, index_max)
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
    :param v_vals: The array of values at which to calculate the minimum and maximum.
    :type v_vals: numpy.ndarray
    :param model_params_list: The list of model parameter arrays from which to calculate the minimum and maximum.
    :type model_params_list: List[numpy.ndarray]
    :param name_stem: The prefix for the name of the output files (optional).
    :type name_stem: str
    :param theta: The viewing angle of the model (optional).
    :type theta: float
    :param redshift: The redshift of the source (optional).
    :type redshift: float
    :param min_freq: The minimum frequency of the SED (optional).
    :type min_freq: float
    :param max_freq: The maximum frequency of the SED (optional).
    :type max_freq: float
    :param torus_temp: The temperature of the torus (optional).
    :type torus_temp: float
    :param torus_luminosity: The luminosity of the torus (optional).
    :type torus_luminosity: float
    :param torus_frac: The fraction of the disk luminosity contributed by the torus (optional).
    :type torus_frac: float
    :param data_folder: The folder where the input and output files are stored (optional).
    :type data_folder: str
    :param executable: The path to the executable for the model calculation (optional).
    :type executable: str
    :param command_params_full: The command-line parameters for the full model calculation (optional).
    :type command_params_full: str
    :param command_params_1: The command-line parameters for the first stage model calculation (optional).
    :type command_params_1: str
    :param command_params_2: The command-line parameters for the second stage model calculation (optional).
    :type command_params_2: str
    :param verbose: Whether to display verbose output (default is False).
    :type verbose: bool
    :param eic: Whether to enable emission-inverse Compton (EIC) calculations (default is False).
    :type eic: bool
    :param fixed_params: The array of fixed model parameters to be inserted into the model calculations (optional).
    :type fixed_params: List[float]
    :return: The arrays of minimum and maximum values per point calculated from the given model parameter arrays.
    :rtype: Tuple[numpy.ndarray, numpy.ndarray]
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
    This function generates a residual plot for given data and model parameters.

    :param data: The data used for generating the residual plot.
    :type data: tuple or list
    :param best_model: The best-fit model used for generating the plot.
    :type best_model: tuple or list
    :param lowest_per_point: The lower error bound for each data point.
    :type lowest_per_point: tuple or list
    :param highest_per_point: The upper error bound for each data point.
    :type highest_per_point: tuple or list
    :return: None
    :rtype: None

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
    Calculate the effective number of galaxies, N_e, per bin for a broken power law luminosity function.

    :param gamma: The value of gamma parameter.
    :type gamma: float
    :param K: The value of K parameter.
    :type K: float
    :param n_1: The value of n_1 parameter.
    :type n_1: float
    :param n_2: The value of n_2 parameter.
    :type n_2: float
    :param gamma_break: The value of gamma_break parameter.
    :type gamma_break: float
    :return: The calculated value of N_e using the broken power law formula.
    :rtype: float
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

    Array of best parameters found by the MCMC fit.
    For simple SSC without fixed parameters, the order is: [delta, K, n_1, n_2, gamma_min, gamma_max, gamma_break, B, R]
    Additional description:

    .. list-table::
       :header-rows: 1

       * - Parameter
         - Description
         - Scale
       * - ``delta``
         - Doppler factor
         - Linear
       * - ``K``
         - Particle density [cm^-3]
         - Log
       * - ``n1``
         - alpha_1 (first index)
         - Linear
       * - ``n2``
         - alpha_2 (second index)
         - Linear
       * - ``gamma_min``
         - Low-energy cutoff
         - Log
       * - ``gamma_max``
         - High-energy cutoff
         - Log
       * - ``gamma_break``
         - Energy break
         - Log
       * - ``B``
         - Magnetic field strength [G]
         - Log
       * - ``R``
         - Blob radius [cm]
         - Log

    Additional params for EIC

    .. list-table::
       :header-rows: 1

       * - Parameter
         - Description
         - Scale
       * - ``bb_temp``
         - Black body temp of disk [K]
         - Log
       * - ``l_nuc``
         - Nucleus luminosity [ergs/s]
         - Log
       * - ``tau``
         - Fraction of luminosity scattered
         - Log
       * - ``blob_dist``
         - Distance of blob [cm]
         - Log

    Note the errors of gamma_min and gamma_max are not included in the contour


    :param best_params: The best-fit parameters for the particle spectrum.
    :type best_params: numpy array (shape: (6,))
    :param min_1sigma_params: Array of 1 sigma lower boundary of free parameters found by the MCMC fit. The order should match the one of best_params
    :type min_1sigma_params: numpy array (shape: (6,))
    :param max_1sigma_params:  Array of 1 sigma upper boundary of free parameters found by the MCMC fit. The order should match the one of best_params
    :type max_1sigma_params: numpy array (shape: (6,))
    :param fixed_params: List of user-fixed parameters. If no fixed parameters in simple SSC model: fixed_params = [-inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf]
    :type fixed_params: numpy array (shape: (9,))
    :param file_name: Absolute path and name for the saved plot.The default is BASE_PATH + RESULTS_FOLDER + "/particle_spectrum.svg".
    :type file_name: str
    :param save: Flag to save the plot of the particle spectrum. Default is False.
    :type save: bool
    :param show: Flag to show the plot of the particle spectrum. Default is True.
    :type show: bool
    :return: None
    :rtype: None

    This function plots the particle spectrum based on the given best-fit parameters, lower and upper bounds of the 1-sigma confidence interval, and fixed parameters. The plot is saved as an SVG file if the "save" flag is set to True. The plot is also displayed if the "show" flag is set to True.
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
    Calculates the cooling time for a given particle energy using the Thomson scattering process. (see e.g. [Inoue]_)

        :param gamma: The Lorentz factor of the particle.
        :type gamma: float
        :param U_B: The energy density of the magnetic field.
        :type U_B: float
        :param U_syn: The energy density of the synchrotron radiation.
        :type U_syn: float
        :param U_blr: The energy density of the broad-line region.
        :type U_blr: float
        :return: The cooling time of the particle according to Thomson scattering.
        :rtype: float

    .. [Inoue] Inoue & Takahara 1996
    """

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
    Plots the observed cooling times in the Thomson limit.

    :param logfile: The path to the log file.
    :type logfile: str
    :param best_params: A list of best fit parameters.
    :type best_params: list
    :param fixed_params: A list of fixed parameters.
    :type fixed_params: list
    :param file_name: The name of the output file.
    :type file_name: str
    :param save: Whether to save the plot as an SVG file.
    :type save: bool
    :param show: Whether to display the plot.
    :type show: bool
    :param eic: Whether to grep the energy density from bjet.log.
    :type eic: bool
    :param redshift: The redshift value.
    :type redshift: float or None
    :return: None
    :rtype: None
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
    :param flat_samples: The flattened samples of the parameters obtained from the MCMC sampling.
    :type flat_samples: array-like
    :param flat_log_probs: The flattened log-probability values corresponding to the samples.
    :type flat_log_probs: array-like
    :param best_params: The best-fit parameter values obtained from the MCMC sampling.
    :type best_params: array-like
    :param min_1sigma_params: The parameter values at the lower 1-sigma confidence interval obtained from the MCMC sampling.
    :type min_1sigma_params: array-like
    :param max_1sigma_params: The parameter values at the upper 1-sigma confidence interval obtained from the MCMC sampling.
    :type max_1sigma_params: array-like
    :param save: Flag indicating whether to save the plotted likelihood profiles. Default is False.
    :type save: bool, optional
    :param show: Flag indicating whether to display the plotted likelihood profiles. Default is True.
    :type show: bool, optional
    :param fixed_params: The fixed parameter values that are not varied during the MCMC sampling. Default is None.
    :type fixed_params: dict, optional
    :param eic: Flag indicating whether the model is an EIC model. Default is False.
    :type eic: bool, optional
    :param folder_path: The path to the folder where the plot will be saved. Default is BASE_PATH + RESULTS_FOLDER.
    :type folder_path: str, optional
    :return: None
    :rtype: None

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
        #automatically rescale x axis if the parameter space is too narrow
        custom_low = best_params[j]-5*(best_params[j]-min_1sigma_params[j])
        custom_high = best_params[j]+5*(max_1sigma_params[j]-best_params[j])
        if param_array[0] < custom_low:
            min_val = custom_low
        else:
            min_val = param_array[0]
        if param_array[-1] > custom_high:
            max_val = custom_high
        else:
            max_val = param_array[-1]
        Drange = max_val - min_val
        binsize = Drange / (nbins - 1)
        bin_min = min_val
        bin_max = bin_min + binsize
        prob_tmp = []
        prob_max = []
        ln_prob_max = []
        param_binned_mid = []
        for i in range(len(param_array)):
            if min_val <= param_array[i] <= max_val:
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
        plt.xlim(min_val, max_val)
        plt.ylim(0,max(prob_max)*1.2)

        if save:
            plt.savefig(folder_path + "/" + figure_name + ".svg")
        if show:
            plt.show()
