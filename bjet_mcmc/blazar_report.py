#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: blazar_report.py

Contains functions for showing results of the MCMC, given the results. It
calculates indices and parameters within 1sigma, finds the best parameters,
creates a written text report, and creates all necessary plots.

Functions:
Function save_plots_and_info (given the configurations, the data, and
the minima/maxima for the parameters, and a source of MCMC data) creates a folder
with an info.txt file and plots:
    - Corner plot
    - Plots of chi squared (med, best, all by step)
    - Plots of the best model with data and models within 1 sigma (random
      models, extreme models, and a combination of both)
    - Plot of the chain

Parameters are always listed in the following order:
[delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R]
All parameters are the logarithm of the true value except for delta, n1, and n2
---------------------------------------------------
delta       doppler factor                  linear
K           particle density [cm^-3]        log
n1          alpha_1 (first index)           linear
n2          alpha_2 (second index)          linear
gamma_min   low-energy cutoff               log
gamma_max   high-energy cutoff              log
gamma_break energy break                    log
B           magnetic field strength [G]     log
R           blob radius (cm)                log
---------------------------------------------------
"""
import glob
import os
import random

import emcee
import numpy as np
from scipy import stats

from bjet_mcmc import blazar_model
from bjet_mcmc import blazar_plots
from bjet_mcmc import blazar_utils
from bjet_mcmc.blazar_properties import *

__all__ = [
    "get_best_log_prob_and_params",
    "get_indices_within_1sigma",
    "load_from_backend",
    "make_text_file",
    "min_max_params_1sigma",
    "parse_info_doc",
    "save",
    "save_plots_and_info",
    "show_results",
    "text_report",
]
# utils ----------------------------------------------------------------------------------------------------------------


def get_indices_within_1sigma(
    values, use_log_probs=True, alpha=0.318, eic=False, ndim=1
):
    """
    Args:
        values: An array of values for which chi-squared values will be calculated.
        use_log_probs: A boolean indicating whether to use log probabilities for calculating chi-squared values. Defaults to True.
        alpha: The significance level used for calculating delta chi-squared. Defaults to 0.318.
        eic: A boolean indicating whether to use modelProperties(eic).NUM_DIM for calculating delta chi-squared. Defaults to False.
        ndim: The number of dimensions used for calculating delta chi-squared. Defaults to 1.

    Returns:
        indices_within_1sigma: An array containing the indices of values that lie within 1 sigma (delta chi-squared) of the minimum chi-squared value.
    """
    if use_log_probs:
        chi_squared = -2 * values
    else:
        chi_squared = 1 * values
    min_chi_squared_index = np.argmin(chi_squared)
    min_chi_squared = chi_squared[min_chi_squared_index]
    delta_chi_squared = stats.chi2.ppf(1 - alpha, ndim)  # modelProperties(eic).NUM_DIM)
    indices_within_1sigma = np.where(chi_squared < min_chi_squared + delta_chi_squared)[
        0
    ]
    return indices_within_1sigma


def min_max_params_1sigma(flat_data, flat_log_probs, eic=False, ndim=1):
    """
    Args:
        flat_data (List[List[float]]): A list of data points, where each data point is a list of floats.
        flat_log_probs (List[float]): A list of log probabilities corresponding to each data point in flat_data.
        eic (bool, optional): Whether to consider events inside the 1-sigma contour. Default is False.
        ndim (int, optional): Number of dimensions of the data points. Default is 1.

    Returns:
        Tuple[List[float], List[float]]: A tuple containing the minimum values and maximum values within the 1-sigma contour.

    """
    n_dim = len(flat_data[0])
    indices_within_1sigma = get_indices_within_1sigma(
        flat_log_probs, eic=eic, ndim=ndim
    )
    mins = np.copy(flat_data[indices_within_1sigma[0]])
    maxes = np.copy(flat_data[indices_within_1sigma[0]])
    for i in range(len(indices_within_1sigma)):
        for j in range(n_dim):
            val = flat_data[indices_within_1sigma[i]][j]
            if val < mins[j]:
                mins[j] = val
            if val > maxes[j]:
                maxes[j] = val
    return mins, maxes


def get_best_log_prob_and_params(
    configs=None,
    sampler=None,
    log_probs=None,
    flat_chain=None,
    discard=None,
    fixed_params=None,
):
    """
    Args:
        configs (dict, optional): Configuration options. Default is None.
        sampler (object, optional): The sampler object. Default is None.
        log_probs (array-like, optional): Array of log probabilities. Default is None.
        flat_chain (array-like, optional): Array of parameter samples. Default is None.
        discard (int, optional): The number of discarded samples. Default is None.
        fixed_params (array-like, optional): Array of fixed parameter values. Default is None.

    Returns:
        tuple: A tuple containing the maximum log probability and the corresponding best parameters.

    Raises:
        ValueError: If either the sampler is not provided or both log_probs and flat_chain are not provided.
        ValueError: If the sampler is provided but either the configs or discard is not provided.

    """
    if (sampler is None) and (log_probs is None or flat_chain is None):
        raise ValueError("sampler or log_probs+chain not provided")
    if sampler is not None:
        if configs is None and discard is None:
            raise ValueError("sampler provided but no configs or discard")
        if discard is None:
            discard = configs["discard"]
        log_probs = sampler.get_log_prob(flat=True, discard=discard)
        flat_chain = sampler.get_chain(flat=True, discard=discard)
    maximum_index = np.argmax(log_probs)
    maximum_log_prob = log_probs[maximum_index]
    best_params = flat_chain[maximum_index]
    if fixed_params:
        for i in range(len(fixed_params)):
            if fixed_params[i] != -np.inf:
                best_params[i] = fixed_params[i]
    return maximum_log_prob, best_params


# report displays ------------------------------------------------------------------------------------------------------
def make_text_file(
    output_file,
    configs,
    data_points,
    backend_file=None,
    reader=None,
    sampler=None,
    use_sampler=False,
    samples=None,
    use_samples=False,
    description=None,
    time=None,
    p0_source=None,
    acceptance_frac=None,
    eic=False,
    fixed_params=None,
):
    """
    Args:
        output_file (str): The name of the output file to write the text report to.
        configs (dict): A dictionary containing the configurations for the analysis.
        data_points ():
        backend_file (str): The name of the backend file to read the samples from. Default is None.
        reader (): The reader object to use for reading the samples. Default is None.
        sampler (): The sampler object to use for sampling the model parameters. Default is None.
        use_sampler (bool): Whether to use the provided sampler object. Default is False.
        samples (): The samples to write to the output file. Default is None.
        use_samples (bool): Whether to use the provided samples. Default is False.
        description (): The description for the text report. Default is None.
        time (): The time for the text report. Default is None.
        p0_source (): The source of the initial parameters for the text report. Default is None.
        acceptance_frac (): The acceptance fraction for the samples. Default is None.
        eic (bool): Whether to use EIC (Effective Information Criteria) for min/max params calculation. Default is False.
        fixed_params (): The fixed parameters for the model. Default is None.
    """
    if use_samples and (samples is not None):
        samples, flat_chain, _, log_probs = samples
    # note that samplers and readers both have the functions needed
    else:
        if use_sampler and (sampler is not None):
            data_source = sampler
        elif backend_file is not None:
            data_source = emcee.backends.HDFBackend(
                BASE_PATH + backend_file, read_only=True
            )
        elif reader is not None:
            data_source = reader
        else:
            raise ValueError("backend file or reader not provided")

        samples = data_source.get_chain(discard=configs["discard"])
        flat_chain = data_source.get_chain(flat=True, discard=configs["discard"])
        log_probs = data_source.get_log_prob(flat=True, discard=configs["discard"])

    with open(BASE_PATH + output_file, "w") as f:
        tau = emcee.autocorr.integrated_time(samples, quiet=True)
        maximum_log_prob, best_params = get_best_log_prob_and_params(
            log_probs=log_probs, flat_chain=flat_chain
        )
        print("best parameters:", best_params)
        min_1sigma_params, max_1sigma_params = min_max_params_1sigma(
            flat_chain, log_probs, eic=eic, ndim=1
        )
        min_1sigma_params_SED, max_1sigma_params_SED = min_max_params_1sigma(
            flat_chain, log_probs, eic=eic, ndim=modelProperties(eic).NUM_DIM
        )
        if description is not None:
            f.write("report description: ")
            f.write(str(description))
            f.write("\n")
        if time is not None:
            f.write("time: ")
            f.write(time)
            f.write("\n")
        if p0_source is not None:
            f.write("p0 from ")
            f.write(p0_source)
            f.write("\n\n")

        f.write("configurations: \n")
        f.write(str(configs))
        f.write("\n\n")

        f.write(
            text_report(
                best_params,
                min_1sigma_params,
                max_1sigma_params,
                min_1sigma_params_SED,
                max_1sigma_params_SED,
                -2 * maximum_log_prob,
                data_points,
                eic=eic,
                fixed_params=fixed_params,
            )
        )
        f.write("\n\n")

        f.write("best params: ")
        f.write(str(best_params))
        f.write(", chi squared = ")
        f.write(str(-2 * maximum_log_prob))
        f.write("\n\n")

        f.write("min_1sigma_params: ")
        f.write(str(min_1sigma_params))
        f.write("\nmax_1sigma_params: ")
        f.write(str(max_1sigma_params))
        f.write("\n\n")

        if acceptance_frac is not None:
            f.write("acceptance fraction: avg = ")
            f.write(str(np.average(acceptance_frac)))
            f.write("\n\n")

        f.write("autocorrelation time: avg = ")
        f.write(str(np.mean(tau)) + " steps")
        f.write("\n\n")


def text_report(
    best,
    min_1sigma,
    max_1sigma,
    min_1sigma_SED,
    max_1sigma_SED,
    best_chi_sq,
    data_points,
    param_names=None,
    eic=False,
    fixed_params=None,
):
    """
    Generates a text report based on the given parameters.

    Args:
        best ():
        min_1sigma ():
        max_1sigma ():
        min_1sigma_SED ():
        max_1sigma_SED ():
        best_chi_sq ():
        data_points ():
        param_names ():
        eic ():
        fixed_params ():

    Returns:
        str: The generated text report.
    """
    if param_names is None:
        param_names = modelProperties(eic, fixed_params=fixed_params).PARAM_NAMES
    log = modelProperties(eic).PARAM_IS_LOG
    # remove any frozen parameter from the log list
    if fixed_params:
        fixed_params2 = fixed_params.copy()
        i = 0
        while i < len(fixed_params2):
            if fixed_params2[i] != -np.inf:
                del log[i]
                del fixed_params2[i]
            else:
                i += 1
    format_string = "{:^13}{:^15}{:^20}{:^40}\n"
    num_format = "{:.2e}"
    range_format = num_format + " - " + num_format
    output = format_string.format(
        "Parameter",
        "Best Value",
        "1-Sigma Range",
        "Extremums within 1-Sigma SED contour",
    )
    for i in range(len(best)):
        b = best[i]
        s1 = min_1sigma[i]
        s2 = max_1sigma[i]
        s1_SED = min_1sigma_SED[i]
        s2_SED = max_1sigma_SED[i]
        if log[i]:
            b = np.power(10, b)
            s1 = np.power(10, s1)
            s2 = np.power(10, s2)
            s1_SED = np.power(10, s1_SED)
            s2_SED = np.power(10, s2_SED)
        output += format_string.format(
            param_names[i],
            num_format.format(b),
            range_format.format(s1, s2),
            range_format.format(s1_SED, s2_SED),
        )

    # dims = modelProperties(eic).NUM_DIM
    dims = len(best)
    output += "Reduced chi squared: {:.2f} / {} = {:.2f}\n".format(
        best_chi_sq, data_points - dims, best_chi_sq / (data_points - dims)
    )

    return output


def show_results(sampler, time, configs=None, discard=None):
    """
    This method displays results from an MCMC sampler.

    Args:
        sampler (): The MCMC sampler object.
        time (): The total duration of the sampling process.
        configs (): Optional configuration parameters for the sampler. Default is None.
        discard (): Optional number of discarded samples. Default is None.

    """
    if configs is None and discard is None:
        raise ValueError("Neither configs nor discard provided")
    if discard is None:
        discard = configs["discard"]
    chain = sampler.get_chain(discard=discard)
    tau = emcee.autocorr.integrated_time(chain, quiet=True)
    maximum_log_prob, best_params = get_best_log_prob_and_params(
        sampler=sampler, discard=discard
    )

    print("MCMC Results:")
    print("Best log_prob:", maximum_log_prob)
    print("Best parameters:", best_params)
    print("Time:", time)
    print("Mean auto-correlation time:", np.mean(tau))


def save_plots_and_info(
    configs,
    data,
    param_min_vals,
    param_max_vals,
    folder=None,
    parent_folder=None,
    backend_file=None,
    reader=None,
    sampler=None,
    use_sampler=False,
    samples=None,
    use_samples=False,
    description=None,
    time=None,
    p0_source=None,
    acceptance_frac=None,
    theta=None,
    redshift=None,
    min_freq=None,
    max_freq=None,
    executable=None,
    data_folder=None,
    command_params_full=None,
    command_params_1=None,
    command_params_2=None,
    torus_temp=None,
    torus_luminosity=None,
    torus_frac=None,
    verbose=False,
    eic=False,
):
    """
    Args:
        configs ():
        data ():
        param_min_vals ():
        param_max_vals ():
        folder ():
        parent_folder ():
        backend_file ():
        reader ():
        sampler ():
        use_sampler ():
        samples ():
        use_samples ():
        description ():
        time ():
        p0_source ():
        acceptance_frac ():
        theta ():
        redshift ():
        min_freq ():
        max_freq ():
        executable ():
        data_folder ():
        command_params_full ():
        command_params_1 ():
        command_params_2 ():
        torus_temp ():
        torus_luminosity ():
        torus_frac ():
        verbose ():
        eic ():
    """
    if (
        backend_file is None
        and reader is None
        and (not use_sampler or sampler is None)
        and not (use_samples or samples is None)
    ):  # source not given
        raise ValueError("backend file or reader not provided")
    if reader is None and not use_sampler and not use_samples:
        reader = emcee.backends.HDFBackend(BASE_PATH + backend_file, read_only=True)
    if folder is None and parent_folder is None:
        raise ValueError("folder or parent folder not provided")
    if folder is None:
        folder = parent_folder + "/new_report"
        os.mkdir(folder)

    if use_sampler:
        chain = sampler.get_chain(discard=configs["discard"])
        flat_chain = sampler.get_chain(flat=True, discard=configs["discard"])
        log_probs = sampler.get_log_prob(discard=configs["discard"])
        flat_log_probs = sampler.get_log_prob(flat=True, discard=configs["discard"])
        acceptance_frac = sampler.acceptance_fraction
    elif use_samples:
        chain, flat_chain, log_probs, flat_log_probs = samples
    else:
        chain = reader.get_chain(discard=configs["discard"])
        flat_chain = reader.get_chain(flat=True, discard=configs["discard"])
        log_probs = reader.get_log_prob(discard=configs["discard"])
        flat_log_probs = reader.get_log_prob(flat=True, discard=configs["discard"])

    best_log_prob, best_params = get_best_log_prob_and_params(
        log_probs=flat_log_probs, flat_chain=flat_chain
    )
    min_1sigma_params, max_1sigma_params = min_max_params_1sigma(
        flat_chain, flat_log_probs, eic=eic
    )
    indices_within_1sigma = get_indices_within_1sigma(
        flat_log_probs, eic=eic, ndim=modelProperties(eic).NUM_DIM
    )
    name_stem = NAME_STEM + str(random.getrandbits(20))

    command_params_1, command_params_2 = blazar_model.command_line_sub_strings(
        name_stem=name_stem, redshift=redshift, prev_files=False, eic=eic
    )
    command_params_2[3] = "99"  # number of points used to make SED

    model = blazar_model.make_model(
        best_params,
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
        prev_files=False,
        use_param_file=False,
        verbose=True,
        eic=eic,
        fixed_params=configs["fixed_params"],
        folder=folder,
    )

    for f in glob.glob(BASE_PATH + DATA_FOLDER + "/" + name_stem + "_*"):
        os.remove(f)
    make_text_file(
        folder + "/info.txt",
        configs,
        len(data[0]),
        backend_file=backend_file,
        reader=reader,
        sampler=sampler,
        use_sampler=use_sampler,
        samples=samples,
        use_samples=use_samples,
        description=description,
        time=time,
        p0_source=p0_source,
        acceptance_frac=acceptance_frac,
        eic=eic,
        fixed_params=configs["fixed_params"],
    )
    blazar_plots.plot_model_and_data(
        model,
        configs["data_file"],
        flat_chain,
        indices_within_1sigma,
        redshift=redshift,
        eic=eic,
        lower_adjust_multiplier=20,
        file_name=folder + "/model_and_data.svg",
        title=None,
        no_title=True,
        save=True,
        show=False,
        fixed_params=configs["fixed_params"],
    )

    # need to run tests with fixed params before releasing it
    blazar_plots.plot_particle_spectrum(
        best_params,
        min_1sigma_params,
        max_1sigma_params,
        configs["fixed_params"],
        file_name=BASE_PATH + folder + "/particle_spectrum.svg",
        save=True,
        show=False,
    )

    blazar_plots.corner_plot(
        flat_chain,
        param_min_vals,
        param_max_vals,
        best_params,
        min_1sigma_params,
        max_1sigma_params,
        file_name=folder + "/corner_plot.svg",
        save=True,
        show=False,
        eic=eic,
        fixed_params=configs["fixed_params"],
    )
    blazar_plots.plot_likelihood_profiles(
        flat_chain,
        flat_log_probs,
        best_params,
        min_1sigma_params,
        max_1sigma_params,
        save=True,
        show=False,
        fixed_params=configs["fixed_params"],
        eic=eic,
        folder_path=BASE_PATH + folder,
    )
    # This output is not easy to read, so not fully relevant. removed for now
    # blazar_plots.plot_chain(chain, file_name=(folder + "/plot_of_chain.jpeg"), save=True, show=False,
    #                         eic=eic)  # too much stuff for svg
    blazar_plots.plot_chi_squared(
        log_probs,
        configs["discard"],
        plot_type="med",
        file_name=(folder + "/chi_squared_plot_med.svg"),
        save=True,
        show=False,
    )
    blazar_plots.plot_chi_squared(
        log_probs,
        configs["discard"],
        plot_type="best",
        file_name=(folder + "/chi_squared_plot_best.svg"),
        save=True,
        show=False,
    )
    blazar_plots.plot_chi_squared(
        log_probs,
        configs["discard"],
        plot_type="all",
        file_name=(folder + "/chi_squared_plot_all.jpeg"),
        save=True,
        show=False,
    )  # svg is big

    blazar_plots.plot_cooling_times(
        folder + "/bjet.log",
        best_params,
        fixed_params=configs["fixed_params"],
        file_name=folder + "/cooling_time_obs(Thomson).svg",
        save=True,
        show=False,
        eic=eic,
        redshift=redshift,
    )


def parse_info_doc(info_doc, info=None):
    """
    Args:
        info_doc (str): The filename of the info documentation file.
        info (dict, optional): Additional information dictionary to populate. Defaults to None.

    Returns:
        dict: Populated information dictionary.

    Description:
    This method parses the content of an info documentation file and extracts relevant information into a dictionary. The parsed information can include the folder name, report description, time, p0 label, configurations, reduced chi squared, best parameters, minimum 1 sigma parameters, maximum 1 sigma parameters, acceptance fraction, tau average, config file, and previous files.

    To use this method, provide the filename of the info documentation file as the info_doc parameter. Optionally, you can provide an existing dictionary to populate or leave it as None to create a new one. The method returns the populated information dictionary.

    Example Usage:
        info = parse_info_doc("info.txt")
    """
    if info is None:
        info = {}

    labels = [
        "folder name:",
        "report description:",
        "time:",
        "p0 ",
        "configurations:",
        "Reduced chi squared:",
        "best params:",
        "min_1sigma_params:",
        "max_1sigma_params:",
        "acceptance fraction:",
        "tau: avg",
        "config file:",
        "prev_files:",
    ]

    def match(to_match):
        for l in labels:
            if to_match.strip().find(l) == 0:
                return l
        else:
            return None

    with open(BASE_PATH + info_doc, "r") as f:
        for line in f:
            matched = match(line)
            if matched is None:
                continue
            elif matched == "folder name:":
                info["folder_name"] = line.strip()[len("folder name:") + 1 :]
            elif matched == "report description:":
                info["description"] = line.strip()[len("report description:") + 1 :]
            elif matched == "time:":
                info["time"] = line.strip()[6:]
            elif matched == "p0 ":
                info["p0_label"] = line.strip()[8:]
            elif matched == "configurations:":
                line2 = f.readline()
                info["configs"] = blazar_utils.read_configs(config_string=line2)
            elif matched == "Reduced chi squared:":
                info["reduced_chi_sq"] = float(line.strip().split()[-1].strip())
            elif matched == "best params:":
                best_params_string_list = line[
                    line.strip().find("[") + 1 : line.rfind("]")
                ].split()
                best_params = []
                for elem in best_params_string_list:
                    best_params.append(float(elem.strip()))
                info["best_params"] = np.array(best_params)
                info["best_chi_sq"] = float(line.strip().split()[-1].strip())
            elif matched == "min_1sigma_params:":
                min_1sigma_string_list = line[
                    line.strip().find("[") + 1 : line.rfind("]")
                ].split()
                min_1sigma_params = []
                for elem in min_1sigma_string_list:
                    min_1sigma_params.append(float(elem.strip()))
                info["min_1sigma_params"] = np.array(min_1sigma_params)
            elif matched == "max_1sigma_params:":
                max_1sigma_string_list = line[
                    line.strip().find("[") + 1 : line.rfind("]")
                ].split()
                max_1sigma_params = []
                for elem in max_1sigma_string_list:
                    max_1sigma_params.append(float(elem.strip()))
                info["max_1sigma_params"] = np.array(max_1sigma_params)
            elif matched == "acceptance fraction:":
                info["acceptance_frac"] = float(line.strip().split()[-1].strip())
            elif matched == "tau: avg":
                info["tau_avg"] = float(line.strip().split()[-1].strip())
            elif matched == "config file:":
                info["config_file"] = line.strip()[len("config file:") + 1 :]
            elif matched == "prev_files:":
                temp = line.strip().split()
                info["prev_files"] = temp[1].strip()[0] == "T"
                if len(temp) >= 4:
                    info["use_param_files"] = temp[3].strip()[0] == "T"
    return info


# will create all results from just a folder w/ a backend file in it
def save(
    folder,
    description=None,
    time=None,
    p0_source=None,
    acceptance_frac=None,
    data_file=None,
    configs=None,
    text_file_name="info.txt",
    text_only=False,
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
        folder (str): The folder where the data and plots will be saved.
        description (str): The description of the data.
        time (str): The time of the data.
        p0_source (str): The p0 label of the data.
        acceptance_frac (str): The acceptance fraction of the data.
        data_file (str): The path to the data file.
        configs (dict): The configurations for the data.
        text_file_name (str): The name of the text file to be generated.
        text_only (bool): If True, only a text file will be generated.
        theta (float): The theta value.
        redshift (float): The redshift value.
        min_freq (float): The minimum frequency value.
        max_freq (float): The maximum frequency value.
        torus_temp (float): The torus temperature value.
        torus_luminosity (float): The torus luminosity value.
        torus_frac (float): The torus fraction value.
        data_folder (str): The folder where the data is stored.
        executable (str): The path to the executable file.
        command_params_full (str): The full command parameters.
        command_params_1 (str): The first command parameters.
        command_params_2 (str): The second command parameters.
        verbose (bool): If True, additional information will be displayed.
        eic (bool): If True, additional EIC information will be considered.
        fixed_params (dict): The fixed parameters for the data.

    """
    # WILL USE CURRENT CONFIGURATIONS IF NONE FOUND
    info = {}
    if os.path.exists(BASE_PATH + folder + "/info.txt"):
        info = parse_info_doc(folder + "/info.txt")
    elif os.path.exists(BASE_PATH + folder + "/basic_info.txt"):
        info = parse_info_doc(folder + "/basic_info.txt")
    # get data
    if "configs" in info:
        configs = info["configs"]
        data = blazar_utils.read_data(configs["data_file"])
        if "eic" in configs:  # not true for some old files
            eic = configs["eic"]
    else:
        if configs is None:
            configs = blazar_utils.read_configs()
        if data_file is not None:
            data = blazar_utils.read_data(data_file)
        else:
            data = blazar_utils.read_data(configs["data_file"])

    if "time" in info:
        time = info["time"]
    if "p0_label" in info:
        p0_source = info["p0_label"]
    if "acceptance_frac" in info:
        acceptance_frac = info["acceptance_frac"]
    if "description" in info:
        description = info["description"]
    param_min_vals, param_max_vals = blazar_utils.min_max_parameters(
        alpha2_limits=configs["alpha2_limits"], eic=eic
    )

    if text_only:
        make_text_file(
            folder + "/" + text_file_name,
            configs,
            len(data[0]),
            backend_file=folder + "/backend.h5",
            description=description,
            time=time,
            p0_source=p0_source,
            acceptance_frac=acceptance_frac,
            eic=eic,
            fixed_params=fixed_params,
        )
    else:
        save_plots_and_info(
            configs,
            data,
            param_min_vals,
            param_max_vals,
            folder=folder,
            backend_file=folder + "/backend.h5",
            description=description,
            time=time,
            p0_source=p0_source,
            acceptance_frac=acceptance_frac,
            theta=theta,
            redshift=redshift,
            min_freq=min_freq,
            max_freq=max_freq,
            executable=executable,
            data_folder=data_folder,
            command_params_full=command_params_full,
            command_params_1=command_params_1,
            command_params_2=command_params_2,
            torus_temp=torus_temp,
            torus_luminosity=torus_luminosity,
            torus_frac=torus_frac,
            verbose=False,
            eic=eic,
        )


def load_from_backend(folder, flat=False):
    """
    Args:
        folder (str): The name of the folder containing the backend file.
        flat (bool, optional): If True, flatten the chain and log probability arrays.
            Defaults to False.

    Returns:
        tuple: A tuple containing the chain and log probability arrays.
    """
    reader = emcee.backends.HDFBackend(
        BASE_PATH + RESULTS_FOLDER + "/" + folder + "/backend.h5"
    )
    return reader.get_chain(flat=flat), reader.get_log_prob(flat=flat)
