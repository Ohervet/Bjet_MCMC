#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file name: blazar_utils.py

Program Purpose: Implement necessary functions for mcmc (reading in configs
    and data, probability functions, etc.)

Note: All file paths are relative to Bjet_MCMC

Parameters are always listed in the following order:
``[delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R]``

The additional parameters for EIC are ``bb_temp, l_nuc, tau, blob_dist``, in that order.
All parameters are the logarithm of the true value except for delta, n1, and n2.

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
"""
# TODO: for eic
# Contains:
# TODO
# FUNCTIONS:
#
# read_configs
# read_data
# get_random_parameters
# random_defaults
# random_eic_from_std
# min_max_parameters
# ----------------
# log_prior
# chi_squared_from_model
# chi_squared
# log_prob_from_model
# log_probability
#
#
# read_configs(config_file="mcmc_config.txt", config_string=None, verbose=False)
#     Reads configurations from a configuration file or string.
#     Returns a dictionary of configurations
#
# read_data(data_file, cols=(0, 1, 4), use_E=True, verbose=False)
#     Reads frequency or energy data and flux data from a file.
#     Returns a tuple of np arrays of floats: v_data, vFv_data, err_data
#
# get_random_parameters(param_min_vals=None, param_max_vals=None, alpha2_limits=None, redshift=None, tau_var=None,
#                           use_variability=True, eic=False)
#     Get a set of valid parameters.
#     Returns an array of length # of dims with random params.
#
# random_defaults(walkers, param_min_vals=None, param_max_vals=None, alpha2_limits=None, redshift=None, tau_var=None,
#                     use_variability=True, eic=False)
#     Given the minimum and maximum values for the parameters and the # of walkers,
#     generate random default positions.
#     Returns a np array of shape (walkers, dim)
#
# random_eic_from_std(std_values, walkers, param_min_vals=None, param_max_vals=None, redshift=None,
#                         tau_var=None, use_variability=True)
#     Given a current state of a chain for a non-EIC run (with 9 free parameters),
#     fill in the last 4 parameters with random defaults.
#
# min_max_parameters(alpha2_limits=None, eic=False)
#     Get the default minimum and maximum values.
#     Returns the minima and maxima as a tuple of floats.
#
# log_prior(params, param_min_vals=None, param_max_vals=None, redshift=None, tau_var=None,
#               use_variability=True, alpha2_limits=None, eic=False)
#               Given a set of parameters and the boundaries, determine if the parameters are
#     "valid." Returns -np.inf if a parameter is out of range or the variability constraint is
#     not satisfied.
#
# chi_squared_from_model(model_results, v_data, vFv_data, err_data)
#     Given model results (a tuple of 2+ np arrays with v and vFv as the first 2)
#     and the data, return the chi squared value
#
# chi_squared(params, v_data, vFv_data, err_data, name_stem=None, theta=None, redshift=None, min_freq=None,
#                 max_freq=None, torus_temp=None, torus_luminosity=None, torus_frac=None, data_folder=None,
#                 executable=None,
#                 command_params_full=None, command_params_1=None, command_params_2=None, parameter_file=None,
#                 prev_files=False, use_param_file=True, verbose=False, eic=False)
#     Given parameters and data, create a model and then return the chi squared
#     value for it.
#
# log_prob_from_model(model_results, v_data, vFv_data, err_data)
#     Given model results (a tuple of 2+ np arrays with v and vFv as the first 2)
#     and the data, return the value for the probability function
#     (-0.5 * chi squared)
#
# log_probability(params, v_data, vFv_data, err_data, name_stem=None, param_min_vals=None,
#                     param_max_vals=None,
#                     theta=None, redshift=None, tau_var=None, use_variability=True, min_freq=None, max_freq=None,
#                     torus_temp=None, torus_luminosity=None, torus_frac=None, data_folder=None, executable=None,
#                     command_params_full=None, command_params_1=None, command_params_2=None, unique_name=False,
#                     parameter_file=None,
#                     prev_files=False, use_param_file=True, verbose=False, eic=False)
#     Given parameters and data, make a model and return the value for the
#     probability function (-0.5 * chi squared)

import glob
import os
import random
import re

import numpy as np
from scipy import interpolate
from astropy.io import ascii
from bjet_mcmc import blazar_model
from bjet_mcmc.blazar_properties import *

__all__ = [
    "chi_squared",
    "chi_squared_from_model",
    "chi_squared_Limit_to_err",
    "e_to_v",
    "get_random_parameters",
    "log_prior",
    "log_prob_from_model",
    "log_probability",
    "min_max_parameters",
    "random_defaults",
    "random_eic_from_std",
    "read_configs",
    "read_data",
    "v_to_e",
]

# ev Hz conversion
h = 4.135667662e-15


def v_to_e(val):
    return val * h


def e_to_v(val):
    return val / h


def read_configs(config_file=None, config_string=None, verbose=False):
    """
    Given a relative path to a file, read the mcmc configs.
    If custom alpha 2 limits are not specified, the defaults of 1.5 and 7.5 will
    be used.

    See README for format of the configs file.

    This can also be used to parse a string of dictionary values of configs

    :param config_file: Absolute path to file with configurations; default is None, using the relative path to "mcmc_config.txt"
    :type config_file: str or None
    :param config_string: This is a string of the form :code:`{'key': value, ...}`. If a config_string is given, it will be parsed instead of reading from the file. This is useful when the values from a previous configuration dictionary are  read from a file.
    :type config_string: str or None
    :param verbose: Controls whether values are shown; defaults to False
    :type verbose: bool
    :return: The configurations as a dictionary.
    :rtype: dict
    """

    default_alpha2_limits = [1.5, 7.5]
    attributes = [
        "eic",
        "data_file",
        "p0",
        "n_steps",
        "n_walkers",
        "discard",
        "parallel",
        "use_variability",
        "redshift",
        "custom_alpha2_limits",
    ]
    optional_attributes = ["cores", "tau_variability", "description", "folder_label"]
    ssc_parameters = [
        "delta",
        "K",
        "n1",
        "n2",
        "gamma_min",
        "gamma_max",
        "gamma_break",
        "B",
        "R",
    ]
    eic_parameters = ["bb_temp", "l_nuc", "tau", "blob_dist"]
    configurations = {}  # dictionary of parameters
    if config_file is None:
        CONFIG_PATH = FOLDER_PATH + "mcmc_config.txt"
    else:
        CONFIG_PATH = config_file
    if config_string is None:
        # read configurations
        with open(CONFIG_PATH, "r") as file:
            if verbose:
                print("Reading configuration options from:", config_file)
            for line in file:
                elements = re.split("=|# ", line)
                if (
                    elements[0].strip() in attributes
                    or elements[0].strip() in optional_attributes
                    or elements[0].strip() in ssc_parameters
                    or elements[0].strip() in eic_parameters
                ):  # determine if line with an attribute
                    configurations[elements[0].strip()] = elements[1].strip()
    else:
        inf = np.inf
        configurations = eval(config_string)

    # check if all configs present
    for att in attributes:
        if att not in configurations:
            raise Exception("No " + att + " provided!")

    if config_string is None:
        # change int params to ints, bools to bools
        configurations["n_steps"] = int(configurations["n_steps"])
        configurations["n_walkers"] = int(configurations["n_walkers"])
        configurations["discard"] = int(configurations["discard"])
        configurations["parallel"] = (
            configurations["parallel"] == "True" or configurations["parallel"] == "true"
        )
        configurations["use_variability"] = (
            configurations["use_variability"] == "True"
            or configurations["use_variability"] == "true"
        )
        configurations["redshift"] = float(configurations["redshift"])

        if "eic" in configurations:  # will default to false
            configurations["eic"] = (
                configurations["eic"] == "True" or configurations["eic"] == "true"
            )
        else:
            configurations["eic"] = False

        if configurations["use_variability"]:
            configurations["tau_variability"] = float(configurations["tau_variability"])
            if "tau_variability" not in configurations:
                raise Exception("No tau_variability provided!")
        else:
            configurations["tau_variability"] = None

        if "cores" in configurations:
            configurations["cores"] = int(configurations["cores"])

        # set alpha2 limits
        configurations["custom_alpha2_limits"] = configurations[
            "custom_alpha2_limits"
        ].split("=")
        configurations["custom_alpha2_limits"] = configurations["custom_alpha2_limits"][
            -1
        ].split(",")
        if (
            configurations["custom_alpha2_limits"][0] == "True"
            or configurations["custom_alpha2_limits"] == "true"
        ):
            if len(configurations["custom_alpha2_limits"]) < 3:
                raise Exception(
                    "custom_alpha2_limits True but insufficient alpha2 limits provided!"
                )
            alpha2_limits = [
                float(configurations["custom_alpha2_limits"][1]),
                float(configurations["custom_alpha2_limits"][2]),
            ]
            alpha2_limits.sort()
            configurations["custom_alpha2_limits"] = True
            print("alpha2 limits sets at", alpha2_limits)
        else:
            configurations["custom_alpha2_limits"] = False
            alpha2_limits = default_alpha2_limits
        configurations["alpha2_limits"] = alpha2_limits

        # set up  fixed parameters
        if configurations["eic"]:
            all_parameters = ssc_parameters + eic_parameters
        else:
            all_parameters = ssc_parameters
        configurations["fixed_params"] = [-np.inf] * len(all_parameters)
        for i in range(len(all_parameters)):
            if configurations[all_parameters[i]] != "null":
                configurations["fixed_params"][i] = float(
                    configurations[all_parameters[i]]
                )
            configurations.pop(all_parameters[i])

    if verbose:
        # show parameters
        print("Configurations:")
        for att in attributes:
            print("  ", att, "=", configurations[att])
        for att in optional_attributes:
            if att in configurations:
                print("  ", att, "=", configurations[att])

    return configurations


def read_data(data_file, cols=(0, 1, 4), use_E=True, verbose=False, instrument=False):
    """
    Read the data for the SED (energy data, vFv, flux). By default, the
    data file satisfies the following:
    - Frequency or energy data in column 0, vFv data in column 1, and vFv error data in column 4
    - Space-separated (not changeable)
    - First row is a header (data not read from first row)
    - The first (#0) column has energy data in eV (which is then converted into frequency)

    :param data_file: relative path of the file w/ the data
    :type data_file: str
    :param cols: columns with v or E data (v_data * h), vFv_data, and err_data; default is (0, 1, 4)
    :type cols: tuple[int]
    :param use_E: specifies if the first (#0) column has v data or E data; default is True
    :type use_E: bool
    :param verbose: Flag indicating whether to print additional information during execution. Default is False.
    :type verbose: bool
    :param instrument: return the instrument use for each data point in addidion to a better data format for plotting the SED; default is False
    :type instrument: bool
    :return: tuple of 3 1D np arrays of floats v_data, vFv_data, err_data
    :rtype: tuple

    """

    table = ascii.read(
        FOLDER_PATH + data_file, format="csv", data_start=1, delimiter="\s"
    )

    # get v data (may need to convert from E)
    if use_E:
        h = 4.135667662e-15  # eV * s
        v_data = table["!E(eV)"] / h
        v_low = table["delta_E(-)"] / h
        v_high = table["delta_E(+)"] / h
    else:
        v_data = table["!E(eV)"]
        v_low = table["delta_E(-)"]
        v_high = table["delta_E(+)"]
    vFv_data = table["F(ergcm-2s-1)"]
    err_data_down = table["delta_F(-)"]
    err_data_up = table["delta_F(+)"]
    err_data = np.array([err_data_down, err_data_up])

    if instrument:
        instrument_data = table["instrument"]
        nubin_data = np.array([v_low, v_high])
        if verbose:
            print("Data from", data_file)
            print("v_data:", v_data)
            print("vFv_data:", vFv_data)
            print("err_data:", err_data)
            print("instrument_data:", instrument_data)
            print("nubin_data:", nubin_data)
        return v_data, vFv_data, err_data, instrument_data, nubin_data
    else:
        return v_data, vFv_data, err_data


def get_random_parameters(
    param_min_vals=None,
    param_max_vals=None,
    alpha2_limits=None,
    redshift=None,
    tau_var=None,
    use_variability=True,
    eic=False,
    fixed_params=None,
):
    """
    Get a random array of parameters in the parameter space.

    :param param_min_vals: minimum values for the params in the normal order; default is None, and
            the values are set to the defaults in blazar_utils.min_max_parameters If not provided, default_min_max is used.
    :type param_min_vals: numpy.array or None
    :param param_max_vals: maximum values for the params in the normal order; default is None, and
            the values are set to the defaults in blazar_utils.min_max_parameters. If not provided, default_min_max is used.
    :type param_max_vals: numpy.array or None
    :param alpha2_limits: values for alpha2 limits; default is None, and the alpha2 limits will be set to the defaults
    :type alpha2_limits: numpy.array or None
    :param redshift: redshift value; default is None, so the log_prior function will use the default, which is 0.143 (the value for J1010)
    :type redshift: float or None
    :param tau_var: time in hours for tau variability; default is None, so the log_prior function will use the default, which is 24 hours
    :type tau_var: float or None
    :param use_variability: whether variability should be taken into account; default is True
    :type use_variability: bool
    :param eic: states whether the run is eic or std; default is false (std)
    :type eic: bool
    :param fixed_params: Fixed parameters. If not provided, default_min_max is used.
    :type fixed_params: dict or None
    :return: Randomly generated parameters.
    :rtype: numpy.array

    .. note:: random parameters within the min/max bounds. They will be valid (gamma min < gamma break < gamma max, etc.)

    The function first checks if param_min_vals and param_max_vals are provided. If not, it uses default_min_max.
    If param_min_vals is not provided, it assigns default_min_max[0] to param_min_vals.
    If param_max_vals is not provided, it assigns default_min_max[1] to param_max_vals.
    Then it calculates the parameter_size and adds param_min_vals with parameter_size * random numbers between 0 and 1.
    It checks if the generated parameters are finite using log_prior function. If not, it generates new parameters until they are valid.
    Finally, it returns the valid parameters.

    """
    dim = len(param_min_vals)
    if param_min_vals is None or param_max_vals is None:
        default_min_max = min_max_parameters(
            alpha2_limits=alpha2_limits, eic=eic, fixed_params=fixed_params
        )
    if param_min_vals is None:
        param_min_vals = default_min_max[0]
    if param_max_vals is None:
        param_max_vals = default_min_max[1]

    # ensures valid parameters
    parameter_size = param_max_vals - param_min_vals
    parameters = param_min_vals + parameter_size * np.random.rand(dim)

    while not np.isfinite(
        log_prior(
            parameters,
            param_min_vals,
            param_max_vals,
            redshift=redshift,
            tau_var=tau_var,
            use_variability=use_variability,
            fixed_params=fixed_params,
            eic=eic
        )
    ):
        parameters = param_min_vals + parameter_size * np.random.rand(dim)
    return parameters


# need to implement gaussian ball of random parameters around "standard blazar values"
# mu gauss = best params for J1010 for example
# sig gauss = parameter space /10 (100?)

#need docstring
def TwoSided_Multivariate_Normal(means, cov_matrix_down, cov_matrix_up):
    parameter_samples_down = np.random.multivariate_normal(means, cov_matrix_down)
    parameter_samples_up = np.random.multivariate_normal(means, cov_matrix_up)
    parameters = parameter_samples_down.copy()
    diff = means - parameter_samples_down
    sign = diff > 0
    # be sure that parameter_samples_up are aways more than the mean
    # by default from our method parameter_samples_down will always be below the mean
    parameter_samples_up = means + np.abs(parameter_samples_up - means)
    # when a projection is below the mean, use parameter_samples_down, when above use parameter_samples_up
    parameters = (
        parameter_samples_down * sign + parameter_samples_up * ~sign
    )
    return parameters

#function below is very preliminary
def get_gaussian_parameters(
    param_file="p0_file.txt",
    param_min_vals=None,
    param_max_vals=None,
    redshift=None,
    tau_var=None,
    use_variability=True,
    eic=False,
    fixed_params=None,
):
    
    #read from a parameter files that contains 1 sigma boundaries in linear space
    p0 = ascii.read(param_file)
    par_name= p0["Parameter"]
    logpar_name = ["K", "gamma_min", "gamma_max", "gamma_break", "B", "R", "bb_temp", "l_nuc", "tau", "blob_dist"]

    #set some parameters into log10 space
    for i in range(len(par_name)):
        if par_name[i] in logpar_name:
            p0["Best_Value"][i] = np.log10(p0["Best_Value"][i])
            p0["low_bound"][i] = np.log10(p0["low_bound"][i])
            p0["up_bound"][i] = np.log10(p0["up_bound"][i])

    # consider asymmetric errors in parameters (two-sided multivariate normal method)
    means = p0["Best_Value"]
    params_error_down =  means - p0["low_bound"]
    params_error_up   =  p0["up_bound"] - means
    cov_matrix_down = np.diag(params_error_down) ** 2
    cov_matrix_up = np.diag(params_error_up) ** 2
    
    parameters = TwoSided_Multivariate_Normal(means, cov_matrix_down, cov_matrix_up)
    
    
    # ensures valid parameters
    if param_min_vals is None or param_max_vals is None:
        default_min_max = min_max_parameters(
            eic=eic, fixed_params=fixed_params
        )
    if param_min_vals is None:
        param_min_vals = default_min_max[0]
    if param_max_vals is None:
        param_max_vals = default_min_max[1]

    while not np.isfinite(
        log_prior(
            parameters,
            param_min_vals,
            param_max_vals,
            redshift=redshift,
            tau_var=tau_var,
            use_variability=use_variability,
            fixed_params=fixed_params,
            eic=eic
        )
    ):
        parameters = TwoSided_Multivariate_Normal(means, cov_matrix_down, cov_matrix_up)
        
    return parameters


def random_defaults(
    walkers,
    param_min_vals=None,
    param_max_vals=None,
    alpha2_limits=None,
    redshift=None,
    tau_var=None,
    use_variability=True,
    eic=False,
    fixed_params=None,
):
    """
    Get the values used for the initial values for the MCMC. The defaults are random values in the acceptable range that satisfy the log_prior criteria.

    :param walkers: number of walkers (specifies how many defaults to generate)
    :type walkers: int

    :param param_min_vals: Minimum values for the parameters. If None, default values will be used. (in the standard order)
    :type param_min_vals: numpy.array or None

    :param param_max_vals: Maximum values for the parameters. If None, default values will be used. (in the standard order)
    :type param_max_vals: numpy.array or None

    :param alpha2_limits: Limits for the alpha2 parameter. If None, default values will be used.
    :type alpha2_limits: tuple or None

    :param redshift: Redshift value. Default is None, so the log_prior function will use the default, which is 0.143 (the value for J1010)
    :type redshift: float or None

    :param tau_var: Variability value for tau. default is None, so the log_prior function will use the default, which is 24 hours
    :type tau_var: float or None

    :param use_variability: Flag to indicate whether variability should be used. Default is True.
    :type use_variability: bool

    :param eic: Flag to indicate whether EIC (Extended Inverted Chirality) should be used. Default is False.
    :type eic: bool

    :param fixed_params: Values for fixed parameters. If None, default values will be used.
    :type fixed_params: numpy.array or None

    :return: A numpy array with the random default parameter values for each walker. 2D np array of floats with NUM_DIM rows (the NUM_DIM parameters) and <walkers> rows. default values for all NUM_DIM parameters for each walker
    :rtype: numpy.array
    """
    return np.array(
        [
            get_random_parameters(
                param_min_vals=param_min_vals,
                param_max_vals=param_max_vals,
                alpha2_limits=alpha2_limits,
                redshift=redshift,
                tau_var=tau_var,
                use_variability=use_variability,
                eic=eic,
                fixed_params=fixed_params,
            )
            for _ in range(walkers)
        ]
    )

#docstring & inputs needs to be checked
#param files needs to be included
def random_gaussian(
    walkers,
    param_min_vals=None,
    param_max_vals=None,
    alpha2_limits=None,
    redshift=None,
    tau_var=None,
    use_variability=True,
    eic=False,
    fixed_params=None,
):
    """
    Get the values used for the initial values for the MCMC. The defaults are random values in the acceptable range that satisfy the log_prior criteria.

    :param walkers: number of walkers (specifies how many defaults to generate)
    :type walkers: int

    :param param_min_vals: Minimum values for the parameters. If None, default values will be used. (in the standard order)
    :type param_min_vals: numpy.array or None

    :param param_max_vals: Maximum values for the parameters. If None, default values will be used. (in the standard order)
    :type param_max_vals: numpy.array or None

    :param alpha2_limits: Limits for the alpha2 parameter. If None, default values will be used.
    :type alpha2_limits: tuple or None

    :param redshift: Redshift value. Default is None, so the log_prior function will use the default, which is 0.143 (the value for J1010)
    :type redshift: float or None

    :param tau_var: Variability value for tau. default is None, so the log_prior function will use the default, which is 24 hours
    :type tau_var: float or None

    :param use_variability: Flag to indicate whether variability should be used. Default is True.
    :type use_variability: bool

    :param eic: Flag to indicate whether EIC (Extended Inverted Chirality) should be used. Default is False.
    :type eic: bool

    :param fixed_params: Values for fixed parameters. If None, default values will be used.
    :type fixed_params: numpy.array or None

    :return: A numpy array with the random default parameter values for each walker. 2D np array of floats with NUM_DIM rows (the NUM_DIM parameters) and <walkers> rows. default values for all NUM_DIM parameters for each walker
    :rtype: numpy.array
    """
    return np.array(
        [
            get_gaussian_parameters(
                param_min_vals=param_min_vals,
                param_max_vals=param_max_vals,
                redshift=redshift,
                tau_var=tau_var,
                use_variability=use_variability,
                eic=eic,
                fixed_params=fixed_params,
            )
            for _ in range(walkers)
        ]
    )


def random_eic_from_std(
    std_values,
    walkers,
    param_min_vals=None,
    param_max_vals=None,
    redshift=None,
    tau_var=None,
    use_variability=True,
):
    """
    Given a current state of a chain for a non-EIC run (with 9 free parameters),
    fill in the last 4 parameters with random defaults.
    This is used when an EIC run that starts with non-EIC values approximated
    is desired.

    :param std_values: 2D np array with # of walkers rows and 9 columns the current state for the std defaults
    :type std_values: list of lists
    :param walkers: number of walkers (specifies how many defaults to generate)
    :type walkers: int
    :param param_min_vals: 1D np array of NUM_DIM floats minimum values (in the standard order)
    :type param_min_vals: list, optional
    :param param_max_vals:  1D np array of NUM_DIM floats maximum values (in the standard order)
    :type param_max_vals: list, optional
    :param redshift: redshift value; default is None, so the log_prior function will use the default, which is 0.143 (the value for J1010)
    :type redshift: float, optional
    :param tau_var: time in hours for tau variability; default is None, so the log_prior function will use the default, which is 24 hours.
    :type tau_var: float, optional
    :param use_variability: Determines whether to use variability in generating random parameters. Defaults to True.
    :type use_variability: bool, optional
    :return: 2D np array of floats with NUM_DIM rows (the NUM_DIM parameters) and <walkers> rows default values for all NUM_DIM parameters for each walker
    :rtype: numpy.ndarray
    :raises ValueError: If the length of std_values is not equal to the number of walkers or if the length of std_values[0] is not 9.
    :raises ValueError: If the length of param_min_vals is not 13.
    """
    if len(std_values[0]) != 9:
        raise ValueError(
            "Given defaults from a std run must have 9 values. "
            + str(len(std_values[0]))
            + " values found."
        )
    if len(param_min_vals) != 13:
        raise ValueError(
            "The minima and maxima must be for EIC with 13 values. "
            + str(len(param_min_vals))
            + " values found."
        )
    if len(std_values) != walkers:
        raise ValueError(
            "Number of walkers is "
            + str(walkers)
            + " but "
            + str(len(std_values))
            + " defaults given"
        )
    defaults = []
    for params in std_values:
        defaults.append(
            np.array(
                list(params)
                + list(
                    get_random_parameters(
                        param_min_vals=param_min_vals,
                        param_max_vals=param_max_vals,
                        redshift=redshift,
                        tau_var=tau_var,
                        use_variability=use_variability,
                        eic=True,
                        fixed_params=fixed_params,
                    )
                )[9:]
            )
        )
    return np.array(defaults)


def min_max_parameters(alpha2_limits=None, eic=False, fixed_params=None):
    """
    Get the default minimum and maximum values for all the  parameters.

    :param alpha2_limits: A tuple of two values representing the lower and upper limits for the alpha2 parameter. Defaults to (1.5, 7.5).
    :type alpha2_limits: tuple
    :param eic: A boolean value indicating whether to include extra parameters. Defaults to False.
    :type eic: bool
    :param fixed_params: A list of fixed parameter values to be removed from the parameter list. Defaults to None.
    :type fixed_params: list or None
    :return: A tuple containing two arrays representing the minimum and maximum parameter values.
    :rtype: tuple

    .. note::
        both returns are in the standard order:
            - :code:`[delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R]`
            - with eic:  :code:`["delta", "K", "n_1", "n_2", "gamma_min", "gamma_max", "gamma_break", "B", "R", "bb_temp", "l_nuc", "tau", "blob_dist"]`
    """
    if alpha2_limits is None or len(alpha2_limits) != 2:
        alpha2_limits = (1.5, 7.5)
    param_min_vals = [1.0, 0.0, 1.0, float(alpha2_limits[0]), 0.0, 3.0, 2.0, -4.0, 14.0]
    param_max_vals = [100, 8.0, 5.0, float(alpha2_limits[1]), 5.0, 8.0, 7.0, 0.0, 19.0]
    if eic:
        extra_min = [3.5, 40.0, -5.0, 15]
        extra_max = [6.0, 50.0, 0.0, 21.0]
        param_min_vals = param_min_vals + extra_min
        param_max_vals = param_max_vals + extra_max

    # remove any frozen parameter from the parameter list
    if fixed_params:
        fixed_params2 = fixed_params.copy()
        i = 0
        while i < len(fixed_params2):
            if fixed_params2[i] != -np.inf:
                del param_min_vals[i]
                del param_max_vals[i]
                del fixed_params2[i]
            else:
                i += 1
    return np.array(param_min_vals), np.array(param_max_vals)


# probability functions ------------------------------------------------------------------------------------------------
def log_prior(
    params,
    param_min_vals=None,
    param_max_vals=None,
    redshift=None,
    tau_var=None,
    use_variability=True,
    alpha2_limits=None,
    eic=False,
    fixed_params=None,
):
    """
    Using a uniform prior distribution. Return whether input params are valid.
    list parameters with eic: ["delta", "K", "n_1", "n_2", "gamma_min", "gamma_max", "gamma_break", "B", "R", "bb_temp", "l_nuc", "tau", "blob_dist"]

    :param params: The parameter values.
    :type params: numpy.ndarray
    :param param_min_vals: The minimum allowed values for each parameter. If None, the function will use default values based on the specified constraints.
    :type param_min_vals: numpy.ndarray
    :param param_max_vals: The maximum allowed values for each parameter. If None, the function will use default values based on the specified constraints.
    :type param_max_vals: numpy.ndarray
    :param redshift: redshift value; default is None, so the log_prior function will use the default, which is 0.143 (the value for J1010)
    :type redshift: float
    :param tau_var: time in hours for tau variability; default is None, so the log_prior function will use the default, which is 24 hours
    :type tau_var: float
    :param use_variability: whether variability should be taken into account; default is True
    :type use_variability: bool
    :param alpha2_limits: The limits for the alpha2 parameter. If None, the function will use default values based on the specified constraints.
    :type alpha2_limits: tuple
    :param eic: A flag indicating whether to consider EIC constraints. Default is False.
    :type eic: bool
    :param fixed_params: The fixed parameter values. If provided, the function will replace the corresponding parameters in the 'params' array with the fixed values.
    :type fixed_params: numpy.ndarray

    :return: 0 if all parameters are valid: n1 > n2, g_min < g_max, all parameters are in range and otherwise -np.inf
    :rtype: float
    """

    if fixed_params:
        # reimplent fixed params
        for i in range(len(fixed_params)):
            if fixed_params[i] != -np.inf:
                params = np.insert(params, i, fixed_params[i])
                param_min_vals = np.insert(param_min_vals, i, fixed_params[i])
                param_max_vals = np.insert(param_max_vals, i, fixed_params[i])

    if param_min_vals is None or param_max_vals is None:
        minima, maxima = min_max_parameters(alpha2_limits=alpha2_limits, eic=eic)
        if param_min_vals is None:
            param_min_vals = minima
        if param_max_vals is None:
            param_max_vals = maxima

    delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R, *other_params = params

    if (
        n1 > n2
        or gamma_min > gamma_max
        or gamma_break < gamma_min
        or gamma_break > gamma_max
    ):
        return -np.inf
    # check if between min and max
    for i in range(len(params)):
        # gamma break is solely constrained by gamma_min and gamma_max
        if i != 6 and (param_min_vals[i] > params[i] or param_max_vals[i] < params[i]):
            return -np.inf
    if tau_var is None:
        tau_var = 1e10

    # tau variability constraint
    if redshift is None:
        print("here in log prior")
        pass
        # redshift = 0.143
    if use_variability:
        tau_var = tau_var * 60 * 60  # convert to seconds
        c = 2.997924 * 1.0e10
        R = np.power(10, R)
        if tau_var < (1 + redshift) / c * R / delta:
            return -np.inf

    if eic:
        # the intrinsic jet half opening angle cannot be more than 5deg. (e.g. Hervet 2016)
        opening_angle = (
            np.arctan(np.power(10, params[8]) / np.power(10, params[12])) * 180 / np.pi
        )
        if opening_angle > 5:
            return -np.inf
    return 0.0


def chi_squared_Limit_to_err(P, func_nu_data, vFv_Limit_data):
    """
    Method to include upper and lower limits in Chi2 calculation.
    This will return the equivalent 1sigma error to be included in a standard Chi2 formula
    in a way that the resulting Chi2 will be equal to -2*ln(P).
    P being the UL/LL likelihood (above and below)

    :param P: likelihood value associated with a data point. For ULs/LLs P is usually set either at 0.05 or 0.95 (probability above and below UL/LL)
    :type P: float
    :param func_nu_data: Model flux value a nu = nu_data
    :type func_nu_data: numpy.array
    :param vFv_Limit_data: measure flux UL/LL.
    :type vFv_Limit_data: numpy.array
    :return: The error limit calculated. standard deviation to be included in a Chi2 formula
    :rtype: numpy.array
    """
    Chi2 = -2 * np.log(P)
    return np.abs(func_nu_data - vFv_Limit_data) / np.sqrt(Chi2)


def chi_squared_from_model(model_results, v_data, vFv_data, err_data):
    """
    Calculates the chi-squared value based on the given model_results, v_data, vFv_data, and err_data. Consider asymmetric error bars in the dataset

    :param model_results: The results of the model calculation. Should contain logv_all and logvFv_all. Model results (logv, logvFv, v, vFv)
    :type model_results: list
    :param v_data: The observed v values. Data values (NOT LOG)
    :type v_data: numpy array
    :param vFv_data: The observed vFv values. Data values (NOT LOG)
    :type vFv_data: numpy array
    :param err_data: The error values for vFv data. Should contain two values: the lower limit and the upper limit. [err_data_down, err_data_up] Data values (NOT LOG)
    :type err_data: numpy array
    :return: The calculated chi-squared value.
    :rtype: float
    """

    logv_all = model_results[0]
    logvFv_all = model_results[1]

    # calculate chi-squared by plugging in the v_data_values into the interpolation
    # from v_all to vFv_all
    func = interpolate.interp1d(logv_all, logvFv_all, fill_value="extrapolate")
    func_nu_data = np.power(10, func(np.log10(v_data)))
    
    # code can break in the unlikely, but possible case when model = data for ULs
    # add an tiny deviation to the model to avoid this case
    spot_on =  func_nu_data - vFv_data == 0
    func_nu_data += spot_on * np.full(len(func_nu_data),func_nu_data/1.0e4)

    diff = func_nu_data - vFv_data
    # check if model is above or below data flux points, True if above
    sign = diff > 0

    # transform err_data to consider ULs (P=95% when below, P=5% when above)
    if_UL = err_data[0] == 0
    err_data[0] += if_UL * chi_squared_Limit_to_err(0.95, func_nu_data, vFv_data)
    err_data[1] += if_UL * chi_squared_Limit_to_err(0.05, func_nu_data, vFv_data)
    # transform err_data to consider LLs (P=95% when above, P=5% when below)
    if_LL = err_data[0] == -1
    err_data[0] += if_LL * (chi_squared_Limit_to_err(0.05, func_nu_data, vFv_data)+1)
    err_data[1] += if_LL * chi_squared_Limit_to_err(0.95, func_nu_data, vFv_data)
    
    return np.sum((diff / (sign * err_data[1] + ~sign * err_data[0])) ** 2.0)


def chi_squared(
    params,
    v_data,
    vFv_data,
    err_data,
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
    parameter_file=None,
    prev_files=False,
    use_param_file=True,
    verbose=False,
    eic=False,
):
    """
    Given parameters and data, create a model and then return the chi squared value for it.

    Get the chi squared value for a given set of parameters and data.

    The purpose of the option to supply command params is to speed up computation--
    the values for the first portion of the parameters and the second portion will
    be the same for every run in the MCMC, and it speeds up the code significantly
    to pass them as arguments instead of re-creating them every time.
    The entire list of parameters (including the ones that are changed in the MCMC)
    can be supplied with command_params_full. In this case, params are ignored.
    Alternatively, only the constant param values are provided with command_params_1
    and command_params_2.

    :param params: The parameters of the model.
    :type params: Any
    :param v_data: The observed values of the dependent variable (v). Data values for v (NOT LOG)
    :type v_data: List or ndarray
    :param vFv_data: The observed values of the dependent variable (vFv). Data values for vFv (NOT LOG)
    :type vFv_data: List or ndarray
    :param err_data: The observed values of the dependent variable errors.
    :type err_data: List or ndarray
    :param name_stem: The name stem for the output and intermediate files. (optional) Default is none; will be then set to default.
    :type name_stem: str
    :param theta: The angle of the jet with respect to the line of sight in degrees. (optional) Default is none, and it will be set to the default value of 0.57.
    :type theta: float
    :param redshift: The redshift of the object. (optional). Default is None, so the log_prior function will use the default, which is 0.143 (the value for J1010)
    :type redshift: float
    :param min_freq: The minimum frequency for the model calculations. Minimum frequency for the SED model. Default is none, where it will be set to the default value of 5.0e+7 in blazar_model.process_model.(optional)
    :type min_freq: float
    :param max_freq: The maximum frequency for the model calculations. Maximum frequency for the SED model. Default is none, where it will be set to the default value of 1.0e+29 in blazar_model.process_model.(optional)
    :type max_freq: float
    :param torus_temp: The temperature of the torus. (optional) Default is none, and it will be set to the default of 2.0e+4
    :type torus_temp: float
    :param torus_luminosity: The luminosity of the torus. (optional) Default is none, and it will be set to the default of 5.5e+20
    :type torus_luminosity: float
    :param torus_frac: The fraction of the total luminosity from the torus. (optional) Value for the fraction of the torus luminosity reprocessed isotropically.  Default is none, and it will be set to the default of 9.0e-5
    :type torus_frac: float
    :param data_folder: Relative path to the folder where data will be saved to. (optional)
    :type data_folder: str
    :param executable: Where the bjet executable is located. (optional)
    :type executable: str
    :param command_params_full: Full set of parameters to pass to the bjet exec--see the README for information. Should be length 22, 23, 28, or 29.
    :type command_params_full: str
    :param command_params_1: The settings and transformation parameters: [prev files flag, data folder, model type, redshift, hubble constant, theta]
    :type command_params_1: str
    :param command_params_2: The constant and numerical parameters: [length of the emitting region, absorption by EBL, blob distance, # of spectral points, min freq, max freq, file name prefix]
    :type command_params_2: str
    :param parameter_file: Relative path of the parameter file (the file where parameters are written to when modeling, **will be overwritten**) Default is <PARAMETER_FOLDER>/params.txt (PARAMETER_FOLDER is in home directory)
    :type parameter_file: str
    :param prev_files: Whether bjet should create _prev files after each run; default is False
    :type prev_files: bool
    :param use_param_file: Whether bjet should be called with a parameter file or with command line args; default is False
    :type use_param_file: bool
    :param verbose: Whether information on the model should be displayed; default is False
    :type verbose: bool
    :param eic: states whether the run is eic or std; default is false (std)

    :type eic: bool
    :return: The chi-squared value for the model.
    :rtype: float
    """
    model_results = blazar_model.make_model(
        params,
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
        parameter_file=parameter_file,
        prev_files=prev_files,
        use_param_file=use_param_file,
        verbose=verbose,
        eic=eic,
    )

    return chi_squared_from_model(model_results, v_data, vFv_data, err_data)


def log_prob_from_model(model_results, v_data, vFv_data, err_data):
    """
    This is the log_prob for the modeling (bigger = better fit).
    It returns -0.5 * the chi squared value for the v and vFv values from the
    given model.

    :param model_results: The results of the model calculation. model results (logv, logvFv, v, vFv)
    :type model_results: Any
    :param v_data: Data for frequency
    :type v_data: Any
    :param vFv_data: Data for energy flux
    :type vFv_data: Any
    :param err_data: Data for vFv error
    :type err_data: Any
    :return: The log probability. -0.5 * the chi squared value
    :rtype: float

    """
    return -0.5 * chi_squared_from_model(model_results, v_data, vFv_data, err_data)


def log_probability(
    params,
    v_data,
    vFv_data,
    err_data,
    name_stem=None,
    param_min_vals=None,
    param_max_vals=None,
    theta=None,
    redshift=None,
    tau_var=None,
    use_variability=True,
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
    unique_name=False,
    parameter_file=None,
    prev_files=False,
    use_param_file=True,
    verbose=False,
    eic=False,
    fixed_params=None,
):
    """
    This is the log_prob for the modeling (bigger = better fit).
    It returns -0.5 * the chi squared value for the v and vFv values from
    the model with the given parameters.

    The purpose of the option to supply command params is to speed up computation--
    the values for the first portion of the parameters and the second portion will
    be the same for every run in the MCMC, and it speeds up the code significantly
    to pass them as arguments instead of re-creating them every time.
    The entire list of parameters (including the ones that are changed in the MCMC)
    can be supplied with command_params_full. In this case, params are ignored.
    Alternatively, only the constant param values are provided with command_params_1
    and command_params_2.

    :param params: The set of parameters.
    :type params: numpy array
    :param v_data: The observed frequency data. (NOT LOG)
    :type v_data: numpy array
    :param vFv_data: The observed frequency times flux data. (NOT LOG)
    :type vFv_data: numpy array
    :param err_data: The observed error data. (NOT LOG)
    :type err_data: numpy array
    :param name_stem: The name stem used in file naming. (optional)
    :type name_stem: str
    :param param_min_vals: The minimum values for each parameter. minimum values for the params in the normal order; default is None, and the values are set to the defaults in blazar_utils.min_max_parameters
    :type param_min_vals: numpy array
    :param param_max_vals: The maximum values for the params in the normal order; default is None, and the values are set to the defaults in blazar_utils.min_max_parameters
    :type param_max_vals: numpy array
    :param theta: The theta value. (optional). Angle from the line of sight. Default is none, and it will be set to the default value of 0.57.
    :type theta: float
    :param redshift: redshift value; default is None, so the log_prior function will use the default, which is 0.143 (the value for J1010)
    :type redshift: float
    :param tau_var: The tau_var value. (optional). Time in hours for tau variability; default is None, so the log_prior function will use the default, which is 24 hours
    :type tau_var: float
    :param use_variability: True to use variability, False otherwise. (default is True)
    :type use_variability: bool
    :param min_freq: The minimum frequency value. (optional) Default is none, where it will be set to the default value of 5.0e+7 in blazar_model.process_model.
    :type min_freq: float
    :param max_freq: The maximum frequency value. (optional) Default is none, where it will be set to the default value of 1.0e+29 in blazar_model.process_model.
    :type max_freq: float
    :param torus_temp: The torus temperature value. (optional) Default is none, and it will be set to the default of 2.0e+4
    :type torus_temp: float
    :param torus_luminosity: The torus luminosity value. (optional) Default is none, and it will be set to the default of 5.5e+20
    :type torus_luminosity: float
    :param torus_frac: The torus fraction value. (optional) Value for the fraction of the torus luminosity reprocessed isotropically. Default is none, and it will be set to the default of 9.0e-5
    :type torus_frac: float
    :param data_folder: The data folder path. (optional)
    :type data_folder: str
    :param executable: The executable file path. (optional)
    :type executable: str
    :param command_params_full: Full set of parameters to pass to the bjet exec--see the README for information. Should be length 22, 23, 28, or 29.
    :type command_params_full: list
    :param command_params_1: The settings and transformation parameters: [prev files flag, data folder, model type, redshift, hubble constant, theta]
    :type command_params_1: list
    :param command_params_2: The constant and numerical parameters: [length of the emitting region, absorption by EBL, blob distance, # of spectral points, min freq, max freq, file name prefix]
    :type command_params_2: list
    :param unique_name: Specifies if the name stem should be created to be unique. This uses a random number, creating a very low risk of conflicts; default is
    :type unique_name: bool
    :param parameter_file: Name of parameter file. This will be created from name_stem if not provided.
    :type parameter_file: str
    :param prev_files: Whether bjet should create _prev files after each run; default is False
    :type prev_files: bool
    :param use_param_file: Whether bjet should be called with a parameter file or with command line args; default is False
    :type use_param_file: bool
    :param verbose: True to enable verbose mode, False otherwise. (default is False)
    :type verbose: bool
    :param eic: True to enable EIC mode, False otherwise. (default is False)
    :type eic: bool
    :param fixed_params: The fixed parameters. (optional)
    :type fixed_params: numpy array
    :return: The log-likelihood value. -0.5 * chi squared value
    :rtype: float
    """
    # all frozen parameters need to be reimplemented for the likelihood computation
    if fixed_params:
        for i in range(len(fixed_params)):
            if fixed_params[i] != -np.inf:
                params = np.insert(params, i, fixed_params[i])
                param_min_vals = np.insert(param_min_vals, i, fixed_params[i])
                param_max_vals = np.insert(param_max_vals, i, fixed_params[i])

    if not np.isfinite(
        log_prior(
            params,
            param_min_vals,
            param_max_vals,
            redshift=redshift,
            tau_var=tau_var,
            use_variability=use_variability,
            eic=eic
        )
    ):
        return -np.inf
    # not infinite
    if name_stem is None:
        name_stem = NAME_STEM
    if unique_name:
        name_stem = name_stem + "_" + str(random.getrandbits(60))
    if use_param_file and parameter_file is None:
        parameter_file = PARAMETER_FOLDER + "/" + name_stem + ".txt"
    if command_params_2 is not None:
        command_params_2[-1] = name_stem  # change the file stem
    result = -0.5 * chi_squared(
        params,
        v_data,
        vFv_data,
        err_data,
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
        parameter_file=parameter_file,
        prev_files=prev_files,
        use_param_file=use_param_file,
        verbose=verbose,
        eic=eic,
    )
    if use_param_file:
        os.remove(parameter_file)

    for f in glob.glob(BASE_PATH + DATA_FOLDER + "/" + name_stem + "_*"):
        os.remove(f)
    if prev_files:
        for f in glob.glob(BASE_PATH + DATA_FOLDER + "/*_prev*.dat"):
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
    return result
