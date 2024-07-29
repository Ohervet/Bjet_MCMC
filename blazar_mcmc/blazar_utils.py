#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file name: blazar_utils.py

Program Purpose: Implement necessary functions for mcmc (reading in configs
    and data, probability functions, etc.)

Note: All file paths are relative to Bjet_MCMC

Parameters are always listed in the following order:
[delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R]
The additional parameters for EIC are bb_temp, l_nuc, tau, blob_dist, in that order.
All parameters are the logarithm of the true value except for delta, n1, and n2
---------------------------------------------------
delta       doppler factor                  linear
K           particle density [cm^-3]        log
n1          first index                     linear
n2          second index                    linear
gamma_min   low-energy cutoff               log
gamma_max   high-energy cutoff              log
gamma_break energy break                    log
B           magnetic field strength [G]     log
R           blob radius [cm]                log
---------------------------------------------------
*Additional params for EIC*
bb_temp     Black body temp of disk [K]     log
l_nuc       Nucleus luminosity [ergs/s]     log
tau         Frac of luminosity scattered    log
blob_dist   Distance of blob [cm]           log
---------------------------------------------------

TODO: for eic
Contains:
TODO
FUNCTIONS:

read_configs
read_data
get_random_parameters
random_defaults
random_eic_from_std
min_max_parameters
----------------
log_prior
chi_squared_from_model
chi_squared
log_prob_from_model
log_probability


read_configs(config_file="mcmc_config.txt", config_string=None, verbose=False)
    Reads configurations from a configuration file or string.
    Returns a dictionary of configurations

read_data(data_file, cols=(0, 1, 4), use_E=True, verbose=False)
    Reads frequency or energy data and flux data from a file.
    Returns a tuple of np arrays of floats: v_data, vFv_data, err_data

get_random_parameters(param_min_vals=None, param_max_vals=None, alpha2_limits=None, redshift=None, tau_var=None,
                          use_variability=True, eic=False)
    Get a set of valid parameters.
    Returns an array of length # of dims with random params.

random_defaults(walkers, param_min_vals=None, param_max_vals=None, alpha2_limits=None, redshift=None, tau_var=None,
                    use_variability=True, eic=False)
    Given the minimum and maximum values for the parameters and the # of walkers,
    generate random default positions.
    Returns a np array of shape (walkers, dim)

random_eic_from_std(std_values, walkers, param_min_vals=None, param_max_vals=None, redshift=None,
                        tau_var=None, use_variability=True)
    Given a current state of a chain for a non-EIC run (with 9 free parameters),
    fill in the last 4 parameters with random defaults.

min_max_parameters(alpha2_limits=None, eic=False)
    Get the default minimum and maximum values.
    Returns the minima and maxima as a tuple of floats.

log_prior(params, param_min_vals=None, param_max_vals=None, redshift=None, tau_var=None,
              use_variability=True, alpha2_limits=None, eic=False)
              Given a set of parameters and the boundaries, determine if the parameters are
    "valid." Returns -np.inf if a parameter is out of range or the variability constraint is
    not satisfied.

chi_squared_from_model(model_results, v_data, vFv_data, err_data)
    Given model results (a tuple of 2+ np arrays with v and vFv as the first 2)
    and the data, return the chi squared value

chi_squared(params, v_data, vFv_data, err_data, name_stem=None, theta=None, redshift=None, min_freq=None,
                max_freq=None, torus_temp=None, torus_luminosity=None, torus_frac=None, data_folder=None,
                executable=None,
                command_params_full=None, command_params_1=None, command_params_2=None, parameter_file=None,
                prev_files=False, use_param_file=True, verbose=False, eic=False)
    Given parameters and data, create a model and then return the chi squared
    value for it.

log_prob_from_model(model_results, v_data, vFv_data, err_data)
    Given model results (a tuple of 2+ np arrays with v and vFv as the first 2)
    and the data, return the value for the probability function
    (-0.5 * chi squared)

log_probability(params, v_data, vFv_data, err_data, name_stem=None, param_min_vals=None,
                    param_max_vals=None,
                    theta=None, redshift=None, tau_var=None, use_variability=True, min_freq=None, max_freq=None,
                    torus_temp=None, torus_luminosity=None, torus_frac=None, data_folder=None, executable=None,
                    command_params_full=None, command_params_1=None, command_params_2=None, unique_name=False,
                    parameter_file=None,
                    prev_files=False, use_param_file=True, verbose=False, eic=False)
    Given parameters and data, make a model and return the value for the
    probability function (-0.5 * chi squared)
"""
import glob
import os
import random
import re

import numpy as np
from scipy import interpolate
from astropy.io import ascii
import blazar_model
from blazar_properties import *

#ev Hz conversion
h = 4.135667662E-15
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

    This can also be used to parse a string of dictionary values of configs.

    Args:
        config_file (optional): str
            Absolute path to file with configurations; default is None, using the 
            relative path to "mcmc_config.txt"
        config_string (optional): str
            This is a string of the form {'key': value, ...}. If a config_string
            is given, it will be parsed instead of reading from the file. This is
            useful when the values from a previous configuration dictionary are
            read from a file.
        verbose (optional): bool
            Controls whether values are shown; defaults to False

    Returns:
        (dictionary, tuple: (float, float))
        configurations, (lower alpha 2 limit, higher alpha 2 limit)

    Dictionary format:
    Keys are strings, value format is mixed
    Key (str)               Value
    "description"           (str) description--optional
    "folder_label"          (str) prefix for results folder name--optional
    "eic"                   (bool) eic or not
    "data_file"             (str) relative path to data
    "n_steps"               (int) number of steps to run
    "n_walkers"             (int) number of walkers
    "discard"               (int) how many steps are discarded before analysis
    "parallel"              (bool) whether parallel processing will be used
    "cores"                 (int) NOT PRESENT IF PARALLEL IS FALSE; # of cores
                                if absent and parallel TRUE,
                                (# cores on machine - 1) used
    "use_variability"       (bool) if parameters should be constrained based on variability constraint
    "tau_variability"       (float) time in hours for variability constraint
    "redshift"              (float) value for the redshift
    "custom_alpha2_limits"  (bool) whether custom alpha2 limits are used
    "alpha2_limits"         ([float, float]) alpha2 limits (set to default if none given)
    All model parameters to be frozen (float)
    """

    default_alpha2_limits = [1.5, 7.5]
    attributes = ["eic", "data_file", "n_steps", "n_walkers", "discard",
                  "parallel", "use_variability", "redshift", "custom_alpha2_limits"]
    optional_attributes = ["cores", "tau_variability", "description", "folder_label"]
    ssc_parameters = ["delta", "K","n1","n2", "gamma_min","gamma_max", "gamma_break","B","R"]
    eic_parameters = ["bb_temp", "l_nuc", "tau", "blob_dist"]
    configurations = {}  # dictionary of parameters
    if config_file is None:
        CONFIG_PATH = FOLDER_PATH+"mcmc_config.txt"
    else:
        CONFIG_PATH = config_file
    if config_string is None:
        # read configurations
        with open(CONFIG_PATH, 'r') as file:
            if verbose:
                print("Reading configuration options from:", config_file)
            for line in file:
                elements = re.split('=|# ', line)
                if elements[0].strip() in attributes or elements[0].strip() in optional_attributes \
                    or elements[0].strip() in ssc_parameters or elements[0].strip() in eic_parameters:  # determine if line with an attribute
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
        configurations["parallel"] = (configurations["parallel"] == "True" or configurations["parallel"] == "true")
        configurations["use_variability"] = (
                configurations["use_variability"] == "True" or configurations["use_variability"] == "true")
        configurations["redshift"] = float(configurations["redshift"])
    
        if "eic" in configurations:  # will default to false
            configurations["eic"] = (configurations["eic"] == "True" or configurations["eic"] == "true")
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
        configurations["custom_alpha2_limits"] = configurations["custom_alpha2_limits"].split('=')
        configurations["custom_alpha2_limits"] = configurations["custom_alpha2_limits"][-1].split(',')
        if configurations["custom_alpha2_limits"][0] == "True" or configurations["custom_alpha2_limits"] == "true":
            if len(configurations["custom_alpha2_limits"]) < 3:
                raise Exception("custom_alpha2_limits True but insufficient alpha2 limits provided!")
            alpha2_limits = [float(configurations["custom_alpha2_limits"][1]),
                             float(configurations["custom_alpha2_limits"][2])]
            alpha2_limits.sort()
            configurations["custom_alpha2_limits"] = True
            print("alpha2 limits sets at",alpha2_limits)
        else:
            configurations["custom_alpha2_limits"] = False
            alpha2_limits = default_alpha2_limits
        configurations["alpha2_limits"] = alpha2_limits
        
        #set up  fixed parameters
        if configurations["eic"]:
            all_parameters = ssc_parameters+eic_parameters
        else:
            all_parameters = ssc_parameters
        configurations["fixed_params"] = [-np.inf]*len(all_parameters)
        for i in range(len(all_parameters)):
            if configurations[all_parameters[i]] != 'null':
                configurations["fixed_params"][i] = float(configurations[all_parameters[i]])
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
    - Frequency or energy data in column 0, vFv data in column 1, and vFv error
        data in column 4
    - Space-separated (not changeable)
    - First row is a header (data not read from first row)
    - The first (#0) column has energy data in eV (which is then
        converted into frequency)

    Args:
        data_file: str
            relative path of the file w/ the data
        cols (optional): tuple of 3 ints
            columns with v or E data (v_data * h), vFv_data, and err_data; default
            is (0, 1, 4)
        use_E (optional): bool
            specifies if the first (#0) column has v data or E data; default
            is True
        verbose (optional): bool
            specifies if data are displayed; default is False
        instrument (optional): bool
            return the instrument use for each data point in addidion to a better data format for plotting the SED; default is False

    Returns:
        tuple of 3 1D np arrays of floats
        v_data, vFv_data, err_data
    """

    
    table = ascii.read(FOLDER_PATH + data_file, format='csv',data_start = 1,delimiter='\s')

                                             	 

    # get v data (may need to convert from E)
    if use_E:
        h = 4.135667662E-15  # eV * s
        v_data = table['!E(eV)']/h
        v_low = table['delta_E(-)']/h
        v_high = table['delta_E(+)']/h
    else:
        v_data = table['!E(eV)']
        v_low = table['delta_E(-)']
        v_high = table['delta_E(+)']
    vFv_data = table['F(ergcm-2s-1)']
    err_data_down = table['delta_F(-)']
    err_data_up = table['delta_F(+)']
    err_data = np.array([err_data_down, err_data_up])


    # # remove values where the error is 0 (upper limits)
    # if instrument == False:
    #     indices_to_remove = np.where(err_data == 0)[0]
    #     err_data_filtered = np.array([np.delete(err_data, indices_to_remove), np.delete(err_data_up, indices_to_remove)])
    #     v_data_filtered = np.delete(v_data, indices_to_remove)
    #     #nubin_data_filtered = np.array([np.delete(v_low, indices_to_remove), np.delete(v_high, indices_to_remove)])
    #     vFv_data_filtered = np.delete(vFv_data, indices_to_remove)

    #     if verbose:
    #         print("Data from", data_file)
    #         print("v_data_filtered:", v_data_filtered)
    #         print("vFv_data_filtered:", vFv_data_filtered)
    #         print("err_data_filtered:", err_data_filtered)
        
    if instrument:
        instrument_data = table['instrument']
        nubin_data = np.array([v_low,v_high])
        #instrument_data_filtered = np.delete(instrument_data, indices_to_remove)
        if verbose:
            print("Data from", data_file)
            print("v_data:", v_data)
            print("vFv_data:", vFv_data)
            print("err_data:", err_data)
            print("instrument_data:", instrument_data)
            print("nubin_data:", nubin_data)
        return v_data, vFv_data, err_data, instrument_data, nubin_data
    else:
        #return v_data_filtered, vFv_data_filtered, err_data_filtered
        return v_data, vFv_data, err_data



def get_random_parameters(param_min_vals=None, param_max_vals=None, alpha2_limits=None, redshift=None, tau_var=None,
                          use_variability=True, eic=False, fixed_params=None):
    """
    Get a random array of parameters in the parameter space.
    Args:
        param_min_vals (optional): np array of NUM_DIM floats
            minimum values for the params in the normal order; default is None, and
            the values are set to the defaults in blazar_utils.min_max_parameters
        param_max_vals (optional): np array of NUM_DIM floats
            maximum values for the params in the normal order; default is None, and
            the values are set to the defaults in blazar_utils.min_max_parameters
        alpha2_limits (optional): list/tuple of 2 floats
            values for alpha2 limits; default is None, and the alpha2 limits will be
            set to the defaults
        redshift (optional): float
            redshift value; default is None, so the log_prior function will use
            the default, which is 0.143 (the value for J1010)
        tau_var (optional): float
            time in hours for tau variability; default is None, so the log_prior
            function will use the default, which is 24 hours
        use_variability (optional): bool
            whether variability should be taken into account; default is True
        eic: bool
            states whether the run is eic or std; default is false (std)
    Returns:
        a numpy array of floats of length NUM_DIM

        random parameters within the min/max bounds. They will be valid
            (gamma min < gamma break < gamma max, etc.)
    """
    dim = len(param_min_vals)
    if param_min_vals is None or param_max_vals is None:
        default_min_max = min_max_parameters(alpha2_limits=alpha2_limits, eic=eic, fixed_params = fixed_params)
    if param_min_vals is None:
        param_min_vals = default_min_max[0]
    if param_max_vals is None:
        param_max_vals = default_min_max[1]

    # ensures valid parameters
    parameter_size = param_max_vals - param_min_vals
    parameters = param_min_vals + parameter_size * np.random.rand(dim)
            
    while not np.isfinite(log_prior(parameters, param_min_vals, param_max_vals, redshift=redshift, tau_var=tau_var,
                                    use_variability=use_variability, fixed_params=fixed_params)):
        parameters = param_min_vals + parameter_size * np.random.rand(dim)
    return parameters

 #need to implement gaussian ball of random parameters around "standard blazar values"
 # mu gauss = best params for J1010 for example
 #sig gauss = parameter space /10 (100?)


def random_defaults(walkers, param_min_vals=None, param_max_vals=None, alpha2_limits=None, redshift=None, tau_var=None,
                    use_variability=True, eic=False, fixed_params=None):
    """
    Get the values used for the initial values for the MCMC. The defaults are
    random values in the acceptable range that satisfy the log_prior criteria.

    Args:
        walkers: int
            number of walkers (specifies how many defaults to generate)
        param_min_vals (optional): 1D np array of NUM_DIM floats
            minimum values (in the standard order)
        param_max_vals (optional): 1D np array of NUM_DIM floats
            maximum values (in the standard order)
        alpha2_limits (optional): list/tuple of 2 floats
            values for alpha2 limits; default is None, and the alpha2 limits will be
            set to the default values
        redshift (optional): float
            redshift value; default is None, so the log_prior function will use
            the default, which is 0.143 (the value for J1010)
        tau_var (optional): float
            time in hours for tau variability; default is None, so the log_prior
            function will use the default, which is 24 hours
        use_variability (optional): bool
            whether variability should be taken into account; default is True
        eic: bool
            states whether the run is eic or std; default is false (std)
    Returns:
        2D np array of floats with NUM_DIM rows (the NUM_DIM parameters) and <walkers> rows

        default values for all NUM_DIM parameters for each walker
    """
    return np.array(
        [get_random_parameters(param_min_vals=param_min_vals, param_max_vals=param_max_vals,
                               alpha2_limits=alpha2_limits, redshift=redshift,
                               tau_var=tau_var,
                               use_variability=use_variability, eic=eic, fixed_params=fixed_params) for _ in range(walkers)])


def random_eic_from_std(std_values, walkers, param_min_vals=None, param_max_vals=None, redshift=None,
                        tau_var=None,
                        use_variability=True):
    """
    Given a current state of a chain for a non-EIC run (with 9 free parameters),
    fill in the last 4 parameters with random defaults.
    This is used when an EIC run that starts with non-EIC values approximated
    is desired.

    Args:
        std_values: 2D np array with # of walkers rows and 9 columns
            the current state for the std defaults
        walkers: int
            number of walkers (specifies how many defaults to generate)
        param_min_vals (optional): 1D np array of NUM_DIM floats
            minimum values (in the standard order)
        param_max_vals (optional): 1D np array of NUM_DIM floats
            maximum values (in the standard order)
        redshift (optional): float
            redshift value; default is None, so the log_prior function will use
            the default, which is 0.143 (the value for J1010)
        tau_var (optional): float
            time in hours for tau variability; default is None, so the log_prior
            function will use the default, which is 24 hours
        use_variability (optional): bool
            whether variability should be taken into account; default is True

    Returns:
        2D np array of floats with NUM_DIM rows (the NUM_DIM parameters) and <walkers> rows

        default values for all NUM_DIM parameters for each walker
    """
    if len(std_values[0]) != 9:
        raise ValueError(
            "Given defaults from a std run must have 9 values. " + str(len(std_values[0])) + " values found.")
    if len(param_min_vals) != 13:
        raise ValueError(
            "The minima and maxima must be for EIC with 13 values. " + str(len(param_min_vals)) + " values found.")
    if len(std_values) != walkers:
        raise ValueError("Number of walkers is " + str(walkers) + " but " + str(len(std_values)) + " defaults given")
    defaults = []
    for params in std_values:
        defaults.append(np.array(list(params) + list(
            get_random_parameters(param_min_vals=param_min_vals, param_max_vals=param_max_vals, redshift=redshift,
                                  tau_var=tau_var,
                                  use_variability=use_variability, eic=True, fixed_params=fixed_params))[9:]))
    return np.array(defaults)


def min_max_parameters(alpha2_limits=None, eic=False, fixed_params=None):
    """
    Get the default minimum and maximum values for all the  parameters.
    Arguments:
        alpha2_limits (optional): list/tuple of 2 floats
            values for alpha2 limits; default is None, and the alpha2 limits will be
            set to the defaults
        eic: bool
            states whether the run is eic or std; default is false (std)
    Returns:
        a tuple of two np arrays of NUM_DIM floats
        (param_min_vals, param_max_vals)

        both are in the standard order:
            [delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R]
        with eic:
            ["delta", "K", "n_1", "n_2", "gamma_min", "gamma_max", "gamma_break", "B", "R",
             "bb_temp", "l_nuc", "tau", "blob_dist"]
    """
    if alpha2_limits is None or len(alpha2_limits) != 2:
        alpha2_limits = (1.5, 7.5)
    param_min_vals = [1., 0., 1., float(alpha2_limits[0]), 0., 3., 2., -4., 14.]
    param_max_vals = [100, 8., 5., float(alpha2_limits[1]), 5., 8., 7.0, 0., 19.]
    if eic:
        extra_min = [3.5, 40.0, -5.0, 15]
        extra_max = [6.0, 50.0, 0.0, 21.0]
        param_min_vals = param_min_vals + extra_min
        param_max_vals = param_max_vals + extra_max
        
    #remove any frozen parameter from the parameter list
    if fixed_params:
        fixed_params2 = fixed_params.copy()
        i = 0
        while i < len(fixed_params2):
          if fixed_params2[i] != -np.inf:       
            del param_min_vals[i]
            del param_max_vals[i]
            del fixed_params2[i]
          else:
            i+=1
    return np.array(param_min_vals), np.array(param_max_vals)


# probability functions ------------------------------------------------------------------------------------------------
def log_prior(params, param_min_vals=None, param_max_vals=None, redshift=None, tau_var=None,
              use_variability=True, alpha2_limits=None, eic=False, fixed_params = None):
    """
    Using a uniform prior distribution. Return whether input params are valid.
    list parameters with eic: ["delta", "K", "n_1", "n_2", "gamma_min", "gamma_max", "gamma_break", "B", "R",
     "bb_temp", "l_nuc", "tau", "blob_dist"]
    Args:
        params: 1D np array of NUM_DIM floats
            Current position
        param_min_vals (optional): 1D np array of NUM_DIM floats
            Min param values
        param_max_vals (optional): 1D np array of NUM_DIM floats
            Max param values
        redshift (optional): float
            redshift value; default is None, so the log_prior function will use
            the default, which is 0.143 (the value for J1010)
        tau_var (optional): float
            time in hours for tau variability; default is None, so the log_prior
            function will use the default, which is 24 hours
        use_variability (optional): bool
            whether variability should be taken into account; default is True
        alpha2_limits (optional): tuple/list of 2 floats
            alpha2_limits (specify only if custom)
        eic: bool
            states whether the run is eic or std; default is false (std)
    Returns:
        float

        0 if all parameters are valid: n1 > n2, g_min < g_max,
        all parameters are in range and otherwise -np.inf
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
            
    
    if n1 > n2 or gamma_min > gamma_max or gamma_break < gamma_min or gamma_break > gamma_max:
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
        c = 2.997924 * 1.0e+10
        R = np.power(10, R)
        if tau_var < (1 + redshift) / c * R / delta:
            return -np.inf

    if eic:
        #the intrinsic jet half opening angle cannot be more than 5deg. (e.g. Hervet 2016)
        opening_angle = np.arctan(np.power(10, params[8])/np.power(10, params[12])) *180/pi
        if opening_angle > 5:
            return -np.inf
    return 0.0


def chi_squared_UL_to_err(P, func_nu_data, vFv_UL_data):
    """
    Method to include upper limits in Chi2 calculation.
    This will return the equivalent 1sigma error to be included in a standard Chi2 formula
    in a way that the resulting Chi2 will be equal to -2*ln(P).
    P being the UL likelihood (above and below)

    Parameters
    ----------
    P : float
        likelihood value associated with a data point.
        For ULs P is usually set either at 0.05 or 0.95 (probability above and below UL)
    func_nu_data : float
        model flux value a nu = nu_data
    vFv_UL_data : float
        measure flux UL

    Returns
    -------
    float: standard deviation to be included in a Chi2 formula

    """
    Chi2 = -2*np.log(P)
    return (func_nu_data-vFv_UL_data)/np.sqrt(Chi2)


def chi_squared_from_model(model_results, v_data, vFv_data, err_data):
    """
    Take the results of a model and computes the chi squared value.
    Consider asymmetric error bars in the datast

    Args:
        model_results: tuple of 4 numpy 1D arrays of floats
            Model results (logv, logvFv, v, vFv)
        v_data: 1D numpy arrays of floats
            Data values (NOT LOG)
        vFv_data: 1D numpy arrays of floats
            Data values (NOT LOG)
        err_data: 2D numpy arrays of floats [err_data_down, err_data_up]
            Data values (NOT LOG)

    Returns:
        float: the chi squared value
    """

    logv_all = model_results[0]
    logvFv_all = model_results[1]
    
    # calculate chi-squared by plugging in the v_data_values into the interpolation
    # from v_all to vFv_all
    func = interpolate.interp1d(logv_all, logvFv_all, fill_value='extrapolate')
    func_nu_data = np.power(10, func(np.log10(v_data)))
    
    diff = func_nu_data - vFv_data
    #check if model is above or below data flux points, True if above
    sign = diff > 0
    
    #transform err_data to consider ULs (P=95% when below, P=5% when above)
    if_UL = err_data[0] == 0
    err_data[0] += if_UL * chi_squared_UL_to_err(0.95, func_nu_data, vFv_data)
    err_data[1] += if_UL * chi_squared_UL_to_err(0.05, func_nu_data, vFv_data)
    
    return np.sum((diff / (sign*err_data[1] + ~sign*err_data[0]))**2.)




def chi_squared(params, v_data, vFv_data, err_data, name_stem=None, theta=None, redshift=None, min_freq=None,
                max_freq=None, torus_temp=None, torus_luminosity=None, torus_frac=None, data_folder=None,
                executable=None,
                command_params_full=None, command_params_1=None, command_params_2=None, parameter_file=None,
                prev_files=False, use_param_file=True, verbose=False, eic=False):
    """
    Get the chi squared value for a given set of parameters and data.

    The purpose of the option to supply command params is to speed up computation--
    the values for the first portion of the parameters and the second portion will
    be the same for every run in the MCMC, and it speeds up the code significantly
    to pass them as arguments instead of re-creating them every time.
    The entire list of parameters (including the ones that are changed in the MCMC)
    can be supplied with command_params_full. In this case, params are ignored.
    Alternatively, only the constant param values are provided with command_params_1
    and command_params_2.

    Args:
        params: 1D np array of NUM_DIM floats
            Current parameters
        v_data: 1D np array of floats
            Data values for v (NOT LOG)
        vFv_data: 1D np array of floats
            Data values for vFv (NOT LOG)
        err_data: 2D np array of floats
            Data values for error
        name_stem (optional): str
            Name stem for make_model. Default is none; will be then set to default.
        theta (optional): float
            Angle from the line of sight. Default is none, and it will be set to
            the default value of 0.57.
        redshift (optional): float
            redshift value; default is None, so the log_prior function will use
            the default, which is 0.143 (the value for J1010)
        min_freq (optional): float
            Minimum frequency for the SED model. Default is none, where it will be
            set to the default value of 5.0e+7 in blazar_model.process_model.
        max_freq (optional): float
            Maximum frequency for the SED model. Default is none, where it will be
            set to the default value of 1.0e+29 in blazar_model.process_model.
        torus_temp (optional): float
            Value for the torus temperature. Default is none, and it will be set
            to the default of 2.0e+4
        torus_luminosity (optional): float
            Value for the torus luminosity. Default is none, and it will be set
            to the default of 5.5e+20
        torus_frac (optional): float
            Value for the fraction of the torus luminosity reprocessed isotropically.
            Default is none, and it will be set to the default of 9.0e-5
        data_folder (optional): str
            Relative path to the folder where data will be saved to
        executable (optional): str
            Where the bjet executable is located
        command_params_full (optional): numpy array of floats
            Full set of parameters to pass to the bjet exec--see the README
            for information. Should be length 22, 23, 28, or 29.
        command_params_1 (optional): numpy array of floats
            The settings and transformation parameters:
            [prev files flag, data folder, model type, redshift, hubble constant,
            theta]
        command_params_2 (optional): numpy array of floats
            The constant and numerical parameters:
            [length of the emitting region, absorption by EBL, blob distance,
            # of spectral points, min freq, max freq, file name prefix]
        parameter_file (optional): str
            Relative path of the parameter file (the file where parameters
            are written to when modeling, **will be overwritten**) Default is
            <PARAMETER_FOLDER>/params.txt (PARAMETER_FOLDER is in home directory)
        prev_files (optional): bool
            Whether bjet should create _prev files after each run; default is False
        use_param_file (optional): bool
            Whether bjet should be called with a parameter file or with command
            line args; default is False
        verbose (optional): bool
            Whether information on the model should be displayed; default is
            False
        eic: bool
            states whether the run is eic or std; default is false (std)

    Returns:
        float
        The chi_squared value
    """
    model_results = blazar_model.make_model(params, name_stem=name_stem, theta=theta, redshift=redshift,
                                            min_freq=min_freq, max_freq=max_freq, torus_temp=torus_temp,
                                            torus_luminosity=torus_luminosity, torus_frac=torus_frac,
                                            data_folder=data_folder, executable=executable,
                                            command_params_full=command_params_full, command_params_1=command_params_1,
                                            command_params_2=command_params_2, parameter_file=parameter_file,
                                            prev_files=prev_files, use_param_file=use_param_file, verbose=verbose,
                                            eic=eic)
    
    return chi_squared_from_model(model_results, v_data, vFv_data, err_data)


def log_prob_from_model(model_results, v_data, vFv_data, err_data):
    """
    This is the log_prob for the modeling (bigger = better fit).
    It returns -0.5 * the chi squared value for the v and vFv values from the
    given model.

    Args:
        model_results: tuple of 4 numpy 1D arrays of floats
            model results (logv, logvFv, v, vFv)
        v_data: 1D numpy arrays of floats
            Data for frequency
        vFv_data: 1D numpy arrays of floats
            Data for energy flux
        err_data: 2D numpy arrays of floats
            Data for vFv error

    Returns:
        float
        -0.5 * the chi squared value
    """
    return -0.5 * chi_squared_from_model(model_results, v_data, vFv_data, err_data)


def log_probability(params, v_data, vFv_data, err_data, name_stem=None, param_min_vals=None,
                    param_max_vals=None,
                    theta=None, redshift=None, tau_var=None, use_variability=True, min_freq=None, max_freq=None,
                    torus_temp=None, torus_luminosity=None, torus_frac=None, data_folder=None, executable=None,
                    command_params_full=None, command_params_1=None, command_params_2=None, unique_name=False,
                    parameter_file=None,
                    prev_files=False, use_param_file=True, verbose=False, eic=False, fixed_params=None):
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

    Args:
        params: 1D np array of NUM_DIM floats
            Current parameters
        v_data: 1D np array of floats
            Data values for v (NOT LOG)
        vFv_data: 1D np array of floats
            Data values for vFv (NOT LOG)
        err_data: 1D np array of floats
            Data values for error
        name_stem (optional): str
            Name stem for make_model. Default is none; will be then set to default.
        param_min_vals (optional): np array of NUM_DIM floats
            minimum values for the params in the normal order; default is None, and
            the values are set to the defaults in blazar_utils.min_max_parameters
        param_max_vals (optional): np array of NUM_DIM floats
            maximum values for the params in the normal order; default is None, and
            the values are set to the defaults in blazar_utils.min_max_parameters
        theta (optional): float
            Angle from the line of sight. Default is none, and it will be set to
            the default value of 0.57.
        redshift (optional): float
            redshift value; default is None, so the log_prior function will use
            the default, which is 0.143 (the value for J1010)
        tau_var (optional): float
            time in hours for tau variability; default is None, so the log_prior
            function will use the default, which is 24 hours
        use_variability (optional): bool
            whether variability should be taken into account; default is True
        min_freq (optional): float
            Minimum frequency for the SED model. Default is none, where it will be
            set to the default value of 5.0e+7 in blazar_model.process_model.
        max_freq (optional): float
            Maximum frequency for the SED model. Default is none, where it will be
            set to the default value of 1.0e+29 in blazar_model.process_model.
        torus_temp (optional): float
            Value for the torus temperature. Default is none, and it will be set
            to the default of 2.0e+4
        torus_luminosity (optional): float
            Value for the torus luminosity. Default is none, and it will be set
            to the default of 5.5e+20
        torus_frac (optional): float
            Value for the fraction of the torus luminosity reprocessed isotropically.
            Default is none, and it will be set to the default of 9.0e-5
        data_folder (optional): str
            Relative path to the folder where data will be saved to
        executable (optional): str
            Where the bjet executable is located
        command_params_full (optional): numpy array of floats
            Full set of parameters to pass to the bjet exec--see the README
            for information. Should be length 22, 23, 28, or 29.
        command_params_1 (optional): numpy array of floats
            The settings and transformation parameters:
            [prev files flag, data folder, model type, redshift, hubble constant,
            theta]
        command_params_2 (optional): numpy array of floats
            The constant and numerical parameters:
            [length of the emitting region, absorption by EBL, blob distance,
            # of spectral points, min freq, max freq, file name prefix]
        unique_name (optional): string
            Specifies if the name stem should be created to be unique. This uses
            a random number, creating a very low risk of conflicts; default is
            False.
        parameter_file (optional): string
            Name of parameter file. This will be created from name_stem if not
            provided.
        prev_files (optional): bool
            Whether bjet should create _prev files after each run; default is False
        use_param_file (optional): bool
            Whether bjet should be called with a parameter file or with command
            line args; default is False
        verbose (optional): bool
            Whether information on the model should be displayed; default is
            False
        eic: bool
            states whether the run is eic or std; default is false (std)

    Returns:
        float
        -0.5 * chi squared value
    """
    #all frozen parameters need to be reimplemented for the likelihood computation
    if fixed_params:
        for i in range(len(fixed_params)):
            if fixed_params[i]!= -np.inf:
                params = np.insert(params, i, fixed_params[i])
                param_min_vals = np.insert(param_min_vals, i, fixed_params[i])
                param_max_vals = np.insert(param_max_vals, i, fixed_params[i])
        
    
    if not np.isfinite(log_prior(params, param_min_vals, param_max_vals, redshift=redshift, tau_var=tau_var,
                                 use_variability=use_variability)):
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
    result = -0.5 * chi_squared(params, v_data, vFv_data, err_data, name_stem=name_stem, theta=theta,
                                redshift=redshift, min_freq=min_freq,
                                max_freq=max_freq, torus_temp=torus_temp, torus_luminosity=torus_luminosity,
                                torus_frac=torus_frac, data_folder=data_folder, executable=executable,
                                command_params_full=command_params_full,
                                command_params_1=command_params_1, command_params_2=command_params_2,
                                parameter_file=parameter_file, prev_files=prev_files, use_param_file=use_param_file,
                                verbose=verbose, eic=eic)
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