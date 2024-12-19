#!/usr/bin/env python3

import datetime
import multiprocessing
import shutil

import emcee

from bjet_mcmc import blazar_utils
from bjet_mcmc import blazar_report
from bjet_mcmc import blazar_properties
from bjet_core import bj_core

import glob
import numpy as np
#import pathlib
import subprocess
import os
import random
from scipy import interpolate
from ctypes import *

__all__ = [
    "RESULTS_FOLDER",
    "DATA_FOLDER",
    "EXECUTABLE",
    "BASE_PATH",
    "configs",
    "EIC",
    "DIM",
    "redshift",
    "use_variability",
    "tau",
    "alpha2_limits",
    "model_type",
    "settings_and_transformation",
    "constant_and_numerical",
    "param_min_vals",
    "param_max_vals",
    "param_min_vals",
    "param_max_vals",
    "make_model",
    "log_prior",
    "random_params",
    "log_prob",
    "mcmc",
]

# TODO Removed PROGRAM_NAME, so now results directory must go elsewhere (since that was the reference to it)
RESULTS_FOLDER = "local_results"
DATA_FOLDER = "sed_calculations"
EXECUTABLE = "bjet_core/bj_core"
BASE_PATH = blazar_properties.BASE_PATH

configs = blazar_utils.read_configs()
v_data, vFv_data, err_data = blazar_utils.read_data(configs["data_file"])

EIC = configs["eic"]
DIM = 13 if EIC else 9
if EIC:
    PARAM_IS_LOG = [
        False,
        True,
        False,
        False,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
    ]
else:
    PARAM_IS_LOG = [False, True, False, False, True, True, True, True, True]

redshift = configs["redshift"]
use_variability = configs["use_variability"]
tau = configs["tau_variability"]
alpha2_limits = configs["alpha2_limits"]

model_type = "1" if EIC else "0"
settings_and_transformation = [
    "3",
    BASE_PATH + DATA_FOLDER,
    model_type,
    str(redshift),
    "69.6",
    "0.57",
]
constant_and_numerical = ["0", "1", "9.0e+17", "99", "5.0e+7", "1e+29", "temp_stem"]

# mins maxes
param_min_vals = [1.0, 0.0, 1.0, float(alpha2_limits[0]), 0.0, 3.0, 2.0, -4.0, 14.0]
param_max_vals = [100.0, 8.0, 5.0, float(alpha2_limits[1]), 5.0, 8.0, 6.699, 0.0, 19.0]
if EIC:
    extra_min = [1.0, 20.0, -10.0, 10.0]
    extra_max = [6.0, 50.0, 0.0, 21.0]
    param_min_vals = param_min_vals + extra_min
    param_max_vals = param_max_vals + extra_max
param_min_vals = np.array(param_min_vals)
param_max_vals = np.array(param_max_vals)


# def make_model_prev(params, name_stem="run"):
#     """
#     This function builds a model based on the provided parameters and name stem. It converts the parameters to logarithmic values where necessary. The model is built by executing a command using subprocess.run() and the output is saved to a file. The logarithmic wavelength and flux values are then extracted from the file and converted to linear values. Additional components are added to the model based on the name stem. Lastly, the built model is interpolated and combined with the existing model.

#     :param params: The parameters for building the model.
#     :type params: list[float]
#     :param name_stem: The prefix for the output file names.
#     :type name_stem: str
#     :return: The logarithmic wavelength values, logarithmic flux values, linear wavelength values, and linear flux values.
#     :rtype: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]

#     A tuple containing the following elements:
#         - logv (float[]): The logarithm of the frequency values.
#         - logvFv (float[]): The logarithm of the flux values.
#         - v (float[]): The frequency values.
#         - vFv (float[]): The flux values.

#     .. note:: This implementation uses subprocess.popen and should be rewritten in the future with bj_core.py methods.
#     """
#     # convert to log
#     # TODO rewrite this entire method using bj_core methods and without building a command string.
#     params = params * 1.0
#     for i in range(len(params)):
#         if PARAM_IS_LOG[i]:
#             params[i] = np.power(10, params[i])
#     if not EIC:
#         command = (
#             settings_and_transformation
#             + [str(val) for val in params]
#             + constant_and_numerical[:-1]
#             + [name_stem]
#         )
#     else:
#         command = (
#             settings_and_transformation
#             + [str(val) for val in params[:9]]
#             + constant_and_numerical[:2]
#             + [str(params[-1])]
#         )
#         command = command + [
#             str(params[9]),
#             "2.0e+4",
#             str(params[10]),
#             str(params[11]),
#             "5.5e+20",
#             "9.0e-5",
#         ]  # EIC components
#         command = command + constant_and_numerical[3:-1] + [name_stem]
#     subprocess.run(
#         [BASE_PATH + EXECUTABLE, *command],
#         stderr=open(os.devnull, "wb"),
#         stdout=open(os.devnull, "wb"),
#     )

#     loaded_model = np.loadtxt(
#         BASE_PATH + DATA_FOLDER + "/" + name_stem + "_ss.dat", delimiter=" "
#     )
#     logv = loaded_model[:, 0]
#     logvFv = loaded_model[:, 2]
#     vFv = np.power(10, logvFv)

#     stems = ["cs"] if not EIC else ["cs", "ecs", "cs2", "nuc"]
#     for s in stems:
#         loaded_model = np.loadtxt(
#             BASE_PATH + DATA_FOLDER + "/" + name_stem + "_" + s + ".dat", delimiter=" "
#         )
#         model_logv = loaded_model[:, 0]
#         model_logvFv = loaded_model[:, 2]
#         current_logv, current_logvFv, current_vFv = logv, logvFv, vFv

#         new_lower = np.where(model_logv < logv[0])[0]
#         new_higher = np.where(model_logv > logv[-1])[0]
#         logv = np.concatenate(
#             (model_logv[new_lower], current_logv, model_logv[new_higher])
#         )
#         logvFv = np.concatenate(
#             (model_logvFv[new_lower], current_logvFv, model_logvFv[new_higher])
#         )
#         vFv = np.power(10, logvFv)

#         overlap_start = np.where(logv >= max(model_logv[0], current_logv[0]))[0][0]
#         overlap_end = np.where(logv <= min(model_logv[-1], current_logv[-1]))[0][-1]

#         interpolation = interpolate.interp1d(model_logv, np.power(10, model_logvFv))

#         new_vFv = np.concatenate(
#             (
#                 np.zeros(overlap_start),
#                 interpolation(logv[overlap_start : overlap_end + 1]),
#                 np.zeros(len(logv) - overlap_end - 1),
#             )
#         )

#         vFv = vFv + new_vFv
#         logvFv = np.log10(vFv)
#     return logv, logvFv, np.power(10, logv), vFv



# def make_model(params, name_stem="run"):
#     """
#     This function builds a model based on the provided parameters and name stem. 
#     It converts the parameters to logarithmic values where necessary. 
    
#     new:
#     The model is built by running the c++ function  bj_core.main_swig() that returns a single flat 1D array of all components
#     outputs (log(nu), log(nuFnu)).
#     Each cemission component is then retrived by slicing this array.
#     old:
#     executing a command using subprocess.run() and the output is saved to a file. 
#     The logarithmic wavelength and flux values are then extracted from the file and converted to linear values. 
    
#     Additional components are added to the model based on the name stem. Lastly, 
#     the built model is interpolated and combined with the existing model.

#     :param params: The parameters for building the model.
#     :type params: list[float]
#     :param name_stem: The prefix for the output file names.
#     :type name_stem: str
#     :return: The logarithmic wavelength values, logarithmic flux values, linear wavelength values, and linear flux values.
#     :rtype: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]

#     A tuple containing the following elements:
#         - logv (float[]): The logarithm of the frequency values.
#         - logvFv (float[]): The logarithm of the flux values.
#         - v (float[]): The frequency values.
#         - vFv (float[]): The flux values.

#     .. note:: This implementation uses subprocess.popen and should be rewritten in the future with bj_core.py methods.
#     """
#     # convert to log
#     # TODO rewrite this entire method using bj_core methods and without building a command string.
#     params = params * 1.0
#     for i in range(len(params)):
#         if PARAM_IS_LOG[i]:
#             params[i] = np.power(10, params[i])
#     if not EIC:
#         command = (
#             settings_and_transformation
#             + [str(val) for val in params]
#             + constant_and_numerical[:-1]
#             + [name_stem]
#         )
#     else:
#         command = (
#             settings_and_transformation
#             + [str(val) for val in params[:9]]
#             + constant_and_numerical[:2]
#             + [str(params[-1])]
#         )
#         command = command + [
#             str(params[9]),
#             "2.0e+4",
#             str(params[10]),
#             str(params[11]),
#             "5.5e+20",
#             "9.0e-5",
#         ]  # EIC components
#         command = command + constant_and_numerical[3:-1] + [name_stem]
#     subprocess.run(
#         [BASE_PATH + EXECUTABLE, *command],
#         stderr=open(os.devnull, "wb"),
#         stdout=open(os.devnull, "wb"),
#     )
#     #New method, without subprocess run:
#     #change into strings
#     list_params = list((map(str, params))) + constant_and_numerical[:-1]
#     aa = bj_core.main_swig(*list_params, 3, "/home/olivier/Bjet_MCMC_root/Bjet_MCMC/sed_calculations", name_stem)

#     loaded_model = np.loadtxt(
#         BASE_PATH + DATA_FOLDER + "/" + name_stem + "_ss.dat", delimiter=" "
#     )
#     logv = loaded_model[:, 0]
#     logvFv = loaded_model[:, 2]
#     vFv = np.power(10, logvFv)

#     stems = ["cs"] if not EIC else ["cs", "ecs", "cs2", "nuc"]
#     for s in stems:
#         loaded_model = np.loadtxt(
#             BASE_PATH + DATA_FOLDER + "/" + name_stem + "_" + s + ".dat", delimiter=" "
#         )

#         #New method
#         if s == "cs":
#             NU_DIM = 99
#             n_components = 2
#             full_size = int(NU_DIM*n_components)
            
#             p = (ctypes.c_double * full_size).from_address(int(aa))
#             p = list(p)
#             #print(p)
#             # #this work for 1D array
#             start = 0
#             stop = NU_DIM
#             SSC = [0]*2
#             for i in range(n_components):
#                 SSC[i] = p[start:stop+start]
#                 start += NU_DIM
#                 #print(SSC[i])
                
#             #remove indexes with 0 flux
#             while 0 in SSC[1]:
#                 tmp = SSC[1].index(0)
#                 del SSC[1][tmp]
#                 del SSC[0][tmp]
#             model_logv = SSC[0]
#             model_logvFv = SSC[1]
#             #end new method
            
#         else:
#             model_logv = loaded_model[:, 0]
#             model_logvFv = loaded_model[:, 2]
        
#         current_logv, current_logvFv, current_vFv = logv, logvFv, vFv

#         new_lower = np.where(model_logv < logv[0])[0]
#         new_higher = np.where(model_logv > logv[-1])[0]
#         logv = np.concatenate(
#             (model_logv[new_lower], current_logv, model_logv[new_higher])
#         )
#         logvFv = np.concatenate(
#             (model_logvFv[new_lower], current_logvFv, model_logvFv[new_higher])
#         )
#         vFv = np.power(10, logvFv)

#         overlap_start = np.where(logv >= max(model_logv[0], current_logv[0]))[0][0]
#         overlap_end = np.where(logv <= min(model_logv[-1], current_logv[-1]))[0][-1]

#         interpolation = interpolate.interp1d(model_logv, np.power(10, model_logvFv))

#         new_vFv = np.concatenate(
#             (
#                 np.zeros(overlap_start),
#                 interpolation(logv[overlap_start : overlap_end + 1]),
#                 np.zeros(len(logv) - overlap_end - 1),
#             )
#         )

#         vFv = vFv + new_vFv
#         logvFv = np.log10(vFv)
#     return logv, logvFv, np.power(10, logv), vFv


# def log_prior(params):
#     """
#     This function calculates the log prior probability for a given set of parameters.

#     Args:
#         params (tuple): A tuple containing the following parameters:
#             - delta (float): The value of delta.
#             - K (float): The value of K.
#             - n1 (float): The value of n1.
#             - n2 (float): The value of n2.
#             - gamma_min (float): The value of gamma_min.
#             - gamma_max (float): The value of gamma_max.
#             - gamma_break (float): The value of gamma_break.
#             - B (float): The value of B.
#             - R (float): The value of R.
#             - other_params (tuple): Optional additional parameters.

#     Returns:
#         float: The log prior value. If any of the parameters violate the constraints, it returns negative infinity (-inf). Otherwise, it returns 0.

#     :param params: The input parameters to calculate the log prior probability.
#     :type params: tuple

#     :return: The log prior probability.
#     :rtype: float
#     """
#     delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R, *other_params = params
#     if (
#         n1 > n2
#         or gamma_min > gamma_max
#         or gamma_break < gamma_min
#         or gamma_break > gamma_max
#     ):
#         return -np.inf
#     # testing if between min and max
#     for i in range(len(params)):
#         # gamma break is solely constrained by gamma_min and gamma_max
#         if i != 6 and (param_min_vals[i] > params[i] or param_max_vals[i] < params[i]):
#             return -np.inf
#     if use_variability:
#         tau_var = tau * 60 * 60  # to seconds
#         c = 2.997924 * 1.0e10
#         R = np.power(10, R)
#         if tau_var < (1 + redshift) / c * R / delta:
#             return -np.inf
#     return 0


def random_params():
    """
    Generates random parameters within given range and checks if they have finite log prior values and are within specified minimum and maximum values.

    Returns:
        numpy.ndarray: Array of random parameters.

    Raises:
        None.

    Examples:
        >>> random_params()
        array([0.5, 1.2, 0.8])

    :return: Randomly generated parameters within the specified range.
    :rtype: numpy.ndarray
    """
    parameter_size = param_max_vals - param_min_vals
    parameters = param_min_vals + parameter_size * np.random.rand(DIM)
    while not np.isfinite(blazar_utils.log_prior(parameters)):
        parameters = param_min_vals + parameter_size * np.random.rand(DIM)
    return parameters


# def log_prob(params):
#     """

#     This function calculates the log probability for a given set of parameters.

#     The array should contain the following elements in order:
#         - **delta**: A float representing a prior value
#         - **K**: An integer representing a prior value
#         - **n1**: An integer representing a prior value
#         - **n2**: An integer representing a prior value
#         - **gamma_min**: A float representing a prior value
#         - **gamma_max**: A float representing a prior value
#         - **gamma_break**: A float representing a prior value
#         - **B**: A float representing a prior value
#         - **R**: A float representing a prior value
#         - **\*other_params**: Additional parameters (optional)

#     :param params: A list of parameters.
#     :type params: list
#     :return: The log probability for the given parameters. It returns -inf if any of the following conditions are met:

#         - **n1 > n2**
#         - **gamma_min > gamma_max**
#         - **gamma_break < gamma_min**
#         - **gamma_break > gamma_max**
#         - Any parameter falls outside the specified range (param_min_vals and param_max_vals)

#     :rtype: float
#     """
#     name_stem = "run_" + str(random.getrandbits(60))

#     # prior
#     delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R, *other_params = params
#     if (
#         n1 > n2
#         or gamma_min > gamma_max
#         or gamma_break < gamma_min
#         or gamma_break > gamma_max
#     ):
#         return -np.inf
#     # testing if between min and max
#     for i in range(len(params)):
#         # gamma break is solely constrained by gamma_min and gamma_max
#         if i != 6 and (param_min_vals[i] > params[i] or param_max_vals[i] < params[i]):
#             return -np.inf
#     if use_variability:
#         tau_var = tau * 60 * 60  # to seconds
#         c = 2.997924 * 1.0e10
#         R = np.power(10, R)
#         if tau_var < (1 + redshift) / c * R / delta:
#             return -np.inf

#     model = make_model(params, name_stem)

#     # calculate chi squared
#     func = interpolate.interp1d(model[0], model[1], fill_value="extrapolate")
#     chi_sq = np.sum(
#         ((vFv_data - np.power(10, func(np.log10(v_data)))) / err_data) ** 2.0
#     )
#     for f in glob.glob(BASE_PATH + DATA_FOLDER + "/" + name_stem + "_*"):
#         os.remove(f)

#     return -0.5 * chi_sq


def mcmc(p0=None):
    """
    :param p0: The initial parameter values for each walker. If not provided, random parameter values will be generated.
    :type p0: numpy.array or None
    :return: A tuple containing the sampler object and the directory where the results are saved.
    :rtype: tuple[emcee.EnsembleSampler, str]

    This function performs a Markov Chain Monte Carlo (MCMC) simulation using the emcee library in Python. The MCMC simulation is used to sample from the posterior distribution of a specified likelihood function. The simulation starts from an initial set of parameter values (p0) and iteratively updates these values to explore the parameter space and converge to the posterior distribution.

    If p0 is provided, the function prints the first element of p0.
    The function then checks for "folder_label" and "description" in the "configs" dictionary. If "folder_label" is present, it assigns its value to the variable "folder_label", otherwise it assigns the default value "run". If "description" is present, it assigns its value to the variable "description", otherwise it assigns None.

    The function then creates a directory to store the results using the current date and time.
    It creates a "backend.h5" file for storing the backend data.
    It creates a "basic_info.txt" file to store basic information about the simulation, including the folder name, description (if provided), and configurations. If p0 is provided, it indicates that p0 was given in the file.

    If the "parallel" configuration is set to True in the "configs" dictionary, the function creates a multiprocessing pool to run the simulation in parallel. The number of processes in the pool is determined by the "cores" configuration in "configs" (default value is the number of available CPU cores). If "parallel" is False, the pool is set to None.

    The function then creates an emcee HDFBackend object to store the backend data.
    It resets the backend with the specified number of walkers and dimensions.

    If p0 is not provided, it generates random initial parameter values for each walker.

    The function creates an emcee EnsembleSampler object with the specified number of walkers, dimensions, log probability function, backend, move, and pool.

    It prints a message indicating the start of the MCMC simulation.

    The function records the start time and runs the MCMC simulation using the run_mcmc() method of the sampler.

    If the "parallel" configuration is True, it closes the multiprocessing pool.

    The function records the end time and appends the total simulation time to the "basic_info.txt" file.

    Finally, the function calls the "show_results()" and "save_plots_and_info()" functions from the blazar_report module with various parameters to display and save the simulation results.

    The function returns a tuple containing the MCMC sampler and the directory where the results are saved.
    """
    if p0 is not None:
        print(p0[0])
    if "folder_label" in configs:
        folder_label = configs["folder_label"]
    else:
        folder_label = "run"
    if "description" in configs:
        description = configs["description"]
    else:
        description = None

    now = datetime.datetime.now()
    date_string = now.strftime("%Y-%m-%d-%H:%M:%S")
    directory = RESULTS_FOLDER + "/" + folder_label + "_" + date_string

    os.mkdir(BASE_PATH + directory)

    backend = directory + "/backend.h5"

    # make file with basic info
    with open(BASE_PATH + directory + "/basic_info.txt", "w") as f:
        f.write("folder name: ")
        f.write(directory)
        if description is not None:
            f.write("\nreport description: ")
            f.write(description)
            f.write("\n")
        f.write("\nconfigurations:\n")
        f.write(str(configs))
        if p0 is not None:
            f.write("\np0 given\n")

    if configs["parallel"]:
        if "cores" not in configs:
            pool = multiprocessing.Pool()
        else:
            pool = multiprocessing.Pool(processes=configs["cores"])
    else:
        pool = None
    backend = emcee.backends.HDFBackend(BASE_PATH + backend)
    backend.reset(configs["n_walkers"], DIM)

    if p0 is None:
        p0 = np.array([random_params() for _ in range(configs["n_walkers"])])

    sampler = emcee.EnsembleSampler(
        configs["n_walkers"],
        DIM,
        log_prob,
        backend=backend,
        moves=[(emcee.moves.StretchMove(live_dangerously=True), 1.0)],
        pool=pool,
    )

    print("starting mcmc")
    start = datetime.datetime.now()
    sampler.run_mcmc(p0, configs["n_steps"], progress=True)
    end = datetime.datetime.now()
    if configs["parallel"]:
        pool.close()

    with open(BASE_PATH + directory + "/basic_info.txt", "a") as f:
        f.write("\ntime: ")
        f.write(str(end - start))
        f.write("\n")

    blazar_report.show_results(sampler, str(end - start), configs=configs)
    blazar_report.save_plots_and_info(
        configs,
        (v_data, vFv_data, err_data),
        param_min_vals,
        param_max_vals,
        folder=directory,
        sampler=sampler,
        use_sampler=True,
        description=description,
        time=str(end - start),
        redshift=redshift,
        eic=EIC,
    )

    return sampler, directory


def main_cli():
    """
    This function is the entry point for the command-line interface of the software. It performs the following steps:

    1. It initializes the variable `p0` to `None`.

    2. It calls the `mcmc()` function with the `p0` variable as an argument. The `mcmc()` function returns two values - `sampler` and `directory`. These values are assigned to the variables `sampler` and `directory` respectively.

    3. If the value of `blazar_properties.TMP` is true, it moves the directory specified by `BASE_PATH + directory` to `blazar_properties.FOLDER_PATH + directory` using the `shutil.move()` function. It also removes the temporary directory specified by `blazar_properties.TEMP_DIR` using the `shutil.rmtree()` function.

    :return: None
    :rtype: None
    """
    # p0_file = "local_results/3C66A_b6_eic_2022-06-08-20:17:26/backend.h5"
    # reader = emcee.backends.HDFBackend(BASE_PATH + p0_file, read_only=True)
    # p0 = reader.get_last_sample().coords
    p0 = None
    sampler, directory = mcmc(p0=p0)
    if blazar_properties.TMP:
        shutil.move(BASE_PATH + directory, blazar_properties.FOLDER_PATH + directory)
        shutil.rmtree(blazar_properties.TEMP_DIR)


if __name__ == "__main__":
    main_cli()
