#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file name: blazar_run_mcmc.py

Program Purpose: Implement the mcmc. Run the file to run the MCMC as a script
using a parameter file.

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

Contains:

FUNCTIONS:
run_mcmc
mcmc


run_mcmc(configs, param_min_vals=None, param_max_vals=None, backend_file="<RESULTS_FOLDER>/backend.h5",
             p0=None, min_freq=None, max_freq=None)
    Runs the actual mcmc given a configuration dictionary.
    Returns the sampler and a string of how long the mcmc took.

mcmc(config_file=None, directory=None, p0=None, p0_label=None)
    Function called by main. Can be run with no parameters and runs the mcmc given
    a configuration file.
    Returns the sampler.
"""
import datetime
import multiprocessing
import os
import sys
import platform
import shutil
import numpy as np
import emcee

from bjet_mcmc import blazar_model
from bjet_mcmc import blazar_report
from bjet_mcmc import blazar_utils
from bjet_mcmc import blazar_clean
from bjet_mcmc import blazar_initialize

from bjet_mcmc.blazar_properties import *

__all__ = ["run_mcmc", "mcmc"]


def run_mcmc(
    configs,
    param_min_vals=None,
    param_max_vals=None,
    backend_file=RESULTS_FOLDER + "/backend.h5",
    p0=None,
    p0_file=None,
    eic_p0_from_std=False,
    min_freq=None,
    max_freq=None,
    prev_files=False,
    use_param_file=True,
    eic=False,
):
    """
    Arguments:
        configs: dictionary; keys are strings, value types are mixed
            Configurations usually obtained with read_configs. See below for
                necessary dictionary elements.
        param_min_vals (optional): np array of floats, shape (NUM_DIM,)
            min values for each param in the usual order; default is None; will
            use defaults in blazar_utils.min_max_parameters
        param_max_vals (optional): np array of floats, shape (NUM_DIM,)
            max values for each param in the usual order; default is None; will
            use defaults in blazar_utils.min_max_parameters
        backend_file (optional): str
            location of backend file; default location will be used if none specified
        p0 (optional): 2D np array of floats, shape (n_walkers, NUM_DIM)
            initial locations. Default is none; then, random start parameters will be used.
            Usually only specify this if you would like to continue a previous run from the
            last state of a chain.
        min_freq (optional): float
            minimum frequency for creating SED models; default is None, then default
            will be used.
        max_freq (optional): float
            maximum frequency for creating SED models; default is None, then default
            will be used.
    ------
    Returns:
        tuple: (emcee.EnsembleSampler, str)
        (sampler, time)
        sampler: sampler created
        time: string of how long the MCMC took in hh:mm:ss.ssssss
    ------

    Config dictionary keys and values:
    Keys are strings, value format is mixed
    Key (str)               Value
    "data_file"             (str) relative path to data
    "n_steps"               (int) number of num to run
    "n_walkers"             (int) number of walkers
    "discard"               (int) how many num are discarded before analysis
    "parallel"              (bool) whether parallel processing will be used
    "cores"                 (int) NOT PRESENT IF PARALLEL IS FALSE; # of cores
                                if absent and parallel TRUE,
                                (# cores on machine - 1) used
    "tau_variability"       (float) time in hours for variability constraint
    """
    if configs["parallel"]:
        if "cores" not in configs:
            pool = multiprocessing.Pool()
        else:
            pool = multiprocessing.Pool(processes=configs["cores"])
    else:
        pool = None
    backend = emcee.backends.HDFBackend(BASE_PATH + backend_file)
    backend.reset(
        configs["n_walkers"],
        modelProperties(eic, fixed_params=configs["fixed_params"]).NUM_DIM,
    )
    # set minima and maxima for the parameter values

    if param_min_vals is None or param_max_vals is None:
        default_min_max = blazar_utils.min_max_parameters(
            alpha2_limits=configs["alpha2_limits"],
            eic=eic,
            fixed_params=configs["fixed_params"],
        )
    if param_min_vals is None:
        param_min_vals = default_min_max[0]
    if param_max_vals is None:
        param_max_vals = default_min_max[1]

    bjet_command_1, bjet_command_2 = blazar_model.command_line_sub_strings(
        min_freq=min_freq,
        max_freq=max_freq,
        prev_files=prev_files,
        redshift=configs["redshift"],
        eic=eic,
    )
    # read data
    v_data, vFv_data, err_data = blazar_utils.read_data(configs["data_file"])

    # starting position, a numpy array with # of walkers rows and # of parameters columns
    if eic_p0_from_std:
        if p0_file is None and p0 is None:
            raise ValueError("eic_p0_from_std is True but no p0 source provided")
        if not eic:
            raise ValueError("eic_p0_from_std is True, but mode is not EIC")
        if p0_file is not None:
            reader = emcee.backends.HDFBackend(BASE_PATH + p0_file, read_only=True)
            p0 = reader.get_last_sample().coords
        p0 = blazar_utils.random_eic_from_std(
            p0,
            configs["n_walkers"],
            param_min_vals,
            param_max_vals,
            redshift=configs["redshift"],
            tau_var=configs["tau_variability"],
            use_variability=configs["use_variability"],
        )
        print("p0 eic from std")
    elif p0_file is not None:
        reader = emcee.backends.HDFBackend(BASE_PATH + p0_file, read_only=True)
        p0 = reader.get_last_sample()
        if len(p0.coords) != configs["n_walkers"]:
            raise ValueError(
                "p0 start has the incorrect number of walkers. p0 has "
                + str(len(p0.coords))
                + " but the number of walkers is "
                + configs["n_walkers"]
                + "."
            )
        if configs["eic"] and len(p0.coords[0]) != 13:
            raise ValueError(
                "for an EIC run, p0 coordinates must have 13 values. "
                + str(len(p0.coords[0]))
                + " values present."
            )
        elif not configs["eic"] and len(p0.coords[0]) != 9:
            raise ValueError(
                "for a non-EIC run, p0 coordinates must have 9 values. "
                + str(len(p0.coords[0]))
                + " values present."
            )
        print("p0 from file")
    elif p0 is not None:
        if len(p0) != configs["n_walkers"]:
            raise ValueError(
                "p0 start has the incorrect number of walkers. p0 has "
                + str(len(p0.coords))
                + " but the number of walkers is "
                + configs["n_walkers"]
                + "."
            )
        if configs["eic"] and len(p0[0]) != 13:
            raise ValueError(
                "for an EIC run, p0 coordinates must have 13 values. "
                + str(len(p0[0]))
                + " values present."
            )
        elif not configs["eic"] and len(p0[0]) != 9:
            raise ValueError(
                "for a non-EIC run, p0 coordinates must have 9 values. "
                + str(len(p0[0]))
                + " values present."
            )
        print("p0 from values")
    else:
        p0 = blazar_utils.random_defaults(
            configs["n_walkers"],
            param_min_vals,
            param_max_vals,
            redshift=configs["redshift"],
            tau_var=configs["tau_variability"],
            use_variability=configs["use_variability"],
            eic=eic,
            fixed_params=configs["fixed_params"],
        )
        print("p0 from random")
        for i in range(len(configs["fixed_params"])):
            if configs["fixed_params"][i] != -np.inf and configs["eic"]:
                print(EIC_PARAM_NAMES[i], "is fixed at", configs["fixed_params"][i])
            elif configs["fixed_params"][i] != -np.inf:
                print(SSC_PARAM_NAMES[i], "is fixed at", configs["fixed_params"][i])
    # create sampler
    sampler = emcee.EnsembleSampler(
        configs["n_walkers"],
        modelProperties(eic, fixed_params=configs["fixed_params"]).NUM_DIM,
        blazar_utils.log_probability,
        args=[v_data, vFv_data, err_data],
        kwargs={
            "param_min_vals": param_min_vals,
            "param_max_vals": param_max_vals,
            "unique_name": True,
            "tau_var": configs["tau_variability"],
            "use_variability": configs["use_variability"],
            "redshift": configs["redshift"],
            "min_freq": min_freq,
            "max_freq": max_freq,
            "prev_files": prev_files,
            "use_param_file": use_param_file,
            "command_params_1": bjet_command_1,
            "command_params_2": bjet_command_2,
            "torus_temp": None,
            "torus_luminosity": None,
            "torus_frac": None,
            "eic": eic,
            "fixed_params": configs["fixed_params"],
        },
        backend=backend,
        moves=[(emcee.moves.StretchMove(a=2, live_dangerously=True), 1.0)],
        pool=pool,
    )

    # moves set up by default moves=[(emcee.moves.StretchMove(live_dangerously=True), 1.)]
    # can try different moves to better probe the parameter space, but no real improvement so far after several tries
    # https://emcee.readthedocs.io/en/stable/user/moves/#moves-user

    # make pre-run txt file with basic info
    print("starting mcmc")
    start = datetime.datetime.now()
    sampler.run_mcmc(p0, configs["n_steps"], progress=True)
    end = datetime.datetime.now()
    if configs["parallel"]:
        pool.close()

    return sampler, str(end - start)


def mcmc(
    config_file=None,
    directory=None,
    folder_label=None,
    p0=None,
    p0_label=None,
    p0_file=None,
    eic_p0_from_std=False,
    description=None,
    prev_files=False,
    use_param_file=False,
):
    """
    The MCMC function that will run just given a configuration file--this is the
    function called by the main method.

    Arguments:
        all arguments are optional
        config_file (optional): str
            Absolute path to configuration file; default is None; if none,
            it will be set to "mcmc_config.txt"
        directory (optional): str
            Location for results; default is none; if none, it will
            be in <RESULTS_FOLDER>/run_YYYY-mm-dd-hh:mm:ss. This will be set with
            configs if possible.
        folder_label (optional): str
            If directory is None and folder_label is not None, results will
            be saved in a folder named <folder_label>_YYYY-mm-dd-hh:mm:ss
        p0 (optional): numpy arr of floats, n_walkers x NUM_DIM)
            Initial positions of the walkers; default is None; if none, random
            positions will be used.
            Usually only set this if continuing a run from a past MCMC run.
        p0_label (optional): str
            Label describing where the p0 values came from; default
            is none; if none, the label will be "random"
    """
    configs = blazar_utils.read_configs(config_file=config_file)

    if description is None and "description" in configs:
        description = configs["description"]

    if description is not None:
        print(description)

    if folder_label is None and "folder_label" in configs:
        folder_label = configs["folder_label"]
    data = blazar_utils.read_data(configs["data_file"])
    param_min_vals, param_max_vals = blazar_utils.min_max_parameters(
        alpha2_limits=configs["alpha2_limits"],
        eic=configs["eic"],
        fixed_params=configs["fixed_params"],
    )

    # file is a folder with the data
    if directory is None:
        now = datetime.datetime.now()
        date_string = now.strftime("%Y-%m-%d-%H:%M:%S")
        if folder_label is None:
            folder_label = "run"
        directory = RESULTS_FOLDER + "/" + folder_label + "_" + date_string
        os.mkdir(BASE_PATH + directory)

    backend = directory + "/backend.h5"

    if p0_label is None:
        if p0_file is not None:
            p0_label = p0_file
        elif p0 is not None:
            p0_label = "from given values"
        else:
            p0_label = "random"
        if eic_p0_from_std:
            p0_label = "eic_p0_from_std " + p0_label

    # make file with basic info
    with open(BASE_PATH + directory + "/basic_info.txt", "w") as f:
        f.write("folder name: ")
        f.write(directory)
        if description is not None:
            f.write("\nreport description: ")
            f.write(description)
            f.write("\n")

        f.write("\nconfig file: ")
        if config_file == None:
            f.write("mcmc_config.txt")
        else:
            f.write(config_file)
        f.write("\nprev_files: ")
        f.write(str(prev_files))
        f.write(", use_param_file: ")
        f.write(str(use_param_file))
        f.write("\n\n")
        f.write("configurations:\n")
        f.write(str(configs))
        f.write("\n\np0: ")
        f.write(p0_label)
        f.write("\n")

    results = run_mcmc(
        configs,
        param_min_vals=param_min_vals,
        param_max_vals=param_max_vals,
        backend_file=backend,
        p0=p0,
        p0_file=p0_file,
        eic_p0_from_std=eic_p0_from_std,
        prev_files=prev_files,
        use_param_file=use_param_file,
        eic=configs["eic"],
    )

    with open(BASE_PATH + directory + "/basic_info.txt", "a") as f:
        f.write("\ntime: ")
        f.write(str(results[1]))
        f.write("\n")

    sampler = results[0]

    blazar_report.show_results(sampler, results[1], configs=configs)

    blazar_report.save_plots_and_info(
        configs,
        data,
        param_min_vals,
        param_max_vals,
        folder=directory,
        sampler=sampler,
        use_sampler=True,
        description=description,
        time=results[1],
        p0_source=p0_label,
        redshift=configs["redshift"],
        eic=configs["eic"],
    )

    return sampler, directory


def main_cli():
    """
    Calling blazar_run_mcmc.py without argument will run mcmc using the configurations in mcmc_config.txt
    To run it with a specific configuration file, the absolute path can be given as first argument
    e.g. $ python blazar_run_mcmc.py /home/my_configfile.txt
    and save data to RESULTS_FOLDER/run_yyyy-mm-dd-hh:mm:ss.
    """
    if INITIALIZE:
        blazar_clean.clean()
        blazar_initialize.initialize()

    if len(sys.argv) > 1:
        config_file = sys.argv[1]
    else:
        config_file = None

    # label_p0 = "From local_results/3C66A_b5_no_eic_2022-07-20-19:42:07/backend.h5"
    label_p0 = None
    p0_values = None

    # file_p0 = "local_results/TON_599_no_eic/backend.h5"
    file_p0 = None
    p0_eic_from_std = False
    sampler_result, results_directory = mcmc(
        config_file=config_file,
        use_param_file=False,
        p0_label=label_p0,
        p0=p0_values,
        p0_file=file_p0,
        eic_p0_from_std=p0_eic_from_std,
    )

    # in case of interrupted process, there may be leftovers in DATA_FOLDER, remove them
    os.chdir(FOLDER_PATH + DATA_FOLDER)
    files = [f for f in os.listdir() if f.startswith("run")]
    for f in files:
        os.remove(f)

    if TMP:
        shutil.move(BASE_PATH + results_directory, FOLDER_PATH + results_directory)
        shutil.rmtree(TEMP_DIR)


if __name__ == "__main__":
    main_cli()
