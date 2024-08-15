#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file name: blazar_model.py

Contains functions to interact with the c++ code to create the SED models. In
separate functions, it creates the necessary parameter file, calls the c++ code,
reads the data saved by the c++ code, and processes the model data into a
usable format.
Then, function make_model does all of this when called with parameters, which is
the function that is used by the rest of the code.

Note: All file paths are relative to blazars-mcmc

Parameters are always listed in the following order:
[delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R]
The additional parameters for EIC are bb_temp, l_nuc, tau, blob_dist, in that order.
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
R           blob radius [cm]                log
---------------------------------------------------
*Additional params for EIC*
bb_temp     Black body temp of disk [K]     log
l_nuc       Nucleus luminosity [ergs/s]     log
tau         Frac of luminosity scattered    log
blob_dist   Distance of blob [cm]           log
---------------------------------------------------

Contains:
make_SED(params, name_stem=None, theta=None, redshift=None, min_freq=None, max_freq=None, torus_temp=None,
             torus_luminosity=None, torus_frac=None, data_folder=None, executable=None, command_params_full=None,
             command_params_1=None, command_params_2=None, prev_files=False, verbose=False, eic=False)
file_make_SED(parameter_file=None, data_folder=None, executable=None, prev_files=False, verbose=False)
add_data(current_data, new_data=None, file_suffix=None, name_stem=None, data_folder=None, cols=(0, 2))

FUNCTIONS
"""
import os

import numpy as np
from scipy import interpolate

from blazar_mcmc.blazar_properties import *
from bjet_core import bj_core

__all__ = [
    "add_data",
    "command_line_sub_strings",
    "create_params_file",
    "file_make_SED",
    "make_model",
    "make_SED",
    "params_linear_to_log",
    "params_log_to_linear",
    "process_model",
]


def make_SED(
    params,
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
    prev_files=False,
    verbose=False,
    eic=False,
    folder=None,
):
    """
    Arguments:
        params: 1D np array of NUM_DIM floats
            parameters
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
               prev_files (optional): bool
            Whether bjet should create _prev files after each run; default is False
        verbose (optional): bool
            Whether information on the model should be displayed; default is
            False
        eic: bool
            states whether the run is eic or std; default is false (std)
    """
    # TODO rewrite this entire method using bj_core methods and without building a command string.
    if executable is None:
        executable = CPP_FOLDER + "/" + EXECUTABLE

    if command_params_full is None:
        params = params_log_to_linear(
            params, param_is_log=modelProperties(eic).PARAM_IS_LOG
        )
        # if substrings are not provided
        if command_params_1 is None or command_params_2 is None:
            command_params_1, command_params_2 = command_line_sub_strings(
                name_stem=name_stem,
                theta=theta,
                redshift=redshift,
                min_freq=min_freq,
                max_freq=max_freq,
                data_folder=data_folder,
                prev_files=prev_files,
                eic=eic,
            )
        # now create the params from the substrings
        if not eic:
            command_params_full = (
                command_params_1 + [str(val) for val in params] + command_params_2
            )
        else:
            if torus_temp is None:
                torus_temp = "2.0e+4"
            else:
                torus_temp = str(torus_temp)
            if torus_luminosity is None:
                torus_luminosity = "5.5e+20"
            else:
                torus_luminosity = str(torus_luminosity)
            if torus_frac is None:
                torus_frac = "9.0e-5"
            else:
                torus_frac = str(torus_frac)

            # this gives you the settings commands (command_params_1), the initial 9 parameters + length of emitting region (spherical) + absorption by EBL = yes, blob distance
            command_params_full = (
                command_params_1
                + [str(val) for val in params[:9]]
                + command_params_2[:2]
                + [str(params[-1])]
            )
            # add eic parameter values
            command_params_full = command_params_full + [
                str(params[9]),
                torus_temp,
                str(params[10]),
                str(params[11]),
                torus_luminosity,
                torus_frac,
            ]

            # add numerical params (min/max freq, number of pts)
            command_params_full = command_params_full + command_params_2[3:]

    # command_params_full have been set now

    if not verbose:
        subprocess.run(
            [BASE_PATH + executable, *command_params_full],
            stderr=open(os.devnull, "wb"),
            stdout=open(os.devnull, "wb"),
        )
    else:
        result = subprocess.run(
            [BASE_PATH + executable, *command_params_full],
            capture_output=True,
            text=True,
        )
        print("params", params)
        print(command_params_full)

        # write log file and print
        print(result.stderr)
        with open(FOLDER_PATH + folder + "/bjet.log", "w") as f:
            f.write(result.stderr)


def file_make_SED(
    parameter_file=None,
    data_folder=None,
    executable=None,
    prev_files=False,
    verbose=False,
):
    """
    Calls the C++ code to generate the SED.
    Compton model data is saved to <home directory>/DATA_FOLDER/<file prefix>_cs.dat
    Sync model data is saved to <home directory>/DATA_FOLDER/<file prefix>_ss.dat
    Args:
        parameter_file (optional): string
            Relative path to param file, will be overwritten with parameters; default
            is <PARAMETER_FOLDER>/params.txt.
        executable (optional): string
            Relative path to executable; default is B_jet02/bj
        prev_files (optional): bool
            Whether bjet should create _prev files
        verbose (optional): bool
            Whether the output of the bjet code is shown; default is false
    After function call:
        Data saved to <home directory>/DATA_FOLDER/<file prefix>_*.dat
        Parameter file overwritten or created
    """
    # TODO rewrite this entire method using bj_core methods and without building a command string.
    if parameter_file is None:
        parameter_file = PARAMETER_FOLDER + "/params.txt"
    if executable is None:
        executable = CPP_FOLDER + "/" + EXECUTABLE
    if data_folder is None:
        data_folder = DATA_FOLDER

    # Call C++ code to make the SED and
    # store it in <home directory>/DATA_FOLDER with given file name prefix
    if prev_files:
        p = "2"
    else:
        p = "3"
    # input mode, data file
    settings_params = [p, BASE_PATH + data_folder]
    if verbose:
        subprocess.run(
            [BASE_PATH + executable, *settings_params, BASE_PATH + parameter_file]
        )
    else:
        subprocess.run(
            [BASE_PATH + executable, *settings_params, BASE_PATH + parameter_file],
            stderr=open(os.devnull, "wb"),
            stdout=open(os.devnull, "wb"),
        )


def add_data(
    current_data,
    new_data=None,
    file_suffix=None,
    name_stem=None,
    data_folder=None,
    cols=(0, 2),
):
    if name_stem is None:
        name_stem = NAME_STEM
    if data_folder is None:
        data_folder = DATA_FOLDER

    current_logv, current_logvFv, current_v, current_vFv = current_data
    if new_data is not None:
        model_logv, model_logvFv, model_v, model_vFv = new_data
    elif file_suffix is not None:
        try:
            model_data = np.loadtxt(
                BASE_PATH + data_folder + "/" + name_stem + "_" + file_suffix + ".dat",
                delimiter=" ",
            )
        except IOError:
            raise IOError(
                "add_data: "
                + BASE_PATH
                + data_folder
                + "/"
                + name_stem
                + "_"
                + file_suffix
                + ".dat cannot be read"
            )

        if name_stem == "F_jet":
            model_logv = model_data[:, 3]
            model_logvFv = model_data[:, 5]
            model_vFv = model_data[:, 2]
        else:
            model_logv = model_data[:, cols[0]]
            model_logvFv = model_data[:, cols[1]]
            model_vFv = np.power(10, model_logvFv)
    else:
        raise ValueError("In add_data, data or file suffix not provided")

    """
    Find overlap
    """
    new_lower = np.where(model_logv < current_logv[0])[0]
    new_higher = np.where(model_logv > current_logv[-1])[0]

    logv = np.concatenate((model_logv[new_lower], current_logv, model_logv[new_higher]))
    logvFv = np.concatenate(
        (model_logvFv[new_lower], current_logvFv, model_logvFv[new_higher])
    )
    vFv = np.power(10, logvFv)

    overlap_start = np.where(logv >= max(model_logv[0], current_logv[0]))[0][0]
    overlap_end = np.where(logv <= min(model_logv[-1], current_logv[-1]))[0][-1]

    interpolation = interpolate.interp1d(model_logv, model_vFv)

    new_vFv = np.concatenate(
        (
            np.zeros(overlap_start),
            interpolation(logv[overlap_start : overlap_end + 1]),
            np.zeros(len(logv) - overlap_end - 1),
        )
    )

    vFv = vFv + new_vFv
    logvFv = np.log10(vFv)
    v = np.power(10, logv)
    return logv, logvFv, v, vFv


def process_model(
    name_stem=None, data_folder=None, verbose=False, eic=False, additional_suffixes=None
):
    """
    Read a model from data files and returns arrays of frequencies and
        energy flux.
    This model only uses data from the *compton model* and the *synchrotron
        model*.
    The data we want are an array of frequencies (v) and an array of
        corresponding flux energies (vFv). Flux energy is the sum of the flux
        energy from the synchrotron model and that from the compton model when
        there are data points for both of them.

    Requirements:
        The SED model must have already been created (using make_SED above)
            with the given name_stem.
        By default, data is read from DATA_FOLDER/<name_stem>_*.dat, the location
            that make_SED writes to.

    Args:
        name_stem (str): model data is saved in files named in the form
            <name_stem>_*.dat, specify file name
        data_folder (str): relative path to the folder containing the data files

    Returns:
        a tuple of 4 1D numpy arrays

        all arrays have the same length.
        logv_all, logvFv_all, v_all, vFv_all

        logv_all: a 1D numpy array of the log of all frequencies
        logvFv_all: a 1D numpy array of the log of all the energy flux
        v_all: a 1D numpy array of all frequencies in the data
        vFv_all: a 1D numpy array of corresponding energy fluxes for the
            frequencies

    Raises:
        IOError when data cannot be read

    Note: The frequency values used are the ones used in the synchrotron spectrum
    and the ones in the compton model greater than the max in the synchrotron. The
    frequency values present in the Compton model in the overlap are not used.
    """

    if name_stem is None:
        name_stem = NAME_STEM
    if data_folder is None:
        data_folder = DATA_FOLDER
    # Read in models from the files created by the c++ code
    # This model uses only the data from the compton model and the synchrotron
    # model.
    # Data is read into a numpy 2D array
    try:
        synchrotron_model = np.loadtxt(
            BASE_PATH + data_folder + "/" + name_stem + "_ss.dat", delimiter=" "
        )
    except IOError:
        raise IOError(
            "process_model: "
            + BASE_PATH
            + data_folder
            + "/"
            + name_stem
            + "_*.dat cannot be read"
        )

    # extract data used - only use frequency and energy flux, which are the 1st
    # and 3rd columns
    logv_synchrotron = synchrotron_model[:, 0]
    logvFv_synchrotron = synchrotron_model[:, 2]
    v_synchrotron = np.power(10, logv_synchrotron)
    vFv_synchrotron = np.power(10, logvFv_synchrotron)

    logv, logvFv, v, vFv = add_data(
        (logv_synchrotron, logvFv_synchrotron, v_synchrotron, vFv_synchrotron),
        file_suffix="cs",
        name_stem=name_stem,
        data_folder=data_folder,
    )
    logv, logvFv, v, vFv = add_data(
        (logv, logvFv, v, vFv),
        file_suffix="cs2",
        name_stem=name_stem,
        data_folder=data_folder,
    )
    # Total energy flux is determined by summing the flux from the compton and
    # synchrotron models.

    if eic:
        if additional_suffixes is None:
            additional_suffixes = ["ecs", "nuc"]
        if verbose:
            print("process_model EIC mode: ss cs", *additional_suffixes, "read")
        for s in additional_suffixes:
            logv, logvFv, v, vFv = add_data(
                (logv, logvFv, v, vFv),
                file_suffix=s,
                name_stem=name_stem,
                data_folder=data_folder,
            )
    elif verbose:
        print("process_model SSC mode")

    return logv, logvFv, v, vFv


def make_model(
    params,
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
    use_param_file=False,
    verbose=False,
    eic=False,
    fixed_params=None,
    folder=None,
):
    """
    Given parameters, returns v and vFv for the corresponding model.

    Args:
        params (list of floats): a list of NUM_DIM floats which are the parameters
        name_stem (str): data will be saved in files of the form
            <name_stem>_*.dat
        parameter_file (str): optional; relative path of the parameter file
        theta: optional; params is fixed for modeling
        redshift: optional; redshift is fixed for modeling
        min_freq (float): optional; specifies min frequency model SEDs should use;
            default value used if none provided, which is 5.0e+7.
        max_freq (float): optional; specifies max frequency model SEDs should use;
            default value used if none provided, which is 1.0e+29.
        executable (str): optional; where the bjet executable is located
        data_folder (str): relative path to the folder where data will be saved to
        use_param_file (bool): specifies if bjet should be called with parameter file or with command line args
        verbose (bool): specifies if parameters and output of bjet should be shown

    After function call: parameter_files is over-written with parameter information.
    Files of the form <name_stem>_*.dat are created inside <path>

    Returns:
        a tuple of 4 1D numpy arrays

        all arrays have the same length.
        logv_all, logvFv_all, v_all, vFv_all

        logv_all: a 1D numpy array of the log of all frequencies
        logvFv_all: a 1D numpy array of the log of all the energy flux
        v_all: a 1D numpy array of all frequencies in the data
        vFv_all: a 1D numpy array of corresponding energy fluxes for the
            frequencies

    Raises:
        IOError if parameter_files cannot be written to (likely, a folder in the path does not exist)
    """
    if use_param_file:
        if eic:
            raise ValueError("EIC mode cannot be used with parameter file")
        create_params_file(
            params,
            name_stem,
            parameter_file,
            theta=theta,
            redshift=redshift,
            min_freq=min_freq,
            max_freq=max_freq,
            verbose=verbose,
        )
    # reimplement fixed parameters for the model computation
    if fixed_params:
        for i in range(len(fixed_params)):
            if fixed_params[i] != -np.inf:
                params = np.insert(params, i, fixed_params[i])
    make_SED(
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
        prev_files=prev_files,
        verbose=verbose,
        eic=eic,
        folder=folder,
    )
    return process_model(name_stem=name_stem, data_folder=data_folder, eic=eic)


def command_line_sub_strings(
    name_stem=None,
    theta=None,
    redshift=None,
    min_freq=None,
    max_freq=None,
    data_folder=None,
    prev_files=False,
    eic=False,
):
    if name_stem is None:
        name_stem = NAME_STEM
    if theta is None:
        theta = 0.57
    if redshift is None:
        print("here in command_line_sub_strings")
        pass
    if min_freq is None:
        min_freq = 1.0e08
    if max_freq is None:
        max_freq = 1.0e28
    if data_folder is None:
        data_folder = DATA_FOLDER
    if prev_files:
        p = "2"
    else:
        p = "3"

    # input type, data file, model specification, redshift, hubble const, angle
    if eic:
        model_type = "1"
    else:
        model_type = "0"
    settings_and_transformation = [
        p,
        BASE_PATH + data_folder,
        model_type,
        str(redshift),
        "69.6",
        str(theta),
    ]

    # length of emitting region, absorption by EBL, distance of blob (host galaxy frame)
    # number of spectral points, minimal frequency, maximal frequency, file name prefix
    constant_and_numerical = [
        "0",
        "1",
        "9.0e+17",
        "99",
        str(min_freq),
        str(max_freq),
        name_stem,
    ]
    return settings_and_transformation, constant_and_numerical


def create_params_file(
    params,
    name_stem=None,
    parameter_file=None,
    theta=None,
    redshift=None,
    min_freq=None,
    max_freq=None,
    verbose=False,
):
    """
    NOTE: Cannot be used with EIC
    Given parameters and a file, put the data in the file in the format
    understandable by the C++ code.

    Args:
        params: List of NUM_DIM floats
            Form: [delta (linear), K (log), n1 (linear), n2 (linear), gamma_min (log),
                gamma_max (log), gamma_break (log), B (log), R (log)]
        name_stem: str
            Data will be saved in files named <name_stem>_*.dat
        parameter_file (optional): str
            String with the relative path of the parameter file; default is None;
            will be set to "<PARAMETER_FOLDER>/params.txt"
        min_freq (optional): float
            Specifies min frequency model SEDs should use; default is None;
            will be set to 5.0e+7
        max_freq (optional): float
            Specifies max frequency model SEDs should use; default is None;
            will be set to 1.0e+29
        redshift (optional): float
            Redshift is fixed for modeling; default is None; will be set to 0.57
        theta (optional): float
            params is fixed for modeling; default is None; will be set to 0.143
        verbose (optional):
            Specifies if parameters are shown; default is False
    After function call:
        parameter_files is over-written (and createed if necessary) with parameter information.

    Raises:
        IOError if parameter_files cannot be written to (likely, a folder in the path does not exist)
    """
    if parameter_file is None:
        parameter_file = PARAMETER_FOLDER + "/params.txt"
    if name_stem is None:
        name_stem = NAME_STEM
    if theta is None:
        theta = 0.57
    if redshift is None:
        redshift = 0.143
    if min_freq is None:
        min_freq = 5.0e7
    if max_freq is None:
        max_freq = 1.0e29
    if verbose:
        print(params)

    # Set parameters from list, convert those in log space to linear space
    (
        doppler_param,
        K_param,
        n1_param,
        n2_param,
        gamma_min_param,
        gamma_max_param,
        gamma_break_param,
        B_param,
        R_param,
    ) = params_log_to_linear(params, eic=False)

    if verbose:
        print("z (fixed) =", redshift)
        print("params (fixed) =", theta)
        print("doppler =", doppler_param)
        print("K =", K_param)
        print("n1 =", n1_param)
        print("n2 =", n2_param)
        print("Gamma min =", gamma_min_param)
        print("Gamma max =", gamma_max_param)
        print("Gamma break =", gamma_break_param)
        print("B =", B_param)
        print("R =", R_param)

    # write data with the parameters in the format the C++ code can understand
    # into the given address for parameter_files
    try:
        with open(BASE_PATH + parameter_file, "w") as file:
            file.write(
                "------------------------------------------------------------------------------ \n"
            )
            file.write("Transformation_parameters \n")
            file.write(
                "------------------------------------------------------------------------------ \n"
            )
            file.write(str(redshift) + "          redshift \n")
            file.write(
                "69.6           Hubble_constant______________________________________[km/(sMpc)] \n"
            )
            file.write(
                str(theta)
                + "            angle_to_the_line_of_sight___________________________[degrees] \n"
            )
            file.write(
                "------------------------------------------------------------------------------ \n"
            )
            file.write("Blob_parameters \n")
            file.write(
                "------------------------------------------------------------------------------ \n"
            )
            file.write(str(doppler_param) + "             Doppler_factor \n")
            file.write(
                str(K_param)
                + "        Particle_density_____________________________________[1/cm^3] \n"
            )
            file.write(
                str(n1_param) + "            First_slope_of_particle_energy_spectrum \n"
            )
            file.write(
                str(n2_param)
                + "            Second_slope_of_particle_energy_spectrum \n"
            )
            file.write(str(gamma_min_param) + "            Minimum_electrons_energy \n")
            file.write(str(gamma_max_param) + "          Maximum_electrons_energy \n")
            file.write(
                str(gamma_break_param) + "        Break_in_electrons_energy_spectrum \n"
            )
            file.write(
                str(B_param)
                + "            Magnetic_field_______________________________________[G] \n"
            )
            file.write(
                str(R_param)
                + "        Radius_of_emitting_region____________________________[cm] \n"
            )
            file.write(
                "0              length_of_emitting_region_(0_for_spherical_geometry)_[cm] \n"
            )
            file.write("1              absorption_by_EBL_(0=NO___1=YES) \n")
            file.write(
                "9.0e+17        distance_blob_SMBH_(host_galaxy_frame)_______________[cm] \n"
            )
            file.write(
                "------------------------------------------------------------------------------ \n"
            )
            file.write("Extern_Inverse_Compton_parameter \n")
            file.write(
                "------------------------------------------------------------------------------ \n"
            )
            file.write("0              compute_EIC_(0=NO___1=YES) \n")
            file.write("0              compute_X_corona_(0=NO___1=YES) \n")
            file.write(
                "4.0e+4         Disk_black_body_temperature__________________________[K] \n"
            )
            file.write(
                "2.0e+4         Torus_black_body_temperature_________________________[K] \n"
            )
            file.write(
                "3.0e+43        Luminosity_of_the_disk_______________________________[erg/s] \n"
            )
            file.write(
                "9.0e-5         Tau___fraction_of_L_disk_reprocessed_isotropically \n"
            )
            file.write(
                "5.5e+20        Luminosity_of_the_torus______________________________[erg/s] \n"
            )
            file.write(
                "9.0e-5         Tau___fraction_of_L_tor_reprocessed_isotropically \n"
            )
            file.write(
                "------------------------------------------------------------------------------ \n"
            )
            file.write("Jet_parameters  \n")
            file.write(
                "------------------------------------------------------------------------------ \n"
            )
            file.write("0              compute_JET_(0=NO___1=YES) \n")
            file.write("3.0            Doppler_factor \n")
            file.write(
                "1.0e+2         Initial_particle_density_____________________________[1/cm^3] \n"
            )
            file.write("2.1            Slope_of_particle_energy_spectrum \n")
            file.write("1.0e+3         Minimum_electrons_energy \n")
            file.write("2.0e+4         Maximum_electrons_energy \n")
            file.write(
                "0.08           Initial_magnetic_field_______________________________[G] \n"
            )
            file.write(
                "1.2e+17        Inner_radius_(host_galaxy_frame)_____________________[cm] \n"
            )
            file.write(
                "300            Jet_length_(host_galaxy_frame)_______________________[pc] \n"
            )
            file.write(
                "1.0            Half-opening_angle_of_jet_(host_galaxy_frame)________[deg] \n"
            )
            file.write("50             number_of_slices \n")
            file.write(
                "------------------------------------------------------------------------------ \n"
            )
            file.write("Numerical_parameters \n")
            file.write(
                "------------------------------------------------------------------------------ \n"
            )
            file.write("99             number_of_spectral_points \n")
            file.write(
                str(min_freq)
                + "        minimal_frequency____________________________________[Hz] \n"
            )
            file.write(
                str(max_freq)
                + "        maximal_frequency____________________________________[Hz] \n"
            )
            file.write(name_stem + "        prefix_of_file_names \n")
    except IOError:
        raise IOError(
            "create_params_file: " + parameter_file + " could not be written to"
        )


def params_log_to_linear(params, param_is_log=None, eic=False):
    if param_is_log is None:
        param_is_log = modelProperties(eic).PARAM_IS_LOG

    ret = params.copy()
    for i in range(len(params)):
        if param_is_log[i]:
            ret[i] = np.power(10, params[i])
    return ret


def params_linear_to_log(params, param_is_log=None, eic=False):
    if param_is_log is None:
        param_is_log = modelProperties(eic).PARAM_IS_LOG

    ret = params.copy()
    for i in range(len(params)):
        if param_is_log[i]:
            ret[i] = np.log10(params[i])
    return ret
