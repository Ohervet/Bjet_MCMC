#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file name: blazar_model.py
Contains functions to interact with the C++ code to create the SED models. In
separate functions, it creates the necessary parameter file, calls the C++ code,
reads the data saved by the C++ code, and processes the model data into a
usable format.

Then, the function ``make_model`` does all of this when called with parameters, which is
the function that is used by the rest of the code.

.. note::
    All file paths are relative to blazars-mcmc

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
import os
import subprocess
import ctypes

import numpy as np
from scipy import interpolate

from bjet_mcmc.blazar_properties import *
from bjet_core import bj_core

#for debugging only
#import emcee
#import h5py


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
    :param params: 1D numpy array of NUM_DIM floats
    :type params: numpy.ndarray
    :param name_stem: Name stem for make_model. Default is None; will then be set to default.
    :type name_stem: str, optional
    :param theta: Angle from the line of sight. Default is None, and it will be set to the default value of 0.57.
    :type theta: float, optional
    :param redshift: Redshift value; default is None, so the log_prior function will use the default, which is 0.143 (the value for J1010).
    :type redshift: float, optional
    :param min_freq: Minimum frequency for the SED model. Default is None, where it will be set to the default value of 5.0e+7 in blazar_model.process_model.
    :type min_freq: float, optional
    :param max_freq: Maximum frequency for the SED model. Default is None, where it will be set to the default value of 1.0e+29 in blazar_model.process_model.
    :type max_freq: float, optional
    :param torus_temp: Value for the torus temperature. Default is None, and it will be set to the default of 2.0e+4.
    :type torus_temp: float, optional
    :param torus_luminosity: Value for the torus luminosity. Default is None, and it will be set to the default of 5.5e+20.
    :type torus_luminosity: float, optional
    :param torus_frac: Value for the fraction of the torus luminosity reprocessed isotropically. Default is None, and it will be set to the default of 9.0e-5.
    :type torus_frac: float, optional
    :param data_folder: Relative path to the folder where data will be saved to.
    :type data_folder: str, optional
    :param executable: Where the bjet executable is located.
    :type executable: str, optional
    :param command_params_full: Full set of parameters to pass to the bjet exec--see the README for information. Should be length 22, 23, 28, or 29.
    :type command_params_full: numpy.ndarray, optional
    :param command_params_1: The settings and transformation parameters: [prev files flag, data folder, model type, redshift, hubble constant, theta].
    :type command_params_1: numpy.ndarray, optional
    :param command_params_2: The constant and numerical parameters: [length of the emitting region, absorption by EBL, blob distance, # of spectral points, min freq, max freq, file name prefix].
    :type command_params_2: numpy.ndarray, optional
    :param prev_files: Whether bjet should create _prev files after each run; default is False.
    :type prev_files: bool, optional
    :param verbose: Whether information on the model should be displayed; default is False.
    :type verbose: bool, optional
    :param eic: States whether the run is eic or std; default is False (std).
    :type eic: bool
    :param folder: Additional folder information; default is None.
    :type folder: str, optional
    """

    if executable is None:
        executable = CPP_FOLDER + "/" + EXECUTABLE
    #for debug only
    #params_original = params[:]
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

    if not verbose:
        #new method, without subprocess
        
        #create a single string of parameters for ints implentation in C++ with swig
        model_type = 0
        params_char = [0]*len(params)
        for i in range(len(params)):
            params_char[i] = f'{params[i]:0.4e}'
        list_params = command_params_1[3:] + params_char + command_params_2[:-1]
        if eic:
            model_type = 1
            list_params = command_params_1[3:] + params_char[:9] + command_params_2[:2] + [params_char[-1]] + [
               params_char[9],
               str(torus_temp),
               params_char[10],
               params_char[11],
               str(torus_luminosity),
               str(torus_frac)] + command_params_2[3:-1] 
            #print("list_params",list_params)
            
        string_params =""
        for i in range(len(list_params)):
            string_params += list_params[i]+" "
        #print("string_params = ["+str(string_params)+"]")
        #print("bjet_inputs:", string_params, model_type, FOLDER_PATH+DATA_FOLDER, name_stem)
        output = bj_core.main_swig(string_params, model_type, FOLDER_PATH+DATA_FOLDER, name_stem)
        NU_DIM = int(command_params_2[3])
        n_components = 6
        if eic:
            n_components = 10
        full_size = int(NU_DIM*n_components)
        output = (ctypes.c_double * full_size).from_address(int(output))
        output = list(output)
        
        
        #this following part is for debug only
        #I want to be sure that the hdf5 file write out the good proba and parameters
        #It seems there is an issue with the swig version
        #check last sample of the HDF5 file
        # results_directory = "/test_debug"
        # backend_path = BASE_PATH+RESULTS_FOLDER + results_directory + "/backend.h5"
        # if os.path.isfile(backend_path):
        #     print("bjet params:", params_original[1:],'\n')
        #     with h5py.File(backend_path, "r") as f:
        #         a = list(f['mcmc']['chain'])
        #         if a:
        #             for j in range(len(a)-1):
        #                 if sum(sum(a[j])) and not sum(sum(a[j+1])): #a bunch of 10 walkers has begun
        #                     for k in range(len(a[j])-1):
        #                         # if np.log10(params[1]) in a[j][k]:
        #                         #     print("backend K:", a[j][k][0])
        #                         if sum(a[j][k]):
        #                             #if abs(1 - np.log10(params[1])/a[j][k-2][0]) <= 1e-3:
        #                             print("backend params:", a[j][k])       
        
        return output
    else:
        # print("params", params)
        # print(command_params_full)
        result = subprocess.run(
            [BASE_PATH + executable, *command_params_full],
            capture_output=True,
            text=True,
        )

        #this following part is for debug only
        #I want to be sure that the hdf5 file write out the good proba and parameters
        #It seems there is an issue with the swig version
        #check last sample of the HDF5 file
        # results_directory = "/test_debug"
        # backend_path = BASE_PATH+RESULTS_FOLDER + results_directory + "/backend.h5"
        # if os.path.isfile(backend_path):
        #     print("bjet params:", params_original[1:],'\n')
        #     with h5py.File(backend_path, "r") as f:
        #         a = list(f['mcmc']['chain'])
        #         if a:
        #             for j in range(len(a)-1):
        #                 if sum(sum(a[j])) and not sum(sum(a[j+1])): #a bunch of 10 walkers has begun
        #                     for k in range(len(a[j])-1):
        #                         # if np.log10(params[1]) in a[j][k]:
        #                         #     print("backend K:", a[j][k][0])
        #                         if sum(a[j][k]):
        #                             #if abs(1 - np.log10(params[1])/a[j][k-2][0]) <= 1e-3:
        #                             print("backend params:", a[j][k])

        
        #write log file and print
        print(result.stderr)
        print(FOLDER_PATH, folder)
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

    :param parameter_file: Relative path to param file, will be overwritten with parameters; default is <PARAMETER_FOLDER>/params.txt.
    :type parameter_file: str
    :param data_folder: The path to the data folder. Default is None.
    :type data_folder: str
    :param executable: The path to the executable file. Default is None.
    :type executable: str
    :param prev_files: Flag indicating whether to create previous files. Default is False.
    :type prev_files: bool
    :param verbose: Flag indicating whether to display verbose output. Default is False.
    :type verbose: bool
    :return: None
    :rtype: None

    This function generates a system energy dissipation (SED) file by calling a C++ code. The SED file is stored in the <home directory>/DATA_FOLDER with a given file name prefix. Data saved to <home directory>/DATA_FOLDER/<file prefix>_*.dat. Parameter file overwritten or created.

    If the parameter_file is not provided, it will default to the "params.txt" file in the PARAMETER_FOLDER.
    If the data_folder is not provided, it will default to the DATA_FOLDER.
    If the executable is not provided, it will default to the EXECUTABLE file in the CPP_FOLDER.

    If prev_files is True, the C++ code is called with a "2" input mode. If prev_files is False, the C++ code is called with a "3" input mode.

    If verbose is True, the C++ code is called with verbose output. If verbose is False, the C++ code is called without displaying the output.

    .. note:: This method uses older implementation and it is recommended to rewrite the method using bj_core methods, avoiding building a command string.
    """
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
        print(BASE_PATH)
        print([BASE_PATH + executable, *settings_params, BASE_PATH + parameter_file])
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
    """
    Add new data to the current data.

    :raises IOError: If the specified data file cannot be read.
    :raises ValueError: If neither new data nor file suffix is provided.

    :param current_data: Tuple containing current data of the format (current_logv, current_logvFv, current_v, current_vFv).
    :type current_data: tuple
    :param new_data: Optional. Tuple containing new data to be added of the format (model_logv, model_logvFv, model_v, model_vFv).
    :type new_data: tuple
    :param file_suffix: Optional. String representing the file suffix for loading model data.
    :type file_suffix: str
    :param name_stem: Optional. String representing the name stem for the model data file.
    :type name_stem: str
    :param data_folder: Optional. String representing the folder path for the model data file to load.
    :type data_folder: str
    :param cols: Optional. Tuple representing the column indices to extract from the loaded data. Defaults to (0, 2).
    :type cols: tuple
    :return: Tuple containing the interpolated data of the format (logv, logvFv, v, vFv).
    :rtype: tuple
    """
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


def add_SED_component(current_component, new_component):
    """
    Add new model components of the SED

    """
    current_logv, current_logvFv, current_v, current_vFv = current_component
    model_logv, model_logvFv, model_v, model_vFv = new_component
    
    if len(model_logv):
    
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
        
    else:
        #This can happen for extreme  parameters when some emisison components are so weak they are not even computed by Bjet
        logv, logvFv, v, vFv = current_component
    return logv, logvFv, v, vFv


def process_model(
    output, name_stem=None, data_folder=None, verbose=False, eic=False, additional_suffixes=None,
):
    """
    Read a model from data files and returns arrays of frequencies and energy flux.
    This model only uses data from the *compton model* and the *synchrotron model*.
    The data we want are an array of frequencies (v) and an array of corresponding flux energies (vFv). Flux energy is the sum of the flux energy from the synchrotron model and that from the compton model when there are data points for both of them.

    .. warning::
        Requirements:
            - The SED model must have already been created (using `make_SED above`) with the given name_stem.
            - By default, data is read from `DATA_FOLDER/<name_stem>_*.dat`, the location that make_SED writes to.

    Returns:
        - a tuple of 4 1D numpy arrays
        - all arrays have the same length.
        - logv_all, logvFv_all, v_all, vFv_all
        - logv_all: a 1D numpy array of the log of all frequencies
        - logvFv_all: a 1D numpy array of the log of all the energy flux
        - v_all: a 1D numpy array of all frequencies in the data
        - vFv_all: a 1D numpy array of corresponding energy fluxes for the frequencies

    :raises IOError: when data cannot be read

    .. note:: The frequency values used are the ones used in the synchrotron spectrum and the ones in the compton model greater than the max in the synchrotron. The frequency values present in the Compton model in the overlap are not used.

    :param name_stem: The stem name that will be used to construct the file names of the models. If not provided, it will use the default value `NAME_STEM`. Model data is saved in files named in the form <name_stem>_*.dat, specify file name
    :type name_stem: str, optional

    :param data_folder: The folder where the model files are located. If not provided, it will use the default value `DATA_FOLDER`.
    :type data_folder: str, optional

    :param verbose: If True, it runs the full swig wrapper (no I/O). If False it runs the classical version

    :param eic: Determines whether to enable Extended Inverse Compton (EIC) mode. Default is `False`.
    :type eic: bool, optional

    :param additional_suffixes: A list of additional suffixes to be used for reading data in EIC mode. Default is `None`.
    :type additional_suffixes: list, optional

    :return: Four numpy arrays representing the logarithm of frequency, logarithm of energy flux, frequency, and energy flux respectively.
    :rtype: tuple of numpy arrays

    """
    if name_stem is None:
        name_stem = NAME_STEM
    if data_folder is None:
        data_folder = DATA_FOLDER

    if verbose: 
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
                + "_ss.dat cannot be read"
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
        
    else:
        #new method
        NU_DIM = 99
        n_components = 6
        if eic:
            n_components = 8
        start = 0
        stop = NU_DIM
        sliced_output = [0]*n_components
        for i in range(n_components):
            sliced_output[i] = output[start:stop+start]
            start += NU_DIM
            #remove indexes with 0 frequency
            if i%2:
                while 0 in sliced_output[i-1]:
                    tmp = sliced_output[i-1].index(0)
                    del sliced_output[i][tmp]
                    del sliced_output[i-1][tmp]   
                if i == 1:
                    logv = np.array(sliced_output[0])
                    logvFv = np.array(sliced_output[1])
                    v = 10**logv
                    vFv = 10**logvFv
                else:
                    logv_new = np.array(sliced_output[i-1])
                    logvFv_new = np.array(sliced_output[i])
                    v_new = 10**logv_new
                    vFv_new = 10**logvFv_new   
                    logv, logvFv, v, vFv = add_SED_component((logv, logvFv, v, vFv), (logv_new, logvFv_new, v_new, vFv_new))  
    
    # if eic:
    #     if additional_suffixes is None:
    #         additional_suffixes = ["nuc"]
    #     if verbose:
    #         print("process_model EIC mode: ss cs", *additional_suffixes, "read")
    #     for s in additional_suffixes:
    #         logv, logvFv, v, vFv = add_data(
    #             (logv, logvFv, v, vFv),
    #             file_suffix=s,
    #             name_stem=name_stem,
    #             data_folder=data_folder,
    #         )

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
        name_stem (str): data will be saved in files of the form `<name_stem>_*.dat`
        parameter_file (str): optional; relative path of the parameter file
        theta: optional; params is fixed for modeling
        redshift: optional; redshift is fixed for modeling
        min_freq (float): optional; specifies min frequency model SEDs should use; default value used if none provided, which is 5.0e+7.
        max_freq (float): optional; specifies max frequency model SEDs should use; default value used if none provided, which is 1.0e+29.
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
        vFv_all: a 1D numpy array of corresponding energy fluxes for the frequencies

    :raises IOError: if parameter_files cannot be written to (likely, a folder in the path does not exist)

    :param params: The parameters for the model computation. A list of NUM_DIM floats which are the parameters
    :type params: numpy.ndarray
    :param name_stem: The stem name for the model files.
    :type name_stem: str
    :param theta: The angle parameter for the model computation.
    :type theta: float
    :param redshift: The redshift parameter for the model computation.
    :type redshift: float
    :param min_freq: The minimum frequency parameter for the model computation.
    :type min_freq: float
    :param max_freq: The maximum frequency parameter for the model computation.
    :type max_freq: float
    :param torus_temp: The temperature parameter for the torus component in the model computation.
    :type torus_temp: float
    :param torus_luminosity: The luminosity parameter for the torus component in the model computation.
    :type torus_luminosity: float
    :param torus_frac: The fraction parameter for the torus component in the model computation.
    :type torus_frac: float
    :param data_folder: The folder containing the data files.
    :type data_folder: str
    :param executable: The path to the executable file.
    :type executable: str
    :param command_params_full: Additional command line parameters for the model computation.
    :type command_params_full: str
    :param command_params_1: Additional command line parameters for the model computation.
    :type command_params_1: str
    :param command_params_2: Additional command line parameters for the model computation.
    :type command_params_2: str
    :param parameter_file: The path to the parameter file.
    :type parameter_file: str
    :param prev_files: Flag indicating whether to use previous files.
    :type prev_files: bool
    :param use_param_file: Flag indicating whether to use the parameter file.
    :type use_param_file: bool
    :param verbose: Flag indicating whether to display verbose output.
    :type verbose: bool
    :param eic: Flag indicating whether to use EIC mode.
    :type eic: bool
    :param fixed_params: The fixed parameters for the model computation.
    :type fixed_params: numpy.ndarray
    :param folder: The folder for the model files.
    :type folder: str
    :return: The computed model.
    :rtype: numpy.ndarray
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
    output = make_SED(
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
    return process_model(output, name_stem=name_stem, data_folder=data_folder, verbose=verbose, eic=eic)


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
    """

    :param name_stem: The prefix for the file name. If not provided, it defaults to NAME_STEM.
    :type name_stem: str
    :param theta: The angle value. If not provided, it defaults to 0.57.
    :type theta: float
    :param redshift: The redshift value.
    :type redshift: float
    :param min_freq: The minimum frequency value. If not provided, it defaults to 1.0e08.
    :type min_freq: float
    :param max_freq: The maximum frequency value. If not provided, it defaults to 1.0e28.
    :type max_freq: float
    :param data_folder: The folder where the data files are located. If not provided, it defaults to DATA_FOLDER.
    :type data_folder: str
    :param prev_files: Flag indicating whether to use previous files. If set to True, it uses "2" as the value of p, otherwise it uses "3".
    :type prev_files: bool
    :param eic: Flag indicating the model type. If set to True, it uses "1" as the value of model_type, otherwise it uses "0".
    :type eic: bool
    :return: A tuple containing the settings and transformation values, and the constant and numerical values.
    :rtype: tuple, tuple
    """
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
    Given parameters and a file, put the data in the file in the format
    understandable by the C++ code.

    .. note:: Cannot be used with EIC

    :raises IOError: if parameter_files cannot be written to (likely, a folder in the path does not exist)

    :param params: List of parameters
        Form: [delta (linear), K (log), n1 (linear), n2 (linear), gamma_min (log), gamma_max (log), gamma_break (log), B (log), R (log)]
    :type params: list
    :param name_stem: Name stem of the output file. Data will be saved in files named <name_stem>_*.dat
    :type name_stem: str
    :param parameter_file: Output file name. String with the relative path of the parameter file; default is None; will be set to "<PARAMETER_FOLDER>/params.txt"
    :type parameter_file: str
    :param theta: Angle to the line of sight
    :type theta: float
    :param redshift: Redshift
    :type redshift: float
    :param min_freq: Minimal frequency. Specifies min frequency model SEDs should use; default is None; will be set to 5.0e+7
    :type min_freq: float
    :param max_freq: Maximal frequency. Specifies max frequency model SEDs should use; default is None; will be set to 1.0e+29
    :type max_freq: float
    :param verbose: Verbose output flag
    :type verbose: bool
    :return: None parameter_files is over-written (and createed if necessary) with parameter information.
    :rtype: None
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
        min_freq = 1.0e8
    if max_freq is None:
        max_freq = 1.0e28
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
    """
    Converts parameters from logarithmic to linear scale.

    :param params: The parameters to be converted.
    :type params: list[float]
    :param param_is_log: A list indicating whether each parameter is in logarithmic scale. If None, it is determined by modelProperties.
    :type param_is_log: list[bool]
    :param eic: Whether to use EIC for determining logarithmic scale. Default is False.
    :type eic: bool
    :return: A new list of parameters converted from log scale to linear scale if necessary.
    :rtype: list[float]
    """
    if param_is_log is None:
        param_is_log = modelProperties(eic).PARAM_IS_LOG

    ret = params.copy()
    for i in range(len(params)):
        if param_is_log[i]:
            ret[i] = np.power(10, params[i])
    return ret


def params_linear_to_log(params, param_is_log=None, eic=False):
    """
    Converts linear scale parameters to logarithmic scale.

    .. note::
        - If `param_is_log` is not provided, it will be inferred from the model properties based on `eic`.
        - The conversion to log scale is performed using the base 10 logarithm (np.log10).

    :param params: The linear scale parameters to be converted.
    :type params: List or numpy array
    :param param_is_log: Optional. A list indicating whether each parameter should be converted to log scale. If not provided, it will be determined based on the model properties.
    :type param_is_log: List or numpy array, default None
    :param eic: Optional. A flag indicating whether the model properties are based on EIC. Defaults to False. If True, the model properties will be used to determine param_is_log. If False, param_is_log will be used as is.
    :type eic: bool, default False
    :return: The parameters converted to logarithmic scale.
    :rtype: List or numpy array
    """
    if param_is_log is None:
        param_is_log = modelProperties(eic).PARAM_IS_LOG

    ret = params.copy()
    for i in range(len(params)):
        if param_is_log[i]:
            ret[i] = np.log10(params[i])
    return ret
