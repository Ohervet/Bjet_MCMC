#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file name: blazar_initialize.py

Purpose: Implement functions to make necessary directories and compile the c++ code

"""
import os
import shutil
import subprocess

from bjet_mcmc.blazar_properties import *

__all__ = ["make_dirs", "compile_bjet", "initialize"]


def make_dirs(
    data_folder=None, results_folder=None, parameter_folder=None, parameter_file=False
):
    """
    This function, `make_dirs`, creates folders based on the input parameters.

    The function first checks if any of the input folder parameters are None. If they are, it assigns default values based on predefined constants.

    Then, it iterates through each of the folders (data_folder, results_folder, parameter_folder) and checks if they are not None. If they are not None, it checks if the folder path does not already exist. If it doesn't exist, it creates the folder.

    Finally, it checks if the current folder being processed is the results_folder and if the folder path doesn't already exist. If it doesn't exist, it creates the folder under a specific path.

    Note: This function relies on the presence of certain constants (e.g. DATA_FOLDER, RESULTS_FOLDER, PARAMETER_FOLDER, BASE_PATH, FOLDER_PATH), but these constants are not defined or included in this documentation.

    Example usage:
        make_dirs(data_folder='data', results_folder='results', parameter_folder='parameters', parameter_file=True)

    :param data_folder: The folder to store data files. Defaults to DATA_FOLDER.
    :type data_folder: str
    :param results_folder: The folder to store results files. Defaults to RESULTS_FOLDER.
    :type results_folder: str
    :param parameter_folder: The folder to store parameter files.  If None and parameter_file is True, it defaults to the value of PARAMETER_FOLDER. If None and parameter_file is False, the parameter folder will not be created.
    :type parameter_folder: str
    :param parameter_file: A boolean indicating whether a parameter file is present or not. Defaults to False.
    :type parameter_file: bool
    :return: None
    :rtype: None


    """
    if data_folder is None:
        data_folder = DATA_FOLDER
    if results_folder is None:
        results_folder = RESULTS_FOLDER
    if not parameter_file:
        parameter_folder = None  # no parameter files --> don't make folder
    elif parameter_folder is None:
        parameter_folder = PARAMETER_FOLDER
    for folder in (data_folder, results_folder, parameter_folder):
        if folder is not None:
            if not os.path.exists(BASE_PATH + folder):
                os.mkdir(BASE_PATH + folder)
        if folder == results_folder and not os.path.exists(FOLDER_PATH + folder):
            os.mkdir(
                FOLDER_PATH + folder
            )  # want results directory in both tmp and main if tmp


def compile_bjet(bjet_folder=None, executable=None, verbose=False):
    """
    Compiles the specified BJET folder using the specified executable.

    :param bjet_folder: The folder containing the BJET files to be compiled. If not provided, the default CPP_FOLDER will be used.
    :type bjet_folder: str
    :param executable: The name of the executable file to be generated. If not provided, the default EXECUTABLE will be used.
    :type executable: str
    :param verbose: If True, the make process will display verbose output. If False, the output will be suppressed. Default is False.
    :type verbose: bool
    :return: None
    :rtype: None
    """
    if bjet_folder is None:
        bjet_folder = CPP_FOLDER
    if executable is None:
        executable = EXECUTABLE

    os.chdir(FOLDER_PATH + bjet_folder)
    files = [f for f in os.listdir() if f.endswith(".o")]
    for f in files:
        os.remove(f)
    if os.path.exists(executable):
        os.remove(executable)
    if verbose:
        subprocess.run("make", executable)
    else:
        subprocess.run(
            "make", stderr=open(os.devnull, "wb"), stdout=open(os.devnull, "wb")
        )


def initialize(
    data_folder=None,
    results_folder=None,
    parameter_folder=None,
    parameter_file=False,
    bjet_folder=None,
    executable=None,
    run_compile=True,
):
    """
    Initializes the software by creating necessary directories, compiling the required files, and copying the executable.

    :param data_folder: The path to the data folder. Defaults to None.
    :type data_folder: str
    :param results_folder: The path to the results folder. Defaults to None.
    :type results_folder: str
    :param parameter_folder: The path to the parameter folder. Defaults to None.
    :type parameter_folder: str
    :param parameter_file: Indicates whether a parameter file is present. Defaults to False.
    :type parameter_file: bool
    :param bjet_folder: The path to the bjet folder. Defaults to None.
    :type bjet_folder: str
    :param executable: The name of the executable file. Defaults to None.
    :type executable: str
    :param run_compile: Indicates whether to compile the bjet executable. Defaults to True.
    :type run_compile: bool
    :return: None
    :rtype: None
    """
    make_dirs(
        data_folder=data_folder,
        results_folder=results_folder,
        parameter_folder=parameter_folder,
        parameter_file=parameter_file,
    )
    if run_compile:
        compile_bjet(bjet_folder=bjet_folder, executable=executable)
    if TMP:
        if not os.path.exists(BASE_PATH + bjet_folder):
            os.mkdir(BASE_PATH + bjet_folder)
            shutil.copy(
                FOLDER_PATH + bjet_folder + "/" + executable,
                BASE_PATH + bjet_folder + "/" + executable,
            )


if __name__ == "__main__":
    initialize()
