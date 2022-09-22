#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file name: blazar_initialize.py

Purpose: Implement functions to make necessary directories and compile the c++ code

Contains:
FUNCTIONS:
make_dirs
compile_bjet
initialize
-----------------
make_dirs(data_folder=None, results_folder=None, parameter_folder=None, parameter_file=False)
    Creates data folder, results folder, and parameter folder. Folders can be
    specified. Others, folders specified in blazar_properties.py will be used.
compile_bjet(bjet_folder=None, executable=None, verbose=False)
    Compile c++ code
initialize(data_folder=None, results_folder=None, parameter_folder=None, parameter_file=False,
               bjet_folder=None, executable=None, run_compile=True)
    Calls make_dirs and compile_bjet. If code is to be run in the temp directory, the executable
    of the C++ code is copied into the temp directory.
"""
import os
import shutil
import subprocess

from blazar_properties import *


def make_dirs(data_folder=None, results_folder=None, parameter_folder=None, parameter_file=False):
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
            os.mkdir(FOLDER_PATH + folder)  # want results directory in both tmp and main if tmp


def compile_bjet(bjet_folder=None, executable=None, verbose=False):
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
        subprocess.run("make", stderr=open(os.devnull, 'wb'),
                       stdout=open(os.devnull, 'wb'))


def initialize(data_folder=None, results_folder=None, parameter_folder=None, parameter_file=False,
               bjet_folder=None, executable=None, run_compile=True):
    make_dirs(data_folder=data_folder, results_folder=results_folder, parameter_folder=parameter_folder,
              parameter_file=parameter_file)
    if run_compile:
        compile_bjet(bjet_folder=bjet_folder, executable=executable)
    if TMP:
        if not os.path.exists(BASE_PATH + bjet_folder):
            os.mkdir(BASE_PATH + bjet_folder)
            shutil.copy(FOLDER_PATH + bjet_folder + "/" + executable, BASE_PATH + bjet_folder + "/" + executable)


if __name__ == '__main__':
    initialize()
