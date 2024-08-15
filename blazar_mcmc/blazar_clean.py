#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file name: blazar_clean.py

The clean function removes all data files and removes parameter files not named
params.txt.
"""
import os

from blazar_properties import *


def clean(data=True, data_folder=None, parameter_files=True, parameter_folder=None):
    """
    Remove data files and remove parameter files not named params.txt

    Args:
        data (optional): bool
            Specifies if data files should be deleted; default is True
        data_folder (optional): str
            Relative path to data folder; default is None, will be set to DATA_FOLDER
        parameter_files (optional): bool
            Specifies if parameter files should be deleted; default is True
        parameter_folder (optional): str
            Relative path to parameter folder; default is None, will be set to
            PARAMETER_FOLDER

    Returns:
        bool
        success: Returns if both the data_folder and the parameter_folder exists.
    """
    success = True
    if data:
        if data_folder is None:
            data_folder = DATA_FOLDER
        if os.path.exists(BASE_PATH + data_folder):
            files = [
                f for f in os.listdir(BASE_PATH + data_folder) if f.endswith(".dat")
            ]
            for f in files:
                os.remove(BASE_PATH + data_folder + "/" + f)
        else:
            success = False

    if parameter_files:
        if parameter_folder is None:
            parameter_folder = PARAMETER_FOLDER
        if os.path.exists(BASE_PATH + parameter_folder):
            files = [
                f
                for f in os.listdir(BASE_PATH + parameter_folder)
                if f.endswith(".txt")
            ]
            for f in files:
                if f != "params.txt":
                    os.remove(BASE_PATH + parameter_folder + "/" + f)
        else:
            success = False
    return success


if __name__ == "__main__":
    clean()
