#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file name: blazar_clean.py

The clean function removes all data files and removes parameter files not named
params.txt.
"""
import os

from bjet_mcmc.blazar_properties import *

__all__ = ["blazar_clean"]


def clean(data=True, data_folder=None, parameter_files=True, parameter_folder=None):
    """
    Clean function documentation. Remove data files and remove parameter files not named params.txt

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


    :param data: A boolean indicating whether to clean data files. Default is True.
    :type data: bool
    :param data_folder: The folder where the data files are stored. If None, the default data folder will be used.
    :type data_folder: str or None
    :param parameter_files: A boolean indicating whether to clean parameter files. Default is True.
    :type parameter_files: bool
    :param parameter_folder: The folder where the parameter files are stored. If None, the default parameter folder will be used.
    :type parameter_folder: str or None
    :return: A boolean indicating whether the clean operation was successful.
    :rtype: bool

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
