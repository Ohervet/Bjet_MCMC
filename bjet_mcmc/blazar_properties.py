"""
This gets the path to .../blazars-mcmc, which then allows for function calls
to use the absolute path."
"""

import os

# import pathlib
import tempfile
import numpy as np

__all__ = [
    "TMP",
    "INITIALIZE",
    "MAIN_FOLDER",
    "DATA_FOLDER",
    "RESULTS_FOLDER",
    "PARAMETER_FOLDER",
    "CPP_FOLDER",
    "EXECUTABLE",
    "TEMP_DIR",
    "BASE_PATH",
    "FOLDER_PATH",
    "NAME_STEM",
    "BlazarProperties",
    "modelProperties",
    "_get_path",
    "EIC_NUM_DIM",
    "EIC_PARAM_NAMES",
    "EIC_FORMATTED_PARAM_NAMES",
    "EIC_DETAILED_PARAM_NAMES",
    "EIC_PARAM_IS_LOG",
    "SSC_NUM_DIM",
    "SSC_PARAM_NAMES",
    "SSC_FORMATTED_PARAM_NAMES",
    "SSC_DETAILED_PARAM_NAMES",
    "SSC_PARAM_IS_LOG",
]

PROGRAM_NAME = "/bjet_mcmc"
TMP = False
INITIALIZE = False
MAIN_FOLDER = "bjet_mcmc"
DATA_FOLDER = "sed_calculations"
RESULTS_FOLDER = "local_results"
PARAMETER_FOLDER = "parameter_files"
CPP_FOLDER = "bjet_core"
EXECUTABLE = "bj_core"
NAME_STEM = "run"


class BlazarProperties(object):
    """
    Class representing properties of a blazar.

    Attributes:
        NUM_DIM (int): The number of dimensions of the blazar.
        PARAM_NAMES (list): List of parameter names.
        FORMATTED_PARAM_NAMES (list): List of formatted parameter names.
        DETAILED_PARAM_NAMES (list): List of detailed parameter names.
        PARAM_IS_LOG (list): List of boolean values indicating whether each parameter is logged.

    Raises:
        Exception: If the dimensions of the arrays conflict.

    """

    def __init__(
        self,
        num_dim,
        param_names,
        formatted_param_names,
        detailed_param_names,
        param_is_log,
    ):
        self.NUM_DIM = num_dim
        self.PARAM_NAMES = param_names
        self.FORMATTED_PARAM_NAMES = formatted_param_names
        self.DETAILED_PARAM_NAMES = detailed_param_names
        self.PARAM_IS_LOG = param_is_log
        if (
            not self.NUM_DIM
            == len(self.PARAM_NAMES)
            == len(self.FORMATTED_PARAM_NAMES)
            == len(self.DETAILED_PARAM_NAMES)
            == len(self.PARAM_IS_LOG)
        ):
            raise Exception("Dimensions of arrays conflict in blazar_properties.")


def modelProperties(is_eic=False, fixed_params=None):
    """
    Args:
        is_eic (bool): Determines whether to use EIC properties or SSC properties. Defaults to False.
        fixed_params (list): List of fixed parameters. Defaults to None.

    Returns:
        BlazarProperties: Object containing the properties of the selected model.

    """
    sscProperties = BlazarProperties(
        9,
        ["delta", "K", "n_1", "n_2", "gamma_min", "gamma_max", "gamma_break", "B", "R"],
        [
            r"$\delta$",
            r"log $K$",
            r"$n_1$",
            r"$n_2$",
            r"log $\gamma_{min}$",
            r"log $\gamma_{max}$",
            r"log $\gamma_{break}$",
            r"log $B$",
            r"log $R$",
        ],
        [
            "doppler factor (delta)",
            "log K [cm^-3]",
            "alpha 1",
            "alpha 2",
            "log gamma_min",
            "log gamma_max",
            "log gamma_break",
            "log B",
            "log R [cm]",
        ],
        [False, True, False, False, True, True, True, True, True],
    )

    eicProperties = BlazarProperties(
        13,
        [
            "delta",
            "K",
            "n_1",
            "n_2",
            "gamma_min",
            "gamma_max",
            "gamma_break",
            "B",
            "R",
            "bb_temp",
            "l_nuc",
            "tau",
            "blob_dist",
        ],
        [
            r"$\delta$",
            r"log $K$",
            r"$n_1$",
            r"$n_2$",
            r"log $\gamma_{min}$",
            r"log $\gamma_{max}$",
            r"log $\gamma_{break}$",
            r"log $B$",
            r"log $R$",
            r"log $T_{BB}$",
            r"log $L_{nuc}$",
            r"log $\tau$",
            r"log $D_{b}$",
        ],
        [
            "doppler factor (delta)",
            "log K [cm^-3]",
            "alpha 1",
            "alpha 2",
            "log gamma_min",
            "log gamma_max",
            "log gamma_break",
            "log B",
            "log R [cm]",
            "log bb_temp",
            "log l_nuc",
            "log tau",
            "log blob_dist",
        ],
        [
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
        ],
    )

    if is_eic:
        parameters = eicProperties
    else:
        parameters = sscProperties
    # remove any frozen parameter from the parameter list
    if fixed_params:
        fixed_params2 = fixed_params.copy()
        i = 0
        while i < len(fixed_params2):
            if fixed_params2[i] != -np.inf:
                parameters.NUM_DIM -= 1
                del parameters.PARAM_NAMES[i]
                del parameters.FORMATTED_PARAM_NAMES[i]
                del parameters.DETAILED_PARAM_NAMES[i]
                del parameters.PARAM_IS_LOG[i]
                del fixed_params2[i]
            else:
                i += 1
    return parameters


# TODO Implement a different scheme for saving output files. Do not place in program install location but in a user
#  defined analysis directory.
def _get_path():
    """

    This method, '_get_path()', returns the base path of the file containing the method.

    Parameters:
        None.

    Returns:
        str: The base path of the file.

    Raises:
        Exception: If the PROGRAM_NAME is not found in the file path.

    """
    base_path = str(os.path.dirname(__file__))
    base_path = f"{base_path[: base_path.find('Bjet_MCMC')]}Bjet_MCMC/bjet_mcmc"
    # base_path = str(pathlib.Path().resolve())
    stop = base_path.find(PROGRAM_NAME)
    if stop == -1:
        raise Exception(PROGRAM_NAME + " is not in file path")
    return base_path[: -len(MAIN_FOLDER)] + "/"


if TMP:
    TEMP_DIR = tempfile.mkdtemp()
    BASE_PATH = TEMP_DIR + "/"
    FOLDER_PATH = _get_path()
    print(BASE_PATH)
else:
    TEMP_DIR = None
    BASE_PATH = _get_path()
    FOLDER_PATH = _get_path()

# with EIC
EIC_NUM_DIM = 13
EIC_PARAM_NAMES = [
    "delta",
    "K",
    "n_1",
    "n_2",
    "gamma_min",
    "gamma_max",
    "gamma_break",
    "B",
    "R",
    "bb_temp",
    "L_nuc",
    "tau",
    "blob_dist",
]
EIC_FORMATTED_PARAM_NAMES = [
    r"$\delta$",
    r"log $K$",
    r"$n_1$",
    r"$n_2$",
    r"log $\gamma_{min}$",
    r"log $\gamma_{max}$",
    r"log $\gamma_{break}$",
    r"log $B$",
    r"log $R$",
    r"log $T_{BB}$",
    r"log $L_{nuc}$",
    r"$\tau$",
    r"log $D_b$",
]
EIC_DETAILED_PARAM_NAMES = [
    "doppler factor (delta)",
    "log K [cm^-3]",
    "alpha 1",
    "alpha 2",
    "log gamma_min",
    "log gamma_max",
    "log gamma_break",
    "log B",
    "log R [cm]",
    "log bb_temp",
    "log L_nuc",
    "log tau",
    "log blob_dist",
]
EIC_PARAM_IS_LOG = [
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

# SSC
SSC_NUM_DIM = 9
SSC_PARAM_NAMES = [
    "delta",
    "K",
    "n_1",
    "n_2",
    "gamma_min",
    "gamma_max",
    "gamma_break",
    "B",
    "R",
]
SSC_FORMATTED_PARAM_NAMES = [
    r"$\delta$",
    r"log $K$",
    r"$n_1$",
    r"$n_2$",
    r"log $\gamma_{min}$",
    r"log $\gamma_{max}$",
    r"log $\gamma_{break}$",
    r"log $B$",
    r"log $R$",
]
SSC_DETAILED_PARAM_NAMES = [
    "doppler factor (delta)",
    "log K [cm^-3]",
    "alpha 1",
    "alpha 2",
    "log gamma_min",
    "log gamma_max",
    "log gamma_break",
    "log B",
    "log R [cm]",
]
SSC_PARAM_IS_LOG = [False, True, False, False, True, True, True, True, True]
