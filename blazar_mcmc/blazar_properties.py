"""
This gets the path to .../blazars-mcmc, which then allows for function calls
to use the absolute path."
"""

import pathlib
import tempfile
import numpy as np

PROGRAM_NAME = "/Bjet_MCMC"
TMP = False
INITIALIZE = False
MAIN_FOLDER = "blazars-mcmc"
DATA_FOLDER = "sed_calculations"
RESULTS_FOLDER = "local_results"
PARAMETER_FOLDER = "parameter_files"
CPP_FOLDER = "bjet_core"
EXECUTABLE = "bj_core"
NAME_STEM = "run"


class BlazarProperties(object):
    def __init__(self, num_dim, param_names, formatted_param_names, detailed_param_names, param_is_log):
        self.NUM_DIM = num_dim
        self.PARAM_NAMES = param_names
        self.FORMATTED_PARAM_NAMES = formatted_param_names
        self.DETAILED_PARAM_NAMES = detailed_param_names
        self.PARAM_IS_LOG = param_is_log
        if not self.NUM_DIM == len(self.PARAM_NAMES) == len(self.FORMATTED_PARAM_NAMES) == len(
                self.DETAILED_PARAM_NAMES) == len(
                self.PARAM_IS_LOG):
            raise Exception("Dimensions of arrays conflict in blazar_properties.")


def modelProperties(is_eic=False, fixed_params=None):
    sscProperties = BlazarProperties(9, ["delta", "K", "n_1", "n_2", "gamma_min", "gamma_max", "gamma_break", "B",
                                         "R"],
                                     [r"$\delta$", r"log $K$", r"$n_1$", r"$n_2$", r"log $\gamma_{min}$",
                                      r"log $\gamma_{max}$", r"log $\gamma_{break}$", r"log $B$", r"log $R$"],
                                     ["doppler factor (delta)", "log K [cm^-3]", "alpha 1", "alpha 2", "log gamma_min",
                                      "log gamma_max", "log gamma_break", "log B", "log R [cm]"],
                                     [False, True, False, False, True, True, True, True, True])
    
    eicProperties = BlazarProperties(13,
                                     ["delta", "K", "n_1", "n_2", "gamma_min", "gamma_max", "gamma_break", "B", "R",
                                      "bb_temp", "l_nuc", "tau", "blob_dist"],
                                     [r"$\delta$", r"log $K$", r"$n_1$", r"$n_2$", r"log $\gamma_{min}$",
                                      r"log $\gamma_{max}$", r"log $\gamma_{break}$", r"log $B$", r"log $R$", r"log $T_{BB}$",
                                      r"log $L_{nuc}$", r"log $\tau$", r"log $D_{b}$"],
                                     ["doppler factor (delta)", "log K [cm^-3]", "alpha 1", "alpha 2", "log gamma_min",
                                      "log gamma_max", "log gamma_break", "log B", "log R [cm]", "log bb_temp", "log l_nuc",
                                      "log tau", "log blob_dist"],
                                     [False, True, False, False, True, True, True, True, True, True, True, True, True]
                                     )    
    
    if is_eic:
        parameters =  eicProperties
    else:
        parameters =  sscProperties
    #remove any frozen parameter from the parameter list
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
              i+=1
    return parameters
    

def _get_path():
    base_path = str(pathlib.Path().resolve())
    stop = base_path.find(PROGRAM_NAME)
    if stop == -1:
        raise Exception(PROGRAM_NAME + " is not in file path")
    return base_path[:-len(MAIN_FOLDER)] + "/"


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
EIC_PARAM_NAMES = ["delta", "K", "n_1", "n_2", "gamma_min", "gamma_max", "gamma_break", "B", "R", "bb_temp",
                   "L_nuc", "tau", "blob_dist"]
EIC_FORMATTED_PARAM_NAMES = [r"$\delta$", r"log $K$", r"$n_1$", r"$n_2$", r"log $\gamma_{min}$",
                             r"log $\gamma_{max}$", r"log $\gamma_{break}$", r"log $B$", r"log $R$", r"log $T_{BB}$",
                             r"log $L_{nuc}$", r"$\tau$", r"log $D_b$"]
EIC_DETAILED_PARAM_NAMES = ["doppler factor (delta)", "log K [cm^-3]", "alpha 1", "alpha 2", "log gamma_min",
                            "log gamma_max", "log gamma_break", "log B", "log R [cm]", "log bb_temp", "log L_nuc",
                            "log tau", "log blob_dist"]
EIC_PARAM_IS_LOG = [False, True, False, False, True, True, True, True, True, True, True, True, True]

# SSC
SSC_NUM_DIM = 9
SSC_PARAM_NAMES = ["delta", "K", "n_1", "n_2", "gamma_min", "gamma_max", "gamma_break", "B", "R"]
SSC_FORMATTED_PARAM_NAMES = [r"$\delta$", r"log $K$", r"$n_1$", r"$n_2$", r"log $\gamma_{min}$",
                             r"log $\gamma_{max}$", r"log $\gamma_{break}$", r"log $B$", r"log $R$"]
SSC_DETAILED_PARAM_NAMES = ["doppler factor (delta)", "log K [cm^-3]", "alpha 1", "alpha 2", "log gamma_min",
                            "log gamma_max", "log gamma_break", "log B", "log R [cm]"]
SSC_PARAM_IS_LOG = [False, True, False, False, True, True, True, True, True]
