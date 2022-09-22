"""

"""

import glob
import os
import random

import emcee
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy import interpolate
import corner


import blazar_model
import blazar_report
import blazar_utils
import blazar_run_mcmc
import blazar_plots
import blazar_initialize
import blazar_clean
from blazar_properties import *

backend_file = input("Enter relative path to backend file: ")
eic = input("Enter eic (T/True, F/False): ")
eic = eic[0].strip().lower()[0] == 't'
discard = input("Enter number of steps to discard: ")
discard = int(discard)

folder = backend_file[:backend_file.rfind('/')]
param_min_vals, param_max_vals = blazar_utils.min_max_parameters(eic=eic)

reader = emcee.backends.HDFBackend(BASE_PATH + backend_file)

chain = reader.get_chain(discard=discard)
flat_samples = reader.get_chain(flat=True, discard=discard)
log_probs = reader.get_log_prob(discard=discard)
flat_log_probs = reader.get_log_prob(flat=True, discard=discard)

best_log_prob, best_params = blazar_report.get_best_log_prob_and_params(
    log_probs=flat_log_probs, flat_chain=flat_samples)

min_1sigma_params, max_1sigma_params = blazar_report.min_max_params_1sigma(flat_samples, flat_log_probs, eic=eic)
indices_within_1sigma = blazar_report.get_indices_within_1sigma(flat_log_probs, eic=eic)

#name_stem = NAME_STEM + str(random.getrandbits(20))
name_stem = "a"
best_model = blazar_model.make_model(best_params, name_stem=name_stem, redshift=0.340, use_param_file=False, eic=eic)
#for f in glob.glob(BASE_PATH + DATA_FOLDER + "/" + name_stem + "_*"):
 #   os.remove(f)

if os.path.exists(FOLDER_PATH + folder + "/info.txt"):
    info = blazar_report.parse_info_doc(folder + "/info.txt")
    configs = info["configs"]
else:
    configs = blazar_utils.read_configs()
data = blazar_utils.read_data(configs["data_file"])
v_data, vFv_data, err_data = data
