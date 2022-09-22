"""

"""

import glob
import os
import random
import sys

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

residual = False
residual2 = False
backend_file = input("Enter relative path to backend file: ")
eic = input("Enter eic (T/True, F/False): ")
block_number = int(input("Enter block number: "))
if block_number not in (1, 2, 6):
    print("Invalid block number.")
    sys.exit(1)
eic = eic[0].strip().lower()[0] == 't'
discard = input("Enter number of steps to discard: ")
discard = int(discard)

folder = backend_file[:backend_file.rfind('/')]
param_min_vals, param_max_vals = blazar_utils.min_max_parameters(eic=eic)

reader = emcee.backends.HDFBackend(FOLDER_PATH + backend_file)

chain = reader.get_chain(discard=discard)
flat_samples = reader.get_chain(flat=True, discard=discard)
log_probs = reader.get_log_prob(discard=discard)
flat_log_probs = reader.get_log_prob(flat=True, discard=discard)

best_log_prob, best_params = blazar_report.get_best_log_prob_and_params(
    log_probs=flat_log_probs, flat_chain=flat_samples)

min_1sigma_params, max_1sigma_params = blazar_report.min_max_params_1sigma(flat_samples, flat_log_probs, eic=eic)
indices_within_1sigma = blazar_report.get_indices_within_1sigma(flat_log_probs, eic=eic)

name_stem = "stem_b6_no_eic"


if os.path.exists(FOLDER_PATH + folder + "/info.txt"):
    info = blazar_report.parse_info_doc(folder + "/info.txt")
    configs = info["configs"]
    print(configs)
else:
    configs = blazar_utils.read_configs()
data = blazar_utils.read_data(configs["data_file"])
v_data, vFv_data, err_data = data

command_params_1, command_params_2 = blazar_model.command_line_sub_strings(name_stem=name_stem, redshift=configs["redshift"], prev_files=False, eic=eic)
command_params_2[3] = "300"  # number of points used to make SED
#command_params_1 = None
#command_params_2 = None

best_model = blazar_model.make_model(best_params, name_stem=name_stem, redshift=.340, command_params_1=command_params_1,
                                     command_params_2=command_params_2, eic=eic)
swift_uvot = np.where(data[0] < np.power(10, 16))[0]
swift_uvot_v = data[0][swift_uvot]
swift_uvot_vFv = data[1][swift_uvot]
swift_uvot_err = data[2][swift_uvot]

swift_xrt = np.where(np.logical_and(data[0] > np.power(10, 17), data[0] < np.power(10, 20)))[0]
swift_xrt_v = data[0][swift_xrt]
swift_xrt_vFv = data[1][swift_xrt]
swift_xrt_err = data[2][swift_xrt]

# block 1, 2
if block_number == 1 or block_number == 2:
    fermi = np.where(np.logical_and(data[0] > 10 ** 20, data[0] < 2 * 10 ** 25))
    veritas = np.where(data[0] > 2 * 10 ** 25)
elif block_number == 6:
    fermi = np.where(np.logical_and(data[0] > 10 ** 20, data[1] > 10 ** -11.1))
    veritas = np.where(np.logical_and(data[0] > 10 ** 20, data[1] < 10 ** -11.1))
fermi_v = data[0][fermi]
fermi_vFv = data[1][fermi]
fermi_err = data[2][fermi]

veritas_v = data[0][veritas]
veritas_vFv = data[1][veritas]
veritas_err = data[2][veritas]


min_per_param, max_per_param = blazar_plots.get_params_1sigma_ranges(flat_samples, indices_within_1sigma,
                                                                     eic=eic)
to_plot = np.concatenate((np.array(min_per_param), np.array(max_per_param)))

logv = best_model[0].copy()
lowest_per_point, highest_per_point = blazar_plots.get_min_max_per_point(logv, to_plot, name_stem=name_stem,
                                                                         redshift=.340,
                                                                         command_params_1=command_params_1,
                                                                         command_params_2=command_params_2, eic=eic)

if residual or residual2:
    fig, (ax, ax1) = plt.subplots(2, sharex='all', gridspec_kw={'height_ratios': [6, 1]})
else:
    fig, ax = plt.subplots()
plt.xlabel(r"Frequency $\nu$ (Hz)")
ax.set_ylabel(r"Energy Flux $\nu F_{\nu}$ (erg cm$^{-2}$ s$^{-1})$")

ax.set_yscale("log")
plt.xscale("log")

ax.fill_between(np.power(10, logv), np.power(10, lowest_per_point), np.power(10, highest_per_point), alpha=.5, label=r"Within 1$\sigma$")

ax.plot(best_model[2], best_model[3], label="Best model")

blazar_model.make_model(best_params, name_stem=name_stem, redshift=.340, command_params_1=command_params_1,
                        command_params_2=command_params_2, eic=eic)

synchrotron_model = np.loadtxt(FOLDER_PATH + "sed_calculations/" + name_stem + "_ss.dat", delimiter=' ')
logv_synchrotron = synchrotron_model[:, 0]
logvFv_synchrotron = synchrotron_model[:, 2]
v_synchrotron = np.power(10, logv_synchrotron)
vFv_synchrotron = np.power(10, logvFv_synchrotron)

compton_model = np.loadtxt(FOLDER_PATH + "sed_calculations/" + name_stem + "_cs.dat", delimiter=' ')
logv_compton = compton_model[:, 0]
logvFv_compton = compton_model[:, 2]
v_compton = np.power(10, logv_compton)
vFv_compton = np.power(10, logvFv_compton)

ax.plot(v_synchrotron, vFv_synchrotron, 'k--', label="Synchrotron", alpha=.5)
ax.plot(v_compton, vFv_compton, 'k-', label="Self Compton", alpha=.5)

if eic:
   cs2 = np.loadtxt(FOLDER_PATH + "sed_calculations/" + name_stem + "_cs2.dat", delimiter=' ')
   logv_cs2 = cs2[:, 0]
   logvFv_cs2 = cs2[:, 2]
   v_cs2 = np.power(10, logv_cs2)
   vFv_cs2 = np.power(10, logvFv_cs2)

   ecs = np.loadtxt(FOLDER_PATH + "sed_calculations/" + name_stem + "_ecs.dat", delimiter=' ')
   logv_ecs = ecs[:, 0]
   logvFv_ecs = ecs[:, 2]
   v_ecs = np.power(10, logv_ecs)
   vFv_ecs = np.power(10, logvFv_ecs)

   nuc = np.loadtxt(FOLDER_PATH + "sed_calculations/" + name_stem + "_nuc.dat", delimiter=' ')
   logv_nuc = nuc[:, 0]
   logvFv_nuc = nuc[:, 2]
   v_nuc = np.power(10, logv_nuc)
   vFv_nuc = np.power(10, logvFv_nuc)

   ax.plot(v_cs2, vFv_cs2, '-', label="2nd Order SC")
   ax.plot(v_ecs, vFv_ecs, '-', label="EIC")
   ax.plot(v_nuc, vFv_nuc, '-', label="nucleus")

l_swift_uvot = ax.errorbar(swift_uvot_v, swift_uvot_vFv, yerr=(swift_uvot_err, swift_uvot_err), fmt='bs', label="Swift-UVOT", markersize=3, elinewidth=1)
l_swift_xrt = ax.errorbar(swift_xrt_v, swift_xrt_vFv, yerr=(swift_xrt_err, swift_xrt_err), fmt='gD', label="Swift-XRT", markersize=3, elinewidth=1)
l_fermi = ax.errorbar(fermi_v, fermi_vFv, yerr=(fermi_err, fermi_err), fmt='r^', label="Fermi-LAT", markersize=3, elinewidth=1)
l_veritas = ax.errorbar(veritas_v, veritas_vFv, yerr=(veritas_err, veritas_err), fmt='mv', label="VERITAS", markersize=3, elinewidth=1)
#ax.errorbar(data[0], data[1], yerr=(data[2], data[2]), fmt='b.', ecolor='k', label="Data points")



#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')
h = 4.135667662E-15

def v_to_e(val):
    return val * h

def e_to_v(val):
    return val / h

secax = ax.secondary_xaxis('top', functions=(v_to_e, e_to_v))
secax.set_xlabel("Energy (eV)")

if residual:
    f_best = interpolate.interp1d(best_model[2], best_model[3])
    f_low = interpolate.interp1d(best_model[2], np.power(10, lowest_per_point))
    f_high = interpolate.interp1d(best_model[2], np.power(10, highest_per_point))
    ax1.fill_between(data[0], data[1] - f_low(data[0]), data[1] - f_high(data[0]), alpha=.5,
                     zorder=1)

    ax1.errorbar(swift_uvot_v, swift_uvot_vFv - f_best(swift_uvot_v), yerr=(swift_uvot_err, swift_uvot_err), fmt='bs', label="Swift-UVOT", markersize=2, elinewidth=.75)
    ax1.errorbar(swift_xrt_v, swift_xrt_vFv - f_best(swift_xrt_v), yerr=(swift_xrt_err, swift_xrt_err), fmt='gD', label="Swift-XRT", markersize=2, elinewidth=.75)
    ax1.errorbar(fermi_v, fermi_vFv - f_best(fermi_v), yerr=(fermi_err, fermi_err), fmt='r^', label="Fermi-LAT", markersize=2, elinewidth=.75)
    ax1.errorbar(veritas_v, veritas_vFv - f_best(veritas_v), yerr=(veritas_err, veritas_err), fmt='mv', label="VERITAS", markersize=2, elinewidth=.75)

    ax1.axhline(0, zorder=0, color='k', alpha=.5)
    ax1.set_ylabel("Residuals")

if residual2:
    f_best = interpolate.interp1d(best_model[2], best_model[3])
    f_low = interpolate.interp1d(best_model[2], np.power(10, lowest_per_point))
    f_high = interpolate.interp1d(best_model[2], np.power(10, highest_per_point))
    ax1.fill_between(data[0], (data[1] - f_low(data[0])) / data[1], (data[1] - f_high(data[0]))/data[1], alpha=.5,
                     zorder=1)

    ax1.errorbar(swift_uvot_v, swift_uvot_vFv - f_best(swift_uvot_v), yerr=(swift_uvot_err, swift_uvot_err), fmt='bs', label="Swift-UVOT", markersize=2, elinewidth=.75)
    ax1.errorbar(swift_xrt_v, (swift_xrt_vFv - f_best(swift_xrt_v)) / swift_xrt_vFv, yerr=(swift_xrt_err, swift_xrt_err), fmt='gD', label="Swift-XRT", markersize=2, elinewidth=.75)
    ax1.errorbar(fermi_v, (fermi_vFv - f_best(fermi_v)) / fermi_vFv, yerr=(fermi_err, fermi_err), fmt='r^', label="Fermi-LAT", markersize=2, elinewidth=.75)
    ax1.errorbar(veritas_v, (veritas_vFv - f_best(veritas_v)) / veritas_vFv, yerr=(veritas_err, veritas_err), fmt='mv', label="VERITAS", markersize=2, elinewidth=.75)

    ax1.axhline(0, zorder=0, color='k', alpha=.5)
    ax1.set_ylabel("Residuals")

new_min, new_max = blazar_plots.scale_to_values(data[1], lower_adjust_multiplier=5000, upper_adjust_multiplier=3)
ax.set_ylim(new_min, new_max)

plt.tight_layout()


plt.savefig(BASE_PATH + folder + "/plot_with_data_sources.pdf")

plt.show()
#blazar_report.save(folder, redshift=configs["redshift"], eic=eic)

#for f in glob.glob(FOLDER_PATH + DATA_FOLDER + "/" + name_stem + "_*"):
#    os.remove(f)

blazar_report.save_plots_and_info(configs, (v_data, vFv_data, err_data), param_min_vals, param_max_vals, folder=folder, samples=(chain, flat_samples, log_probs, flat_log_probs), use_samples=True, redshift=0.340, eic=True, verbose=True)
"""
np.array([36.42218375,  1.14982865,  2.06910652,  4.031,  3.11270207,
        6.97329194,  4.64969522, -2.89173394, 18.28945182,  4.50994885,
       38.8       , -5.4       , 15.3       ])
       Combine 2 data sources

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

residual = False
residual2 = False

backend_1 = "local_results/3C66A_b6_eic_2022-06-08-20:17:26/backend.h5"
backend_2 = "local_results/3C66A_b6_eic_cont_2022-06-09-17:49:20/backend.h5"
eic = True
discard = 500
folder1 = backend_1[:backend_1.rfind('/')]
folder2 = backend_2[:backend_2.rfind('/')]

param_min_vals, param_max_vals = blazar_utils.min_max_parameters(eic=eic)

name_stem = "stem"
reader1 = emcee.backends.HDFBackend(BASE_PATH + backend_1)
reader2 = emcee.backends.HDFBackend(BASE_PATH + backend_2)
chain1 = reader1.get_chain(discard=500)
flat_samples1 = reader1.get_chain(flat=True, discard=500)
log_probs1 = reader1.get_log_prob(discard=500)
flat_log_probs1 = reader1.get_log_prob(flat=True, discard=500)
chain2 = reader2.get_chain()
flat_samples2 = reader2.get_chain(flat=True)
log_probs2 = reader2.get_log_prob()
flat_log_probs2 = reader2.get_log_prob(flat=True)
chain = np.concatenate((chain1, chain2))
flat_chain = np.concatenate((flat_samples1, flat_samples2))
flat_samples = flat_chain
log_probs = np.concatenate((log_probs1, log_probs2))
flat_log_probs = np.concatenate((flat_log_probs1, flat_log_probs2))

best_log_prob, best_params = blazar_report.get_best_log_prob_and_params(
    log_probs=flat_log_probs, flat_chain=flat_samples)

info = blazar_report.parse_info_doc(folder1 + "/info.txt")
configs = info["configs"]
data = blazar_utils.read_data(configs["data_file"])
v_data, vFv_data, err_data = data

command_params_1 = None
command_params_2 = None

best_model = blazar_model.make_model(best_params, name_stem=name_stem, command_params_1=command_params_1, command_params_2=command_params_2, redshift=.340, eic=eic)

min_1sigma_params, max_1sigma_params = blazar_report.min_max_params_1sigma(flat_samples, flat_log_probs, eic=eic)
indices_within_1sigma = blazar_report.get_indices_within_1sigma(flat_log_probs, eic=eic)


backend_1 = "local_results/3C66A_b6_eic_2022-06-08-20:17:26/backend.h5"
backend_2 = "local_results/3C66A_b6_eic_cont_2022-06-09-17:49:20/backend.h5"
"""