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
cmap = plt.get_cmap("tab10")

import blazar_model
import blazar_report
import blazar_utils
import blazar_run_mcmc
import blazar_plots
import blazar_initialize
import blazar_clean
from blazar_properties import *


residual = True
backend_file = input("Enter relative path to backend file: ")
#example: local_results/_2023-01-07-10:20:38/backend.h5

folder = backend_file[:backend_file.rfind('/')]

if os.path.exists(FOLDER_PATH + folder + "/info.txt"):
    info = blazar_report.parse_info_doc(folder + "/info.txt")
    configs = info["configs"]
    print(configs)
else:
    configs = blazar_utils.read_configs()
data = blazar_utils.read_data(configs["data_file"], instrument=True)
v_data, vFv_data, err_data_down, err_data_up, instrument_data = data
redshift = configs["redshift"]
discard = configs["discard"]
eic = configs["eic"]


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

name_stem = "plot"
command_params_1, command_params_2 = blazar_model.command_line_sub_strings(name_stem=name_stem, redshift=redshift, prev_files=False, eic=eic)
command_params_2[3] = "300"  # number of points used to make SED
#command_params_1 = None
#command_params_2 = None

best_model = blazar_model.make_model(best_params, name_stem=name_stem, redshift=redshift, command_params_1=command_params_1,
                                     command_params_2=command_params_2, eic=eic)
min_per_param, max_per_param = blazar_plots.get_params_1sigma_ranges(flat_samples, indices_within_1sigma,eic=eic)
to_plot = np.concatenate((np.array(min_per_param), np.array(max_per_param)))
logv = best_model[0].copy()
lowest_per_point, highest_per_point = blazar_plots.get_min_max_per_point(logv, to_plot, name_stem=name_stem,
                                                                         redshift=redshift,
                                                                         command_params_1=command_params_1,
                                                                         command_params_2=command_params_2, eic=eic)

if residual:
    fig, (ax, ax1) = plt.subplots(2, sharex='all', gridspec_kw={'height_ratios': [6, 1.5]})
    #plt.subplots_adjust(hspace=0)
else:
    fig, ax = plt.subplots()
plt.xlabel(r"Frequency $\nu$ (Hz)")
ax.set_ylabel(r"Energy Flux $\nu F_{\nu}$ (erg cm$^{-2}$ s$^{-1})$")

ax.set_yscale("log")
plt.xscale("log")

ax.fill_between(np.power(10, logv), np.power(10, lowest_per_point), np.power(10, highest_per_point), alpha=.5, label=r"Within 1$\sigma$")

ax.plot(best_model[2], best_model[3], label="Best model")

blazar_model.make_model(best_params, name_stem=name_stem, redshift=redshift, command_params_1=command_params_1,
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
   
   
list_intruments=[instrument_data[0]]
v_data_inst = [v_data[0]]
vFv_data_inst = [vFv_data[0]]
err_data_inst_down = [err_data_down[0]]
err_data_inst_up = [err_data_up[0]]
for i in range(1,len(instrument_data)):
    if instrument_data[i] != list_intruments[-1]:
        ax.errorbar(v_data_inst, vFv_data_inst, yerr=(err_data_inst_down, err_data_inst_up), fmt='o', label=str(list_intruments[-1]), 
                    markersize=3, elinewidth=1, color = cmap(len(list_intruments)))
        list_intruments.append(instrument_data[i])
        v_data_inst = [v_data[i]]
        vFv_data_inst = [vFv_data[i]]
        err_data_inst_down = [err_data_down[i]]
        err_data_inst_up = [err_data_up[i]]
    else:   
        v_data_inst.append(v_data[i])
        vFv_data_inst.append(vFv_data[i])
        err_data_inst_down.append(err_data_down[i])
        err_data_inst_up.append(err_data_up[i])
    if i == len(instrument_data)-1:
        ax.errorbar(v_data_inst, vFv_data_inst, yerr=(err_data_inst_down, err_data_inst_up), fmt='o', label=str(list_intruments[-1]), 
                    markersize=3, elinewidth=1, color = cmap(len(list_intruments)))


ax.legend(loc='best',ncol=2)
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
    #ax1.fill_between(data[0], f_high, (data[1] - f_high(data[0]))/data[1], alpha=.5,zorder=1)
    
    ax1.hlines(0,best_model[2][0], best_model[2][-1])
    ax1.fill_between(best_model[2], f_high(best_model[2])/f_best(best_model[2]) -1, -f_best(best_model[2])/f_low(best_model[2]) +1, alpha=.5,zorder=1)
    
    list_intruments=[instrument_data[0]]
    v_data_inst = [v_data[0]]
    vFv_data_inst = [vFv_data[0]]
    err_data_inst_down = [err_data_down[0]]
    err_data_inst_up = [err_data_up[0]]
    for i in range(1,len(instrument_data)):
        if instrument_data[i] != list_intruments[-1]:
            ax1.errorbar(v_data_inst, vFv_data_inst/f_best(v_data_inst)-1, yerr=(err_data_inst_down/f_best(v_data_inst), err_data_inst_up/f_best(v_data_inst)), fmt='o', 
                        markersize=3, elinewidth=1, color = cmap(len(list_intruments)))
            list_intruments.append(instrument_data[i])
            v_data_inst = [v_data[i]]
            vFv_data_inst = [vFv_data[i]]
            err_data_inst_down = [err_data_down[i]]
            err_data_inst_up = [err_data_up[i]]
        else:   
            v_data_inst.append(v_data[i])
            vFv_data_inst.append(vFv_data[i])
            err_data_inst_down.append(err_data_down[i])
            err_data_inst_up.append(err_data_up[i])
        if i == len(instrument_data)-1:
            ax1.errorbar(v_data_inst, vFv_data_inst/f_best(v_data_inst)-1, yerr=(err_data_inst_down/f_best(v_data_inst), err_data_inst_up/f_best(v_data_inst)), fmt='o', 
                        markersize=3, elinewidth=1, color = cmap(len(list_intruments)))
    
    ax1.set_ylim(min(vFv_data/f_best(v_data)-1)*2, max(vFv_data/f_best(v_data)-1)*2)
    ax1.set_ylabel("Residuals")

new_min, new_max = blazar_plots.scale_to_values(data[1], lower_adjust_multiplier=3000, upper_adjust_multiplier=3)
ax.set_ylim(new_min, new_max)

plt.tight_layout()


plt.savefig(BASE_PATH + folder + "/user_plot_SED.pdf")

#plt.show()

#blazar_report.save_plots_and_info(configs, (v_data, vFv_data, err_data_down), param_min_vals, param_max_vals, folder=folder, 
#                                  samples=(chain, flat_samples, log_probs, flat_log_probs), use_samples=True, redshift=redshift, 
#                                  eic=eic, verbose=True)
