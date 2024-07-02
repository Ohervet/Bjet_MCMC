#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: blazar_report.py

Contains functions for showing results of the MCMC, given the results. It
calculates indices and parameters within 1sigma, finds the best parameters,
creates a written text report, and creates all necessary plots.

Functions:
Function save_plots_and_info (given the configurations, the data, and
the minima/maxima for the parameters, and a source of MCMC data) creates a folder
with an info.txt file and plots:
    - Corner plot
    - Plots of chi squared (med, best, all by step)
    - Plots of the best model with data and models within 1 sigma (random
      models, extreme models, and a combination of both)
    - Plot of the chain

Parameters are always listed in the following order:
[delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R]
All parameters are the logarithm of the true value except for delta, n1, and n2
---------------------------------------------------
delta       doppler factor                  linear
K           particle density [cm^-3]        log
n1          alpha_1 (first index)           linear
n2          alpha_2 (second index)          linear
gamma_min   low-energy cutoff               log
gamma_max   high-energy cutoff              log
gamma_break energy break                    log
B           magnetic field strength [G]     log
R           blob radius (cm)                log
---------------------------------------------------
"""
import glob
import os
import random

import emcee
import numpy as np
from scipy import stats

import blazar_model
import blazar_plots
import blazar_utils
from blazar_properties import *


# utils ----------------------------------------------------------------------------------------------------------------


def get_indices_within_1sigma(values, use_log_probs=True, alpha=0.32, eic=False):
    """
    Arguments:
    values (2D np array, _ rows, NUM_DIM columns): flat samples
    dimensions:
    use_log_probs:
    alpha:
    """
    if use_log_probs:
        chi_squared = -2 * values
    else:
        chi_squared = 1 * values
    min_chi_squared_index = np.argmin(chi_squared)
    min_chi_squared = chi_squared[min_chi_squared_index]
    delta_chi_squared = stats.chi2.ppf(1 - alpha, modelProperties(eic).NUM_DIM)
    indices_within_1sigma = np.where(chi_squared < min_chi_squared + delta_chi_squared)[0]
    return indices_within_1sigma


def min_max_params_1sigma(flat_data, flat_log_probs, eic=False):
    n_dim = len(flat_data[0])
    indices_within_1sigma = get_indices_within_1sigma(flat_log_probs, eic=eic)
    mins = np.copy(flat_data[indices_within_1sigma[0]])
    maxes = np.copy(flat_data[indices_within_1sigma[0]])
    for i in range(len(indices_within_1sigma)):
        for j in range(n_dim):
            val = flat_data[indices_within_1sigma[i]][j]
            if val < mins[j]:
                mins[j] = val
            if val > maxes[j]:
                maxes[j] = val
    return mins, maxes


def get_best_log_prob_and_params(configs=None, sampler=None, log_probs=None, flat_chain=None, discard=None, fixed_params=None):
    if (sampler is None) and (log_probs is None or flat_chain is None):
        raise ValueError("sampler or log_probs+chain not provided")
    if sampler is not None:
        if configs is None and discard is None:
            raise ValueError("sampler provided but no configs or discard")
        if discard is None:
            discard = configs["discard"]
        log_probs = sampler.get_log_prob(flat=True, discard=discard)
        flat_chain = sampler.get_chain(flat=True, discard=discard)
    maximum_index = np.argmax(log_probs)
    maximum_log_prob = log_probs[maximum_index]
    best_params = flat_chain[maximum_index]
    if fixed_params:
        for i in range(len(fixed_params)):
            if fixed_params[i] != -np.inf:
                best_params[i] = fixed_params[i]
    return maximum_log_prob, best_params


# report displays ------------------------------------------------------------------------------------------------------
def make_text_file(output_file, configs, data_points, backend_file=None, reader=None, sampler=None, use_sampler=False,
                   samples=None, use_samples=False, description=None, time=None, p0_source=None, acceptance_frac=None,
                   eic=False, fixed_params=None):
    if use_samples and (samples is not None):
        samples, flat_chain, _, log_probs = samples
    # note that samplers and readers both have the functions needed
    else:
        if use_sampler and (sampler is not None):
            data_source = sampler
        elif backend_file is not None:
            data_source = emcee.backends.HDFBackend(BASE_PATH + backend_file, read_only=True)
        elif reader is not None:
            data_source = reader
        else:
            raise ValueError("backend file or reader not provided")

        samples = data_source.get_chain(discard=configs["discard"])
        flat_chain = data_source.get_chain(flat=True, discard=configs["discard"])
        log_probs = data_source.get_log_prob(flat=True, discard=configs["discard"])

    with open(BASE_PATH + output_file, 'w') as f:
        tau = emcee.autocorr.integrated_time(samples, quiet=True)
        maximum_log_prob, best_params = get_best_log_prob_and_params(log_probs=log_probs, flat_chain=flat_chain)
        print("best parameters:", best_params)
        min_1sigma_params, max_1sigma_params = min_max_params_1sigma(flat_chain, log_probs, eic=eic)
        if description is not None:
            f.write("report description: ")
            f.write(str(description))
            f.write("\n")
        if time is not None:
            f.write("time: ")
            f.write(time)
            f.write("\n")
        if p0_source is not None:
            f.write("p0 from ")
            f.write(p0_source)
            f.write("\n\n")

        f.write("configurations: \n")
        f.write(str(configs))
        f.write("\n\n")

        f.write(
            text_report(best_params, min_1sigma_params, max_1sigma_params, -2 * maximum_log_prob, data_points, eic=eic,
                        fixed_params=fixed_params))
        f.write("\n\n")

        f.write("best params: ")
        f.write(str(best_params))
        f.write(", chi squared = ")
        f.write(str(-2 * maximum_log_prob))
        f.write("\n\n")

        f.write("min_1sigma_params: ")
        f.write(str(min_1sigma_params))
        f.write("\nmax_1sigma_params: ")
        f.write(str(max_1sigma_params))
        f.write("\n\n")

        if acceptance_frac is not None:
            f.write("acceptance fraction: avg = ")
            f.write(str(np.average(acceptance_frac)))
            f.write("\n\n")

        f.write("autocorrelation time: avg = ")
        f.write(str(np.mean(tau))+" steps")
        f.write("\n\n")


def text_report(best, min_1sigma, max_1sigma, best_chi_sq, data_points, param_names=None, eic=False, fixed_params=None):
    if param_names is None:
        param_names = modelProperties(eic, fixed_params=fixed_params).PARAM_NAMES
    log = modelProperties(eic).PARAM_IS_LOG
    #remove any frozen parameter from the log list
    if fixed_params:
        fixed_params2 = fixed_params.copy()
        i = 0
        while i < len(fixed_params2):
          if fixed_params2[i] != -np.inf:       
            del log[i]
            del fixed_params2[i]
          else:
            i+=1
    format_string = "{:^13}{:^15}{:^20}\n"
    num_format = "{:.2e}"
    range_format = num_format + " - " + num_format
    output = format_string.format("Parameter", "Best Value", "1sigma Range")
    for i in range(len(best)):
        b = best[i]
        s1 = min_1sigma[i]
        s2 = max_1sigma[i]
        if log[i]:
            b = np.power(10, b)
            s1 = np.power(10, s1)
            s2 = np.power(10, s2)
        output += format_string.format(param_names[i], num_format.format(b), range_format.format(s1, s2))

    dims = modelProperties(eic).NUM_DIM
    output += "Reduced chi squared: {:.2f} / {} = {:.2f}\n".format(best_chi_sq, data_points - dims,
                                                                   best_chi_sq / (data_points - dims))
    return output


def show_results(sampler, time, configs=None, discard=None):
    if configs is None and discard is None:
        raise ValueError("Neither configs nor discard provided")
    if discard is None:
        discard = configs["discard"]
    chain = sampler.get_chain(discard=discard)
    tau = emcee.autocorr.integrated_time(chain, quiet=True)
    maximum_log_prob, best_params = get_best_log_prob_and_params(sampler=sampler, discard=discard)

    print("MCMC Results:")
    print("Best log_prob:", maximum_log_prob)
    print("Best parameters:", best_params)
    print("Time:", time)
    print("Mean auto-correlation time:", np.mean(tau))


def save_plots_and_info(configs, data, param_min_vals, param_max_vals, folder=None, parent_folder=None,
                        backend_file=None, reader=None, sampler=None, use_sampler=False, samples=None,
                        use_samples=False, description=None, time=None, p0_source=None, acceptance_frac=None,
                        theta=None, redshift=None, min_freq=None, max_freq=None, executable=None, data_folder=None,
                        command_params_full=None, command_params_1=None, command_params_2=None, torus_temp=None,
                        torus_luminosity=None, torus_frac=None, verbose=False, eic=False):
    if backend_file is None and reader is None and (not use_sampler or sampler is None) and not (use_samples or samples is None):  # source not given
        raise ValueError("backend file or reader not provided")
    if reader is None and not use_sampler and not use_samples:
        reader = emcee.backends.HDFBackend(BASE_PATH + backend_file, read_only=True)
    if folder is None and parent_folder is None:
        raise ValueError("folder or parent folder not provided")
    if folder is None:
        folder = parent_folder + "/new_report"
        os.mkdir(folder)

    if use_sampler:
        chain = sampler.get_chain(discard=configs["discard"])
        flat_chain = sampler.get_chain(flat=True, discard=configs["discard"])
        log_probs = sampler.get_log_prob(discard=configs["discard"])
        flat_log_probs = sampler.get_log_prob(flat=True, discard=configs["discard"])
        acceptance_frac = sampler.acceptance_fraction
    elif use_samples:
        chain, flat_chain, log_probs, flat_log_probs = samples
    else:
        chain = reader.get_chain(discard=configs["discard"])
        flat_chain = reader.get_chain(flat=True, discard=configs["discard"])
        log_probs = reader.get_log_prob(discard=configs["discard"])
        flat_log_probs = reader.get_log_prob(flat=True, discard=configs["discard"])

    best_log_prob, best_params = get_best_log_prob_and_params(
        log_probs=flat_log_probs, flat_chain=flat_chain)
    min_1sigma_params, max_1sigma_params = min_max_params_1sigma(flat_chain, flat_log_probs, eic=eic)
    indices_within_1sigma = get_indices_within_1sigma(flat_log_probs, eic=eic)
    name_stem = NAME_STEM + str(random.getrandbits(20))
    
    command_params_1, command_params_2 = blazar_model.command_line_sub_strings(name_stem=name_stem, redshift=redshift, 
                                                                               prev_files=False, eic=eic)
    command_params_2[3] = "99"  # number of points used to make SED
    
    
    
    model = blazar_model.make_model(best_params, name_stem=name_stem, theta=theta, redshift=redshift, min_freq=min_freq,
                                    max_freq=max_freq, torus_temp=torus_temp, torus_luminosity=torus_luminosity,
                                    torus_frac=torus_frac, data_folder=data_folder, executable=executable,
                                    command_params_full=command_params_full, command_params_1=command_params_1,
                                    command_params_2=command_params_2, prev_files=False, use_param_file=False,
                                    verbose=verbose, eic=eic, fixed_params=configs["fixed_params"], folder=folder)
    

    for f in glob.glob(BASE_PATH + DATA_FOLDER + "/" + name_stem + "_*"):
        os.remove(f)
    make_text_file(folder + "/info.txt", configs, len(data[0]), backend_file=backend_file, reader=reader,
                   sampler=sampler, use_sampler=use_sampler, samples=samples, use_samples=use_samples,
                   description=description, time=time, p0_source=p0_source, acceptance_frac=acceptance_frac, eic=eic,
                   fixed_params=configs["fixed_params"])
    blazar_plots.plot_model_and_data(model, configs["data_file"], flat_chain, indices_within_1sigma,
                                     redshift=redshift,eic=eic,lower_adjust_multiplier=20, file_name=folder + "/model_and_data.svg", 
                                     title=None, no_title=True, save=True, show=False, 
                                     fixed_params=configs["fixed_params"])
    
    #need to run tests with fixed params before releasing it
    blazar_plots.plot_particle_spectrum(best_params, min_1sigma_params, max_1sigma_params, configs["fixed_params"],
                                        file_name= folder+"/particle_spectrum.svg", save=True, show=False)
    
    blazar_plots.corner_plot(flat_chain, param_min_vals, param_max_vals, best_params, min_1sigma_params,
                              max_1sigma_params, file_name=folder + "/corner_plot.svg", save=True, show=False, eic=eic,
                              fixed_params=configs["fixed_params"])
    #This output is not easy to read, so not fully relevant. removed for now
    # blazar_plots.plot_chain(chain, file_name=(folder + "/plot_of_chain.jpeg"), save=True, show=False,
    #                         eic=eic)  # too much stuff for svg
    blazar_plots.plot_chi_squared(log_probs, configs["discard"], plot_type='med',
                                  file_name=(folder + "/chi_squared_plot_med.svg"), save=True, show=False)
    blazar_plots.plot_chi_squared(log_probs, configs["discard"], plot_type='best',
                                  file_name=(folder + "/chi_squared_plot_best.svg"), save=True, show=False)
    blazar_plots.plot_chi_squared(log_probs, configs["discard"], plot_type='all',
                                  file_name=(folder + "/chi_squared_plot_all.jpeg"), save=True, show=False)  # svg is big
    
    blazar_plots.plot_cooling_times(folder + "/bjet.log", best_params, fixed_params=configs["fixed_params"], 
                                    file_name= folder + "/cooling_time_obs(Thomson).svg", save=True, show=False, eic=eic, redshift=redshift)

def parse_info_doc(info_doc, info=None):
    if info is None:
        info = {}
    
    labels = ["folder name:", "report description:", "time:", "p0 ", "configurations:", "Reduced chi squared:", "best params:", "min_1sigma_params:",
              "max_1sigma_params:", "acceptance fraction:", "tau: avg", "config file:", "prev_files:"]

    def match(to_match):
        for l in labels:
            if to_match.strip().find(l) == 0:
                return l
        else:
            return None

    with open(BASE_PATH + info_doc, 'r') as f:
        for line in f:
            matched = match(line)
            if matched is None:
                continue
            elif matched == "folder name:":
                info["folder_name"] = line.strip()[len("folder name:") + 1:]
            elif matched == "report description:":
                info["description"] = line.strip()[len("report description:") + 1:]
            elif matched == "time:":
                info["time"] = line.strip()[6:]
            elif matched == "p0 ":
                info["p0_label"] = line.strip()[8:]
            elif matched == "configurations:":
                line2 = f.readline()
                info["configs"] = blazar_utils.read_configs(config_string=line2)
            elif matched == "Reduced chi squared:":
                info["reduced_chi_sq"] = float(line.strip().split()[-1].strip())
            elif matched == "best params:":
                best_params_string_list = line[line.strip().find('[') + 1: line.rfind(']')].split()
                best_params = []
                for elem in best_params_string_list:
                    best_params.append(float(elem.strip()))
                info["best_params"] = np.array(best_params)
                info["best_chi_sq"] = float(line.strip().split()[-1].strip())
            elif matched == "min_1sigma_params:":
                min_1sigma_string_list = line[line.strip().find('[') + 1: line.rfind(']')].split()
                min_1sigma_params = []
                for elem in min_1sigma_string_list:
                    min_1sigma_params.append(float(elem.strip()))
                info["min_1sigma_params"] = np.array(min_1sigma_params)
            elif matched == "max_1sigma_params:":
                max_1sigma_string_list = line[line.strip().find('[') + 1: line.rfind(']')].split()
                max_1sigma_params = []
                for elem in max_1sigma_string_list:
                    max_1sigma_params.append(float(elem.strip()))
                info["max_1sigma_params"] = np.array(max_1sigma_params)
            elif matched == "acceptance fraction:":
                info["acceptance_frac"] = float(line.strip().split()[-1].strip())
            elif matched == "tau: avg":
                info["tau_avg"] = float(line.strip().split()[-1].strip())
            elif matched == "config file:":
                info["config_file"] = line.strip()[len("config file:") + 1:]
            elif matched == "prev_files:":
                temp = line.strip().split()
                info["prev_files"] = temp[1].strip()[0] == "T"
                if len(temp) >= 4:
                    info["use_param_files"] = temp[3].strip()[0] == "T"
    return info


# will create all results from just a folder w/ a backend file in it
def save(folder, description=None, time=None, p0_source=None, acceptance_frac=None, data_file=None, configs=None,
         text_file_name="info.txt", text_only=False, theta=None, redshift=None, min_freq=None, max_freq=None,
         torus_temp=None, torus_luminosity=None, torus_frac=None, data_folder=None, executable=None,
         command_params_full=None, command_params_1=None, command_params_2=None, verbose=False, eic=False, fixed_params=None):
    # WILL USE CURRENT CONFIGURATIONS IF NONE FOUND
    info = {}
    if os.path.exists(BASE_PATH + folder + "/info.txt"):
        info = parse_info_doc(folder + "/info.txt")
    elif os.path.exists(BASE_PATH + folder + "/basic_info.txt"):
        info = parse_info_doc(folder + "/basic_info.txt")
    # get data
    if "configs" in info:
        configs = info["configs"]
        data = blazar_utils.read_data(configs["data_file"])
        if "eic" in configs:  # not true for some old files
            eic = configs["eic"]
    else:
        if configs is None:
            configs = blazar_utils.read_configs()
        if data_file is not None:
            data = blazar_utils.read_data(data_file)
        else:
            data = blazar_utils.read_data(configs["data_file"])

    if "time" in info:
        time = info["time"]
    if "p0_label" in info:
        p0_source = info["p0_label"]
    if "acceptance_frac" in info:
        acceptance_frac = info["acceptance_frac"]
    if "description" in info:
        description = info["description"]
    param_min_vals, param_max_vals = blazar_utils.min_max_parameters(alpha2_limits=configs["alpha2_limits"], eic=eic)

    if text_only:
        make_text_file(folder + "/" + text_file_name, configs, len(data[0]), backend_file=folder + "/backend.h5",
                       description=description, time=time, p0_source=p0_source, acceptance_frac=acceptance_frac,
                       eic=eic, fixed_params=fixed_params)
    else:
        save_plots_and_info(configs, data, param_min_vals, param_max_vals, folder=folder,
                            backend_file=folder + "/backend.h5", description=description, time=time,
                            p0_source=p0_source, acceptance_frac=acceptance_frac, theta=theta, redshift=redshift,
                            min_freq=min_freq, max_freq=max_freq, executable=executable, data_folder=data_folder,
                            command_params_full=command_params_full, command_params_1=command_params_1,
                            command_params_2=command_params_2, torus_temp=torus_temp, torus_luminosity=torus_luminosity,
                            torus_frac=torus_frac, verbose=False, eic=eic)


def load_from_backend(folder, flat=False):
    reader = emcee.backends.HDFBackend(BASE_PATH + RESULTS_FOLDER + '/' + folder + '/backend.h5')
    return reader.get_chain(flat=flat), reader.get_log_prob(flat=flat)
