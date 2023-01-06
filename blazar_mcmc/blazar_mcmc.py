import datetime
import multiprocessing
import shutil

import emcee

import blazar_utils
import blazar_report
import blazar_properties

import glob
import numpy as np
import pathlib
import subprocess
import os
import random
from scipy import interpolate

PROGRAM_NAME = "Bjet_MCMC"
RESULTS_FOLDER = "local_results"
DATA_FOLDER = "sed_calculations"
EXECUTABLE = "bjet_mcmc/bj_mcmc"
BASE_PATH = blazar_properties.BASE_PATH

configs = blazar_utils.read_configs()
v_data, vFv_data, err_data = blazar_utils.read_data(configs["data_file"])

EIC = configs["eic"]
DIM = 13 if EIC else 9
if EIC:
    PARAM_IS_LOG = [False, True, False, False, True, True, True, True, True, True, True, True, True]
else:
    PARAM_IS_LOG = [False, True, False, False, True, True, True, True, True]

redshift = configs["redshift"]
use_variability = configs["use_variability"]
tau = configs["tau_variability"]
alpha2_limits = configs["alpha2_limits"]

model_type = "1" if EIC else "0"
settings_and_transformation = ["3", BASE_PATH + DATA_FOLDER, model_type, str(redshift), "69.6", "0.57"]
constant_and_numerical = ["0", "1", "9.0e+17", "99", "5.0e+7", "1e+29", "temp_stem"]

# mins maxes
param_min_vals = [1., 0., 1., float(alpha2_limits[0]), 0., 3., 2., -4., 14.]
param_max_vals = [100., 8., 5., float(alpha2_limits[1]), 5., 8., 6.699, 0., 19.]
if EIC:
    extra_min = [1.0, 20.0, -10.0, 10.0]
    extra_max = [6.0, 50.0, 0.0, 21.0]
    param_min_vals = param_min_vals + extra_min
    param_max_vals = param_max_vals + extra_max
param_min_vals = np.array(param_min_vals)
param_max_vals = np.array(param_max_vals)


def make_model(params, name_stem="run"):
    # convert to log
    params = params * 1.0
    for i in range(len(params)):
        if PARAM_IS_LOG[i]:
            params[i] = np.power(10, params[i])
    if not EIC:
        command = settings_and_transformation + [str(val) for val in params] + constant_and_numerical[:-1] + [name_stem]
    else:
        command = settings_and_transformation + [str(val) for val in params[:9]] + constant_and_numerical[:2] + [str(params[-1])]
        command = command + [str(params[9]), "2.0e+4", str(params[10]), str(params[11]), "5.5e+20", "9.0e-5"]  # EIC components
        command = command + constant_and_numerical[3: -1] + [name_stem]
    subprocess.run([BASE_PATH + EXECUTABLE, *command], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))

    loaded_model = np.loadtxt(BASE_PATH + DATA_FOLDER + "/" + name_stem + "_ss.dat", delimiter=' ')
    logv = loaded_model[:, 0]
    logvFv = loaded_model[:, 2]
    vFv = np.power(10, logvFv)

    stems = ["cs"] if not EIC else ["cs", "ecs", "cs2", "nuc"]
    for s in stems:
        loaded_model = np.loadtxt(BASE_PATH + DATA_FOLDER + "/" + name_stem + "_" + s + ".dat", delimiter=' ')
        model_logv = loaded_model[:, 0]
        model_logvFv = loaded_model[:, 2]
        current_logv, current_logvFv, current_vFv = logv, logvFv, vFv

        new_lower = np.where(model_logv < logv[0])[0]
        new_higher = np.where(model_logv > logv[-1])[0]
        logv = np.concatenate((model_logv[new_lower], current_logv, model_logv[new_higher]))
        logvFv = np.concatenate((model_logvFv[new_lower], current_logvFv, model_logvFv[new_higher]))
        vFv = np.power(10, logvFv)

        overlap_start = np.where(logv >= max(model_logv[0], current_logv[0]))[0][0]
        overlap_end = np.where(logv <= min(model_logv[-1], current_logv[-1]))[0][-1]

        interpolation = interpolate.interp1d(model_logv, np.power(10, model_logvFv))

        new_vFv = np.concatenate((np.zeros(overlap_start), interpolation(logv[overlap_start:overlap_end + 1]),
                                  np.zeros(len(logv) - overlap_end - 1)))

        vFv = vFv + new_vFv
        logvFv = np.log10(vFv)
    return logv, logvFv, np.power(10, logv), vFv


def log_prior(params):
    delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R, *other_params = params
    if n1 > n2 or gamma_min > gamma_max or gamma_break < gamma_min or gamma_break > gamma_max:
        return -np.inf
    # testing if between min and max
    for i in range(len(params)):
        # gamma break is solely constrained by gamma_min and gamma_max
        if i != 6 and (param_min_vals[i] > params[i] or param_max_vals[i] < params[i]):
            return -np.inf
    if use_variability:
        tau_var = tau * 60 * 60  # to seconds
        c = 2.997924 * 1.0e+10
        R = np.power(10, R)
        if tau_var < (1 + redshift) / c * R / delta:
            return -np.inf
    return 0


def random_params():
    parameter_size = param_max_vals - param_min_vals
    parameters = param_min_vals + parameter_size * np.random.rand(DIM)
    while not np.isfinite(log_prior(parameters)):
        parameters = param_min_vals + parameter_size * np.random.rand(DIM)
    return parameters


def log_prob(params):
    name_stem = "run_" + str(random.getrandbits(60))

    # prior
    delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R, *other_params = params
    if n1 > n2 or gamma_min > gamma_max or gamma_break < gamma_min or gamma_break > gamma_max:
        return -np.inf
    # testing if between min and max
    for i in range(len(params)):
        # gamma break is solely constrained by gamma_min and gamma_max
        if i != 6 and (param_min_vals[i] > params[i] or param_max_vals[i] < params[i]):
            return -np.inf
    if use_variability:
        tau_var = tau * 60 * 60  # to seconds
        c = 2.997924 * 1.0e+10
        R = np.power(10, R)
        if tau_var < (1 + redshift) / c * R / delta:
            return -np.inf

    model = make_model(params, name_stem)

    # calculate chi squared
    func = interpolate.interp1d(model[0], model[1], fill_value='extrapolate')
    chi_sq = np.sum(((vFv_data - np.power(10, func(np.log10(v_data)))) / err_data) ** 2.)
    for f in glob.glob(BASE_PATH + DATA_FOLDER + "/" + name_stem + "_*"):
        os.remove(f)

    return -0.5 * chi_sq


def mcmc(p0=None):
    if p0 is not None:
        print(p0[0])
    if "folder_label" in configs:
        folder_label = configs["folder_label"]
    else:
        folder_label = "run"
    if "description" in configs:
        description = configs["description"]
    else:
        description = None

    now = datetime.datetime.now()
    date_string = now.strftime("%Y-%m-%d-%H:%M:%S")
    directory = RESULTS_FOLDER + "/" + folder_label + "_" + date_string

    os.mkdir(BASE_PATH + directory)

    backend = directory + "/backend.h5"

    # make file with basic info
    with open(BASE_PATH + directory + "/basic_info.txt", 'w') as f:
        f.write("folder name: ")
        f.write(directory)
        if description is not None:
            f.write("\nreport description: ")
            f.write(description)
            f.write("\n")
        f.write("\nconfigurations:\n")
        f.write(str(configs))
        if p0 is not None:
            f.write("\np0 given\n")

    if configs["parallel"]:
        if "cores" not in configs:
            pool = multiprocessing.Pool()
        else:
            pool = multiprocessing.Pool(processes=configs["cores"])
    else:
        pool = None
    backend = emcee.backends.HDFBackend(BASE_PATH + backend)
    backend.reset(configs["n_walkers"], DIM)

    if p0 is None:
        p0 = np.array([random_params() for _ in range(configs["n_walkers"])])

    sampler = emcee.EnsembleSampler(configs["n_walkers"], DIM, log_prob, backend=backend, moves=[(emcee.moves.StretchMove(live_dangerously=True), 1.)], pool=pool)

    print("starting mcmc")
    start = datetime.datetime.now()
    sampler.run_mcmc(p0, configs["n_steps"], progress=True)
    end = datetime.datetime.now()
    if configs["parallel"]:
        pool.close()

    with open(BASE_PATH + directory + "/basic_info.txt", 'a') as f:
        f.write("\ntime: ")
        f.write(str(end - start))
        f.write("\n")

    blazar_report.show_results(sampler, str(end - start), configs=configs)
    blazar_report.save_plots_and_info(configs, (v_data, vFv_data, err_data), param_min_vals, param_max_vals,
                                      folder=directory, sampler=sampler, use_sampler=True, description=description,
                                      time=str(end - start), redshift=redshift, eic=EIC)

    return sampler, directory


if __name__ == '__main__':
    """
    p0_file = "local_results/3C66A_b6_eic_2022-06-08-20:17:26/backend.h5"
    reader = emcee.backends.HDFBackend(BASE_PATH + p0_file, read_only=True)
    p0 = reader.get_last_sample().coords
    """
    p0 = None
    sampler, directory = mcmc(p0=p0)
    if blazar_properties.TMP:
        shutil.move(BASE_PATH + directory, blazar_properties.FOLDER_PATH + directory)
        shutil.rmtree(blazar_properties.TEMP_DIR)
