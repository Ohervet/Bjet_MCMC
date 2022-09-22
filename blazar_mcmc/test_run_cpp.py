"""
Parameters are always listed in the following order:
[delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R]
"""
import os
import timeit

import numpy as np

# import bj_mcmc
import blazar_model
from blazar_properties import *

"""
This file is for testing various methods of calling the C++ code
"""

constant_blob_params = (0, 1, 9.0e+17)

# redshift, hubble const, angle
transformation_params = (0.143, 71.0, 0.57)

# number of spectral points, minimal frequency, maximal frequency, file name prefix
numerical_params = (99, 50000000.0, 1e+29, "cppv2")

"""
command line parameter string:
1 and 0 at the beginning specify no prev files made and we want a model with just the blob
then: <transformation params> <our params> <constant blob params> <numerical params>
for a total of 21 params
length of emitting region, ebl absorption, distance of blob
"""
command_line_params = ("3 " + BASE_PATH + " 0 " + 3 * "{} " + "our_params" + 7 * "{} ").format(*transformation_params,
                                                                                               *constant_blob_params,
                                                                                               *numerical_params)
command_line_params = command_line_params.replace("our_params", 9 * "{} ")[:-1]

c_params1 = [3, BASE_PATH + "data", 0, *transformation_params]
c_params1 = [str(val) for val in c_params1]
c_params2 = [*constant_blob_params, *numerical_params]
c_params2 = [str(val) for val in c_params2]


# concatenate the two with params in the middle
# ^don't want the space at the end
# now, this string can be used to insert the params

# initialization for version calling c++ code
# bj_mcmc.get_transformation_parameters(*transformation_params)
# bj_mcmc.get_numerical_parameters(*numerical_params)
# bj_mcmc.set_input_mode(1)


def get_linear_params(params):
    return params[0], np.power(10, params[1]), params[2], params[3], np.power(10, params[4]), np.power(10, params[
        5]), np.power(10, params[6]), np.power(10, params[7]), np.power(10, params[8])


def original_model(params, verbose=False):
    blazar_model.create_params_file(params, "original")
    blazar_model.make_SED(verbose=verbose, use_param_file=False)
    os.remove(BASE_PATH + "parameter_files/params.txt")


def original_with_prev_model(params, verbose=False):
    blazar_model.create_params_file(params, "original")
    blazar_model.make_SED(verbose=verbose, prev_files=True, use_param_file=False)
    os.remove(BASE_PATH + "parameter_files/params.txt")


def pybind_model(params):
    params = get_linear_params(params)
    # bj_mcmc.get_blob_parameters(*params, *constant_blob_params)
    # bj_mcmc.run_models()


def command_model(params, verbose=False):
    blazar_model.make_SED(params, command_params_1=c_params1, command_params_2=c_params2, verbose=verbose)


def command_different(params, verbose=False):
    # command_params = get_linear_params(params)
    # command_params = command_line_params.format(*command_params)
    # command_params = command_params.split()
    blazar_model.make_SED(params, command_params_1=c_params1, command_params_2=c_params2, verbose=verbose,
                          executable="bjet_mcmc/bj_mcmc2")


def test_original(parameters_list):
    for p in parameters_list:
        original_model(p)


def test_original_with_prev(parameters_list):
    for p in parameters_list:
        original_with_prev_model(p)


def test_pybind(parameters_list):
    for p in parameters_list:
        pybind_model(p)


def test_command(parameters_list):
    for p in parameters_list:
        command_model(p)


def test_other(parameters_list):
    for p in parameters_list:
        command_different(p)


if __name__ == '__main__':
    cerruti_params = [96.83, 2.40, 2.0, 4.0, 2.0, 6.70, 4.725, -1.82, 16.11]  # for testing
    mcmc_params = [36.865367, 3.45330174, 2.95115932, 6.63132403, 3.06158895, 5.17967147,
                   4.83918436, -0.24782528, 16.63055046]
    params_list = [cerruti_params, mcmc_params]
    timeit_setup = "import test_run_cpp; params_list = {}".format(params_list)
    # """
    # timeit_run = "test_run_cpp.test_pybind(params_list)"
    # print("pybind timeit", timeit.timeit(timeit_run, setup=timeit_setup, number=10))
    timeit_run = "test_run_cpp.test_original(params_list)"
    print("param file timeit", timeit.timeit(timeit_run, setup=timeit_setup, number=20))
    timeit_run = "test_run_cpp.test_original_with_prev(params_list)"
    print("param file with prev", timeit.timeit(timeit_run, setup=timeit_setup, number=20))
    timeit_run = "test_run_cpp.test_command(params_list)"
    print("command timeit", timeit.timeit(timeit_run, setup=timeit_setup, number=20))
    timeit_run = "test_run_cpp.test_other(params_list)"
    print("command no print timeit", timeit.timeit(timeit_run, setup=timeit_setup, number=20))
    # """
    # original_model(mcmc_params)
    # command_model(mcmc_params)
