# Bjet_MCMC

Bjet_MCMC is a MCMC python wrapper around the C++ code Bjet. It can be used to model multiwavelength spectral energy distributions of blazars, considering one-zone synchrotron-self-Compton (SSC) model with or without the addition of external inverse-Compton process from the thermal emission of the nucleus.

The main contributors of Bjet_MCMC are: Olivier Hervet, Caitlin Johnson, and Adrian Youngquist

License
-------
The code is licensed under a [BSD-3-Clause License](LICENSE).


Acknowledgments and citation
-------
If you use this code in a scientific publication or communication, please cite the paper [Hervet at al. 2023](https://arxiv.org/abs/2307.08804) (Submitted to ApJ). For any use of the multi-zones SSC option, please also cite [Hervet at al. 2015](https://ui.adsabs.harvard.edu/abs/2015A%26A...578A..69H/abstract).




## Installation and running Bjet_MCMC:
1. Ensure all dependencies are installed, see [dependencies](#dependencies). Recommended: create the conda env "bjet-mcmc" from `environment.yml` using `conda env create -f environment.yml`. Then load your environment with `conda activate bjet-mcmc`.
2. Create a copy of `mcmc_config_template.txt` called `mcmc_config.txt`. For information, see [configuration file](#mcmc-configuration-file). 
3. Ensure the data file is formatted as described below and is at the relative path specified in `mcmc_config.txt`, see [data format](#data-format).
4. Navigate inside folder `blazar_mcmc`. Currently, code **must** be run from inside this folder or the imports won't work.
5. First time use: execute `blazar_initialize.py` with `python blazar_initialize.py` or `python3 blazar_initialize.py` depending on your Python setup. This creates all necessary folders and compiles the C++ code (see [bjet](#using-bjet_02)).
6. Execute `blazar_run_mcmc.py` with `python blazar_run_mcmc.py` or `python3 blazar_run_mcmc.py` depending on your Python setup. 

`console_data_processing.py` can be executed in an interactive Python console. It loads all information necessary for calling functions to create statistics and plots given a path relative to `BASE_PATH` for a results folder. 

By default, configurations are read from `mcmc_config.txt` and results are written to a file in `local_results` entitled `run_yyyy-mm-dd-hh:mm:ss`. 

### Dependencies
###### **Recommended:** create conda env from `environment.yml` using `conda env create -f environment.yml`


`emcee`
- To install with conda: `conda install -c conda-forge emcee`
- To install with pip: `emcee` recommends calling `pip install -U setuptools setuptools_scm pep517` and then `pip install -U emcee`
- `emcee` installation [documentation](https://emcee.readthedocs.io/en/stable/user/install/).

`scipy`
- To install with pip: `python -m pip install scipy`
- To install with conda: `conda install -c conda-forge scipy`
- To install with pip: `python -m pip install -U scipy`
- `scipy` installation [documentation](https://scipy.org/install/)

`invoke`
- To install with pip: `python -m pip install invoke`
- To install with conda: `conda install -c conda-forge invoke`


`matplotlib`
- To install with conda: `conda install matplotlib`
- To install with pip: `python -m pip install -U matplotlib`
- `matplotlib` installation [documentation](https://matplotlib.org/stable/users/installing/index.html)

`corner`
- To install with conda:`conda install -c astropy corner`
- To install with pip: `python -m pip install corner`
- `corner` installation [documentation](https://corner.readthedocs.io/en/latest/install.html)

`tqdm`
- To install with pip: `python -m pip install tqdm`

`h5py`
- `python -m pip install h5py`

**NOTE:** `numpy` is also necessary. The `emcee` installation will install `numpy` if it is not there already. 
`numpy` itself can be installed:
- Using conda: `conda install numpy`
- Using pip: `pip install numpy`

### Necessary Contents:
- Directory `sed_calculations` (do not modify or delete)
- Directory `mcmc_results`, where data will go (do not delete directory)
- Directory `parameter_files` (do not modify or delete)
- `blazar_run_mcmc.py`, `blazar_model.py`, `blazar_plots.py`, `blazar_utils.py`, `blazar_properties.py`, `blazar_report.py`

## Formatting information

### Listing parameters
Parameters are always listed in the following order:
[delta, K, n1, n2, gamma_min, gamma_max, gamma_break, B, R]
The additional parameters for EIC are bb_temp, l_nuc, tau, blob_dist, in that order.
All parameters are the logarithm of the true value except for delta, n1, and n2

---------------------------------------------------
| Variable    | Description                 | Scale  |
|-------------|-----------------------------|--------|
| delta       | Doppler factor              | linear |
| K           | particle density [cm^-3]    | log    |
| n1          | alpha_1 (first index)       | linear | 
| n2          | alpha_2 (second index)      | linear | 
| gamma_min   | low-energy cutoff           | log    | 
| gamma_max   | high-energy cutoff          | log    | 
| gamma_break | energy break                | log    | 
| B           | magnetic field strength [G] | log    | 
| R           | blob radius [cm]            | log    | 
---------------------------------------------------
| Additional for EIC | Description                  | Scale |
|--------------------|------------------------------|-------|
| bb_temp            | Black body temp of disk [K]  | log   | 
| l_nuc              | Nucleus luminosity [ergs/s]  | log   | 
| tau                | Frac of luminosity scattered | log   | 
| blob_dist          | Distance of blob [cm]        | log   | 

### Data format
Data by default is assumed to be in data files with the following columns as labels in the first row:
`!E(eV)		F(ergcm-2s-1)	delta_E(-)	delta_E(+)	delta_F(-)	delta_F(+)	instrument`.
The first row is therefore skipped when reading data.

`E`, `F`, and `delta F(-)` are the columns used by the program, which by default are in columns 0, 1, and 4.
As an optional parameter, these can be changed. The data must be space-delimited.

All rows of data **must have the same number of rows**.

### MCMC configuration file
The default configuration file is `mcmc_config.txt` in the main directory. If you would like to use a different file, in `blazar_run_mcmc.py` in main, replace `config_file=None` and/or `directory=None` with the location of the config file and/or the desired directory for results. (If the defaults are OK, no code should need editing.)

There is a file named `mcmc_config_template.txt` as an example. Make a copy of this file named `mcmc_config.txt`. 
This file is automatically in `.gitignore` since it changes locally.

Rows are:
```
description = <description of run>
folder_label = <prefix for folder name>
eic = <True/False, SSC or SSC + EIC>
data_file=<file name relative to blazars-mcmc>
n_steps=<# of steps>
n_walkers=<# of walkers>
discard=<steps to discard>
parallel=<True or False (use parallel processing)
cores = <# of cores to use if parallel True, optional>
tau_var=<Tau variability in hours>
redshift = <float for redshift>
custom_alpha2_limits=False or <True, lower, upper>

```
Notes: 
- Labels must be exactly as listed
- `data_file` is the *relative path* from the home directory to the file with data.
#### Example mcmc configuration file
(file named `mcmc_config.txt`, can be modified)
```
# Configuration file for running mcmc

# description:
description=3C66A b2 eic

# folder label:
folder_label =3C66A_b2_eic

# eic (True/False):
eic=True

# Data file:
data_file=real_data/3C66A_sed_block2.dat

# Number of steps:
n_steps=4000

# Number of walkers:
n_walkers=300

# Discard number:
discard=500

# Parallel processing (True/False):
parallel=True

# If parallel processing, # of cores (will use # of cores - 1 if not specified)
#cores=9

# use tau variability (boolean):
use_variability= True

# tau variability (in hours):
tau_variability = 168

# redshift (J1010 0.143, 3C66A 0.340)
redshift = 0.340

# Custom alpha2 limits (True/False, <val1>, <val2>) val1 and val2 optional
custom_alpha2_limits=False
```
## Code information
The `blazar_run_mcmc.py` script creates a results folder, by default in local_results named with the folder_prefix and the current date.

At the start, a `basic_info.txt` doc is created in the folder with configuration information saved.

The history of the entire chain from the MCMC and the corresponding log_probs are stored in an h5 file named `backend.h5`. This is updated continuously. This means that if the code crashes, all existing results are still present. h5 files can be re-loaded as the current state for the MCMC at the beginning of another run (see the main method in `blazar_run_mcmc.py`) or loaded to use for statistics or plotting (see `load_backend()` in `blazar_report.py`).
### Using the blazar_mcmc code
Using the functions directly rather than running the `blazar_run_mcmc.py` file as a script allows for greater customizability and increased functionality.

See docstrings of the modules for function information and the docstrings of individual functions for more detailed information. 


General note: paths are generally specified relative to the base directory (`blazars-mcmc`) or in more specific function to the parent folder (`RESULTS_FOLDER` or `DATA_FOLDER`). See function documentation.
#### Modules:
- `blazar_clean.py`: Deletes old model data files
- `blazar_initialize.py`: Compiles the c++ code, ensures all necessary directories are present. 
- `blazar_model.py`: Functions for creating models
- `blazar_plots.py`: Functions for plotting
- `blazar_properties.py`: Advanced configuration options and default values
- `blazar_report.py`: Functions for calculating stats and creating reports
- `blazar_run_mcmc.py`: Main python file for running the MCMC.
- `blazar_utils.py`: Utility functions

#### Notable Options
In the main method in `blazar_run_mcmc.py`, it is possible to specify the starting position (p0) for the MCMC. This is useful for continuing a run after a crash.

The `save()` function in `blazar_report.py` will create all plots and statistics given just the name of the directory that contains an h5 file and an `info.txt` doc. 

*Make models, plot results, and get statistics either in the Python console or in scripts:* 
- `blazar_report.parse_info_doc()` will, given a path to an `info.txt` file or `basic_info.txt` file, parse the file and return a dictionary with all information. Notably, `info["configs"]` contains the dictionary of configurations.
- `blazar_utils.read_configs()` will, given the path to an `mcmc_config.txt` (default is `"mcmc_config.txt"`), return a dictionary of configurations.
- `blazar_utils.read_data()` returns a tuple of numpy arrays with the `v_data`, `vFv_data`, and the `err_data` given a data file.
- `blazar_report.load_from_backend()` will return a tuple of the entire chain from the MCMC and the log_probs given a data folder relative to `RESULTS_FOLDER`
- After using the above in a script or from the console, statistics and plots can be calculated from functions in `blazar_plots.py` and `blazar_report.py`. 
## Using bjet_core

`bjet_core` contains a modified version of the Bjet_02 code written by Olivier Hervet that is used for the MCMC. 

Compile the exec with `make bj_core`.

Code is in the `bjet_core` folder. The necessary files are `bj_core.cpp`, `bj_core.h`, `processes_supp_core.cpp`, and `processes_supp_core.h`
### Usage
#### Options for calling the executable:
- Usage 1: `bj_core --help`
- Usage 2: `bj_core <parameter file>`
- Usage 3: `bj_core 0 <parameter file>` same execution as usage 2, `_prev files` made
- Usage 4: `bj_core 1 <parameter file>`                           no _prev files made
- Usage 5: `bj_core 2 <data folder> <parameter file>`            _prev files made
- Usage 6: `bj_core 3 <data folder> <parameter file>`               no _prev files made
- Usage 7: `bj_core 0/1 <model type> <model params, at least 19 depending on model>`
    ^params from command line, 0 = yes `_prev files`, 1 = no
- Usage 8: `bj_core 2/3 <data folder> <model type> <model params, at least 19 depending on model>`
    ^params from command line, 2 = yes `_prev files`, 3 = no

#### Parameter files
Parameter files specified in the command line should be paths relative to the `bj_core` executable. See `params.txt` in `bjet_core` for an example parameter file.
#### Modes:
- 0: `_prev` files are made, save in default data folder
- 1: no `_prev` files are made, save in default data folder
- 2: `_prev` files are made, save in specified data folder
- 3: no _prev files are made, save in specified data folder

(default data directory is in the same directory as the executable, named "data")

#### Valid values for model type:
- 0: model with just blob
- 1: blob + EIC parameters
- Other models not yet implemented

#### Order of parameters in command line:
##### Model type 0:
 ```
 bj_core <0, 1, 2, 3> <data folder if applicable> 0 (model type)
 z, H_0, THETA, DOP_B, K1, N1, N2, GAMMA_MIN, GAMMA_MAX, GAMMA_BRK, B, R_src,
 L_src, IIR_level, D_b, NU_DIM, NU_STR, NU_END, prefix
 ```
 argc should be 22 or 23 depending on if data folder is listed

##### Model type 1:
 ```
 bj_core <0, 1, 2, 3> <data folder if applicable> 1 (model type)
 z, H_0, THETA, DOP_B, K1, N1, N2, GAMMA_MIN, GAMMA_MAX, GAMMA_BRK, B, R_src,
 L_src, IIR_level, D_b, T_BB, TBB_tor, L_nuc, tau, L_tor, tau, NU_DIM, NU_STR, NU_END, prefix
 *Note that tau is present twice, this is a slight error in the bjet code. The second tau value is not used for
 anything, but it must be inputted.
```
 argc should be 28 or 29 depending on if data folder is listed

Example:
```
 ./bj_core 3 /Users/sed_calculations 1 0.34 69.6 0.57 50.0 612.1 2.28 3.74 2816.9 1803000 44806 0.00236 5.94e+17 0 1 3.8e+15 2013 2.0e+4 1.7e+21 1.5e-10 5.5e+20 9.0e-5 99 50000000.0 1e+29 run
```
 ^ here, the 3 indicates that the data folder is specified and no prev file is made. 1 is the EIC model type. Then 0.34 is z (redshift) and then the rest of the parameters are enumerated.
