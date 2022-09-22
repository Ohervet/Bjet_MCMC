# Bjet mcmc

This code is a slightly modified version of `bj02` by Olivier Hervet. 
Modified by Sarah Youngquist in February 2022, updated June 2022.

## Usage

### Calling the executable:
- Usage 1: `bj_mcmc --help`
- Usage 2: `bj_mcmc <parameter file>`
- Usage 3: `bj_mcmc 0 <parameter file>`:                             same execution as usage 2, _prev files made
- Usage 4: `bj_mcmc 1 <parameter file>`:                             no _prev files made
- Usage 5: `bj_mcmc 2 <data folder> <parameter file>`:               _prev files made
- Usage 6: `bj_mcmc 3 <data folder> <parameter file> `:              no _prev files made
- Usage 7: `bj_mcmc 0/1 <model type> <model params, at least 19 depending on model>`:
    params from command line, 0 = yes _prev files, 1 = no
- Usage 8: `bj_mcmc 2/3 <data folder> <model type> <model params, at least 19 depending on model>`:
    params from command line, 2 = yes _prev files, 3 = no

#### Modes:
- 0: _prev files are made, save in default data folder
- 1: no _prev files are made, save in default data folder
- 2: _prev files are made, save in specified data folder
- 3: no _prev files are, save in specified data folder

#### Valid values for model type:
- 0: model with just blob
- 1: blob + EIC parameters
- Other models not yet implemented

(default data directory is in the same directory as the executable, named "data")


### Order of parameters in command line:

#### Model type 0:
 After `bj_mcmc <0, 1, 2, 3> <data folder if applicable> 0`:

 `z, H_0, THETA, DOP_B, K1, N1, N2, GAMMA_MIN, GAMMA_MAX, GAMMA_BRK, B, R_src,
 L_src, IIR_level, D_b, NU_DIM, NU_STR, NU_END, prefix`

 argc should be 22 or 23 depending on if data folder is listed

#### Model type 1:
 After `bj_mcmc <0, 1, 2, 3> <data folder if applicable> 1`:

 `z, H_0, THETA, DOP_B, K1, N1, N2, GAMMA_MIN, GAMMA_MAX, GAMMA_BRK, B, R_src,
 L_src, IIR_level, D_b, T_BB, TBB_tor, L_nuc, tau, L_tor, tau, NU_DIM, NU_STR, NU_END, prefix`

 *Note that tau is present twice, this is a slight error in the bjet code. The second tau value is not used for
 anything, but it must be inputted.

 argc should be 28 or 29 depending on if data folder is listed

#### Example:

 ```
 ./bj_mcmc 3 /Users/sed_calculations 1 0.34 69.6 0.57 50.0 612.1 2.28 3.74 2816.9 1803000 44806 0.00236 5.94e+17 0 1 3.8e+15 2013 2.0e+4 1.7e+21 1.5e-10 5.5e+20 9.0e-5 99 50000000.0 1e+29 run
 ```

 ^ here, the 3 indicates that the data folder is specified and no prev file is made. 1 is the EIC model type. Then 0.34 is z (redshift) and then the rest of the parameters are enumerated.
 
## Changes
Changes to the code by Sarah Youngquist in Feb 2022:

  - Allow for _prev files not to be created.
        Specified in the command line with the mode.
        Stored in INPUT_MODE
        Now, any time \_prev\_*.dat file would be made, it is conditional on INPUT_MODE being 0 (user entered 0 or 2 for mode)
  - Allow for specifying the data folder
        This is specified in the command line with the mode.
        Note: relative path to any data file
  - Use the model without a parameter file, inputting parameters from the
    command line.
  - Function load_params_from_list(argv, model_type, index, argc) added.
  - Most of main method moved into run_models
  - In main method, parse the input mode and set the input type, set data directory
     if applicable, parse model type if applicable, call either load_params or load_params_from_list
