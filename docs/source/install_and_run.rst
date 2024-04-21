Installation and running Bjet_MCMC
==================================

.. _installation:

Installation
------------
You first need to clone the Bjet_MCMC Github repo in your local computer with

.. code-block:: console
   $ git clone https://github.com/Ohervet/Bjet_MCMC

Ensure all dependencies are installed, see [dependencies](#dependencies). Recommended: create the conda env "bjet-mcmc" from `environment.yml` using `conda env create -f environment.yml`. 

.. code-block:: console

   $ conda env create -f environment.yml

Then load your environment with

.. code-block:: console

   $ conda activate bjet-mcmc



Running Bjet_MCMC for the first time
------------------------------------

1. Create a copy of ``mcmc_config_template.txt`` called `mcmc_config.txt`. For information and customized configuration files, see [configuration file](#mcmc-configuration-file). 
2. Ensure the data file is formatted as described below and is at the relative path specified in ``mcmc_config.txt``, see [data format](#data-format).
3. Navigate inside folder ``blazar_mcmc``. Currently, code **must** be run from inside this folder or the imports won't work.
4. First time use: execute ``blazar_initialize.py`` with ``python blazar_initialize.py`` or ``python3 blazar_initialize.py`` depending on your Python setup. This creates all necessary folders and compiles the C++ code (see [bjet](#using-bjet_02)).
5. Execute ``blazar_run_mcmc.py`` with ``python blazar_run_mcmc.py`` or ``python3 blazar_run_mcmc.py`` depending on your Python setup. 

By default, configurations are read from ``mcmc_config.txt`` and results are written to a file in ``local_results`` entitled ``run_yyyy-mm-dd-hh:mm:ss``.

