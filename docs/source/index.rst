Welcome to Bjet_MCMC
===================================

**Bjet_MCMC** is a tool to automatically model multiwavelength spectral energy distributions of blazars, considering one-zone synchrotron-self-Compton (SSC) model with or without the addition of external inverse-Compton process from the thermal emission of the nucleus. It also contains manual fitting functionalities for multi-zone SSC modeling.
This tool is built as an MCMC python wrapper around the C++ code Bjet.

License
-------
The code is licensed under a :doc:`BSD 3-Clause License <License>`.


Acknowledgments and citation
-------
- If you use this code in a scientific publication or communication, please cite the paper `Hervet et al. 2024 <https://iopscience.iop.org/article/10.3847/1538-4357/ad09c0>`_ .
- For any use of the multi-zones SSC option, please also cite `Hervet at al. 2015 <https://ui.adsabs.harvard.edu/abs/2015A%26A...578A..69H/abstract>`_ .
- To reference a specific version, you can also use Zenodo citation: 

  .. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.10070356.svg
    :target: https://doi.org/10.5281/zenodo.10070356



.. note::

   This project is under active development.

Contents
--------
.. toctree::

   install_and_run
   dependencies
   configuration_file
   data_format
   outputs
   bjet_core
