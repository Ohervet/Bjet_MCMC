
Dependencies
============

**Recommended:** create conda env from ``environment.yml`` using ``conda env create -f environment.yml``


``emcee``
---------

- conda: 

  .. code-block:: console
  
    $ conda install -c conda-forge emcee

- pip: 

``emcee`` recommends calling

.. code-block:: console

  $ pip install -U setuptools setuptools_scm pep517
  $ pip install -U emcee
- ``emcee`` installation `documentation <https://emcee.readthedocs.io/en/stable/user/install/>`_

``scipy``

- pip: 

.. code-block:: console
  $ python -m pip install -U scipy

- conda:

.. code-block:: console

  $ conda install -c conda-forge scipy

- `scipy` installation `documentation <https://scipy.org/install/>`_

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
