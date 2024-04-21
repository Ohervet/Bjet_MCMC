
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
---------

- pip: 

  .. code-block:: console

    $ python -m pip install -U scipy

- conda:

  .. code-block:: console
  
    $ conda install -c conda-forge scipy

- `scipy` installation `documentation <https://scipy.org/install/>`_

``invoke``
----------
- pip: 

  .. code-block:: console
  
    $ python -m pip install invoke

- conda: 

  .. code-block:: console
  
    $ conda install -c conda-forge invoke


``matplotlib``
--------------
- conda: 

.. code-block:: console

  $ conda install matplotlib

- pip: 

.. code-block:: console

  $ python -m pip install -U matplotlib

- `matplotlib` installation `documentation <https://matplotlib.org/stable/users/installing/index.html>`_

``corner``
----------
- conda:

.. code-block:: console

  $ conda install -c astropy corner

- pip: 

.. code-block:: console

  $ python -m pip install corner
- `corner` installation `documentation <https://corner.readthedocs.io/en/latest/install.html>`_

``tqdm``
--------
- pip: 

.. code-block:: console

  $ python -m pip install tqdm

``h5py``
--------
- pip: 

.. code-block:: console

  $ python -m pip install h5py
