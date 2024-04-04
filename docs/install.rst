.. _install:

Quick Installation
==================


To quickly get started with running pyrcel, complete the following steps:

- Set up a new Python environment; we recommend using `mambaforge <https://conda-forge.org/miniforge/>`_:

.. code-block:: shell

    $ mamba create -n pyrcel_quick_start python=3.11

- Activate the new Python environment and install the model and its dependencies. If you install the published version from PyPi (recommended), then you also need to install `Assimulo <http://www.jmodelica.org/assimulo>`_ using the Mamba package manager - but no other manual dependency installation is necessary:

.. code-block:: shell

    $ mamba activate pyrcel_quick_start
    $ pip install pyrcel
    $ mamba install -c conda-forge assimulo

- Run a test simulation using the CLI tool and a sample YAML file from **pyrcel/examples/*.yml** (you may want to clone the repository or download them locally):

.. code-block:: shell

    $ run_parcel simple.yml


Detailed Installation Notes
===========================


From PyPI
---------

If you already have all the dependencies satisfied, then you can install the
latest release from `PyPI <https://badge.fury.io/py/pyrcel>`_ by using
``pip``:

.. code-block:: bash

    $ pip install pyrcel


From source code
----------------

To grab and build the latest bleeding-edge version of the model, you should use
``pip`` and point it to the source code `repository`_ on github:


.. code-block:: bash

    $ pip install git+git://github.com/darothen/pyrcel.git

This should automatically build the necessary Cython modules and export the
code package to your normal package installation directory. If you wish to
simply build the code and run it in place, clone the `repository`_, navigate
to it in a terminal, and invoke the build command by hand:


.. code-block:: bash

    $ python setup.py build_ext --inplace

This should produce the compiled file `parcel_aux.so` in the model package.
You can also install the code from the cloned source directory by invoking
``pip install`` from within it; this is useful if you're updating or
modifying the model, since you can install an "editable" package which
points directly to the git-monitored code:


.. code-block:: bash

    $ cd path/to/pyrcel/
    $ pip install -e .


Dependencies
------------

This code was originally written for Python 2.7, and then
`futurized <http://python-future.org/>`_ to Python 3.3+ with hooks for
backwards compatibility. It should work on modern Python versions, and we recommend
using Python 3.11+ for the greatest compatibility with required dependencies.

The easiest way to manage dependencies is to use a tool like `Mambaforge <https://conda-forge.org/miniforge/>`
to set up an environment. Suitable environment files can be found in the ``pyrcel/ci``
directory.

Necessary dependencies
^^^^^^^^^^^^^^^^^^^^^^

All of these (except for Assimulo; see the note below) can be installed via `pip`:

- `Assimulo <http://www.jmodelica.org/assimulo_home/index.html>`_

- `numba <http://numba.pydata.org>`_

- `numpy <http://www.numpy.org/>`_

- `scipy <http://www.scipy.org/>`_

- `pandas <http://pandas.pydata.org/>`_

.. note::

    As of version 1.2.0, the model integration components are being re-written
    and only the CVODE interface is exposed. As such, Assimulo is temporarily
    a core and required dependency; in the future the other solvers will
    be re-enabled. You should first try to install Assimulo via conda

    .. code-block:: bash

        $ mamba install -c conda-forge assimulo

    since this will automatically take care of obtaining necessary compiled
    dependencies like sundials. However, for best results you may want to
    `manually install Assimulo <http://www.jmodelica.org/assimulo_home/installation.html>`_,
    since the conda-forge recipe may default to a sundials/OpenBLAS combination
    which could degare the performance of the model.

Numerical solver dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **LSODA** - `scipy <http://www.scipy.org/>`_ or
  `odespy <https://github.com/hplgit/odespy/>`_

- **VODE**, **LSODE** - `odespy <https://github.com/hplgit/odespy/>`_

- **CVODE** - `Assimulo <http://www.jmodelica.org/assimulo_home/index.html>`_

Recommended additional packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

    These are not required for the model to run, but are useful for
    post-processing and visualization of the model output. They should be installed
    automatically if you install the model from PyPI or the source code repository.

- `matplotlib <http://matplotlib.sourceforge.net/>`_

- `seaborn <http://stanford.edu/~mwaskom/software/seaborn/index.html>`_

- `PyYAML <http://pyyaml.org/wiki/PyYAMLDocumentation>`_

- `xarray <http://xarray.pydata.org/en/stable/>`_

Testing
-------

A nose test-suite is under construction. To check that your model is configured
and running correctly, you copy and run the notebook corresponding to the
:ref:`basic run example <example_basic>`, or run the command-line interface
version of the model with the pre-packed simple run case:

.. code-block:: bash

    $ cd path/to/pyrcel/
    $ ./run_parcel examples/simple.yml


Bugs / Suggestions
------------------

The code has an
`issue tracker on github <https://github.com/darothen/pyrcel/issues>`_
and I strongly encourage you to note any problems with the model there, such
as typos or weird behavior and results. Furthermore, I'm looking for ways to
expand and extend the model, so if there is something you might wish to see
added, please note it there or `send me an e-mail <mailto:daniel@danielrothenberg.com>`_.
The code was written in such a way that it should be trivial to add physics in a modular fashion.

.. _repository: http://github.com/darothen/pyrcel
