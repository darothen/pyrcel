.. _install:

Installation
============

To grab and build the latest version of the model, you should use ``pip`` and point it
to the source code `repository`_ on github:

    $> pip install git+git://github.com/darothen/parcel_model.git

This should automatically build the necessary Cython modules and export the code
package to your normal package installation directory. If you wish to simply build
the code and run it in place, clone the `repository`_, navigate to it in a terminal,
and invoke the build command by hand:

    $> python setup.py build_ext --inplace

This should produce the compiled file `parcel_aux.so` in the model package. You can
also install the code from the cloned source directory by invoking ``pip install`` from
within it.


Dependencies
------------

This code was written from the Python 2.7 codebase. It probably won't work with Python
3.0+, so don't bother! By far, the simplest way to run this code is to grab a
scientific python distribution, such as
`Anaconda <https://store.continuum.io/cshop/anaconda/>`_. This code should work
out-of-the box with almost all dependencies filled (exception being numerical solvers)
on a recent version (2.0+) of this distribution.

Necessary dependencies
^^^^^^^^^^^^^^^^^^^^^^

- `numpy <http://www.numpy.org/>`_

- `scipy <http://www.scipy.org/>`_

- `pandas <http://pandas.pydata.org/>`_

- `Cython <http://cython.org/>`_

Numerical solver dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **LSODA** - `scipy <http://www.scipy.org/>`_ or
  `odespy <https://github.com/hplgit/odespy/>`_

- **VODE**, **LSODE** - `odespy <https://github.com/hplgit/odespy/>`_

- **CVODE** - `Assimulo <http://www.jmodelica.org/assimulo_home/index.html>`_

Recommended additional packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- `matplotlib <http://matplotlib.sourceforge.net/>`_

- `seaborn <http://stanford.edu/~mwaskom/software/seaborn/index.html>`_

Testing
-------

At the moment no testing package ships with the code. This will change in the future;
however, you may wish to navigate to the :ref:`basic run example <example_basic>` and
attempt to reproduce the result saved there to ensure the model code is working
correctly.

Bugs / Suggestions
------------------

The code has an
`issue tracker on github <https://github.com/darothen/parcel_model/issues>`_ and
I strongly encourage you to note any problems with the model there, such as typos
or weird behavior and results. Furthermore, I'm looking for ways to expand and
extend the model, so if there is something you might wish to see added, please
note it there or `send me an e-mail <mailto:darothen@mit.edu>`_. The code was written
in such a way that it should be trivial to add physics in a modular fashion.

.. _repository: http://github.com/darothen/parcel_model