============
Parcel Model
============

This is an implementation of a simple, adiabatic cloud parcel model for use in aerosol-cloud interaction studies. It is based on the model used by `Nenes et al. (2001) <http://nenes.eas.gatech.edu/Preprints/KinLimitations_TellusPP.pdf>`_, but with several key modifications:

* Implementation of :math:`\kappa`-Kohler theory for droplet physics
* Extension of model to handle arbitrary sectional representations of aerosol populations
* Improved, modular numerical framework for integrating the model

among other details. In addition, packaged with the model is a maturing library of warm-cloud droplet activation parameterizations, which will improve in breadth and documentation over time.

Updated code can be found the project GitHub `repository`_. If you'd like to use this code or have any questions about it, please contact the author at:

    Daniel Rothenberg <darothen@mit.edu>

.. _repository:
    https://github.com/darothen/parcel_model

-----------
Model Setup
-----------

No setup is necessarily required to run the model, although the main derivative function used by the numerical integrator has been Cythonized to improve performance. To use this optimized function, its module must be built by invoking::

    $ python setup.py build_ext --inplace

You may need to over-ride in flags ``setup.py`` depending on the OpenMP flag used by your compiler. Then, when the model runs, the number of threads used to calculate the derivative function is controlled by the environmental variable **OMP_NUM_THREADS**.

A suite of external utilties and functions for studying aerosol activation can be found in the `ext` folder. To run them, navigate to their location and execute the provided `Makefile`. If these components haven't been built, a warning will be passed whenever the parcel model module is imported.

------------
Dependencies
------------

Several basic scientific Python libraries are needed to run the model, including:

    1. NumPy
    2. SciPy
    3. Pandas
    4. Cython (optional)
    5. `odespy <https://github.com/hplgit/odespy>`_ (optional - unifed ODEPACK solver interface)
    6. `Assimulo <http://www.jmodelica.org/assimulo>`_ (optional - wrapper for `SUNDIALS CVODE <http://computation.llnl.gov/casc/sundials/main.html>`_ solver)

It's recommended to use a pre-packaged scientific Python distribution; this code is developed using the `Anaconda <https://store.continuum.io/cshop/anaconda/>`_ distribution, and executes on the most recent stable package.