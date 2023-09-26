.. _reference:

.. currentmodule:: pyrcel

Reference
=========

Main Parcel Model
-----------------

The core of the model has its own documentation page, which you can access :ref:`here <parcel>`.

.. autosummary::
    :toctree: generated/

    ParcelModel

Driver Tools
------------

.. automodule:: pyrcel.driver

.. autosummary::
    :toctree: generated/

    run_model
    iterate_runs

Thermodynamics/Kohler Theory
----------------------------

.. automodule:: pyrcel.thermo

.. autosummary::
    :toctree: generated/

    dv
    ka
    rho_air
    es
    sigma_w
    Seq
    Seq_approx
    kohler_crit
    critical_curve

Aerosols
--------

.. automodule:: pyrcel.aerosol

.. autosummary::
    :toctree: generated/
    :template: class.rst

    AerosolSpecies

The following are utility functions which might be useful in studying
and manipulating aerosol distributions for use in the :ref:`model <parcel>`
or activation routines.

.. autosummary::
    :toctree: generated/

    dist_to_conc

Distributions
-------------

.. automodule:: pyrcel.distributions

.. autosummary::
    :toctree: generated/
    :template: class.rst

    BaseDistribution
    Gamma
    Lognorm
    MultiModeLognorm

The following dictionaries containing (multi) Lognormal aerosol size distributions have also been saved for convenience:

1. ``FN2005_single_modes``: Fountoukis, C., and A. Nenes (2005), Continued development of a cloud droplet formation parameterization for global climate models, J. Geophys. Res., 110, D11212, doi:10.1029/2004JD005591
2. ``NS2003_single_modes``: Nenes, A., and J. H. Seinfeld (2003), Parameterization of cloud droplet formation in global climate models, J. Geophys. Res., 108, 4415, doi:10.1029/2002JD002911, D14.
3. ``whitby_distributions``: Whitby, K. T. (1978), The physical characteristics of sulfur aerosols, Atmos. Environ., 12(1-3), 135–159, doi:10.1016/0004-6981(78)90196-8.
4. ``jaenicke_distributions``: Jaenicke, R. (1993), Tropospheric Aerosols, in *Aerosol-Cloud-Climate Interactions*, P. V. Hobbs, ed., Academic Press, San Diego, CA, pp. 1-31.


Activation
----------

.. automodule:: pyrcel.activation

.. autosummary::
    :toctree: generated/

    lognormal_activation
    binned_activation
    multi_mode_activation
    arg2000
    mbn2014
    shipwayabel2010
    ming2006

.. _constants:

Constants
---------

.. automodule:: pyrcel.constants