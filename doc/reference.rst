.. _reference:

.. currentmodule:: parcel_model

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

.. automodule:: parcel_model.driver

.. autosummary::
    :toctree: generated/

    run_model
    iterate_runs

Thermodynamics/Kohler Theory
----------------------------

.. automodule:: parcel_model.thermo

.. autosummary::
    :toctree: generated/

    dv
    ka
    rho_air
    es
    sigma_w
    Seq
    kohler_crit
    critical_curve

Aerosols
--------

.. automodule:: parcel_model.aerosol

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

.. automodule:: parcel_model.distributions

.. autosummary::
    :toctree: generated/
    :template: class.rst

    Lognorm

Activation
----------

.. automodule:: parcel_model.activation

.. autosummary::
    :toctree: generated/

    multi_mode_activation
    act_fraction

Constants
---------

.. automodule:: parcel_model.constants