.. _reference:

.. currentmodule:: parcel_model

Reference
=========

Main Parcel Model
-----------------

The core of the model has its own documentation page, which you can access :ref:`here <parcel>`.

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

    dist_to_conc
    AerosolSpecies

Distributions
-------------

.. automodule:: parcel_model.distributions

.. autosummary::
    :toctree: generated/

    Lognorm
    MultiModeLognorm

Activation
----------

.. automodule:: parcel_model.activation

.. autosummary::
    :toctree: generated/

    act_fraction