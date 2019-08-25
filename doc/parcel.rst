.. _parcel:

.. currentmodule:: pyrcel

Parcel Model Details
====================

Below is the documentation for the parcel model, which is useful for debugging
and development. For a higher-level overview, see the :ref:`scientific description
<sci_descr>`.

Implementation
--------------

.. autoclass:: ParcelModel
    :members: set_initial_conditions, run

Derivative Equation
-------------------

.. automethod:: parcel.parcel_ode_sys