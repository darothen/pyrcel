Cloud parcel model
==================

This is an implementation of a simple, 0D adiabatic cloud parcel model tool (following `Nenes et al, 2001`_ and `Pruppacher and Klett, 1997`_). It allows flexible descriptions of an initial aerosol population, and simulates the evolution of a proto-cloud droplet population as the parcel ascends adiabatically at either a constant or time/height-dependent updraft speed. Droplet growth within the parcel is tracked on a Lagrangian grid.

.. _Pruppacher and Klett, 1997: http://books.google.com/books?hl=en&lr=&id=1mXN_qZ5sNUC&oi=fnd&pg=PR15&ots=KhdkC6uhB3&sig=PSlNsCeLSB2FvR93Vzo0ptCAnYA#v=onepage&q&f=false
.. _Nenes et al, 2001: http://onlinelibrary.wiley.com/doi/10.1034/j.1600-0889.2001.d01-12.x/abstract

.. image:: figs/model_example.png

You are invited to use the model (in accordance with the `licensing <https://raw.githubusercontent.com/darothen/parcel_model/master/LICENSE>`_) as long as you get in touch with the author via `e-mail <mailto:darothen@mit.edu>`_ or on `twitter <https://twitter.com/darothen>`_. Up-to-date versions can be obtained through the model's `github repository <https://github.com/darothen/parcel_model>`_ or directly from the author. 

Documentation Outline
---------------------

.. toctree::
    :maxdepth: 2

    sci_descr
    install
    examples/basic_run
    parcel
    reference


Current version: |version|

Documentation last compiled: |today|
