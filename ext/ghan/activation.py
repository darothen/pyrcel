"""

.. module:: parcel
    :synopsis: Python wrapper for calling Steve Ghan's model.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""
__docformat__ = 'reStructuredText'

import numpy as np
import pandas as pd

from ghan_interface import ghan_model

def _unpack_aerosols(aerosols):
    """Unpackages aerosols used in the Python parcel model to be used in the
    Ghan code.

    Looks at each :class:`AerosolSpecies` in ``aerosols``, and compiles lists of the
    following aerosol properties:

        * number concentration
        * distribution width parameter
        * hygroscopicity of species
        * geometric mean size
        * species dry density

    If the density is not specified for a given aerosol (contained as an attribute of
    the :class:`AerosolSpecies` object representing it), then a default of 1 g/m^3 is
    assumed.

    **Args**:
        * *aerosols* -- a list of :class:`AerosolSpecies` objects

    """

    ## Unpack the aerosol data for the Ghan code
    na, sig, rhodry, hygro, rad = [], [], [], [], []
    for aerosol in aerosols:
        na.append(aerosol.N)
        sig.append(aerosol.sigma)
        hygro.append(aerosol.kappa)
        rad.append(aerosol.mu)

        if hasattr(aerosol, 'rho'):
            rhodry.append(aerosol.rho*1e-3)
        else:
            rhodry.append(1.)

    return na, sig, rhodry, hygro, rad

def explicit(aerosols, V, T0, S0, P0):
    """Run Steve Ghan's explicit parcel model.

    """
    na, sig, rhodry, hygro, rad = _unpack_aerosols(aerosols)

    ## Convert units
    P0 = P0/100. # Pa -> hPa
    V = V*100. # m/s -> cm/s

    ## Call Ghan model code
    explicit = True
    fn, smax = ghan_model(T0, P0, V, explicit, na, sig, rhodry, hygro, rad)

    ## TODO: It would be nice to spit out a trace of the Ghan model vertical profile,
    ## but this is hard. Could either re-write part of subgrid_mod.f90 (NO), or try to
    ## read the output file "fort.100". Problem with the latter is that the file doesn't
    ## finalize until after Python has exited.

    return fn, smax

def parameterization(aerosols, V, T0, S0, P0):
    """Run the provided parameterization of the Ghan parcel model.

    """
    na, sig, rhodry, hygro, rad = _unpack_aerosols(aerosols)

    ## Convert units
    P0 = P0/100. # Pa -> hPa
    V = V*100. # m/s -> cm/s

    ## Call Ghan model code
    explicit = False
    fn, smax = ghan_model(T0, P0, V, explicit, na, sig, rhodry, hygro, rad)

    return fn, smax
