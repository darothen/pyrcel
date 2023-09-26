"""
Adiabatic Cloud Parcel Model
----------------------------

This module implements a zero-dimensional, constant updraft
adiabatic cloud parcel model, suitable for studying aerosol effects
on droplet activation.

"""

from importlib.metadata import version as _version

try:
    __version__ = _version("pyrcel")
except Exception:
    # This is a local copy, or a copy that was not installed via setuptools
    __version__ = "local"

__author__ = "Daniel Rothenberg <daniel@danielrothenberg.com>"

# TODO: Re-factor module-wide implicit imports
from .activation import *
from .aerosol import *
from .distributions import *
from .driver import *
from .parcel import *
from .thermo import *
