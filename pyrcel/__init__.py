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

# Lightweight modules (NumPy/SciPy only) are imported eagerly.
from .activation import *
from .aerosol import *
from .distributions import *
from .thermo import *

# The legacy numba/Assimulo-backed entry points (``ParcelModel`` and the driver
# helpers) are imported lazily. This lets the v2 JAX modules
# (``pyrcel.thermo_jax``, ``pyrcel.parcel_aux_jax``) and the lightweight utilities
# above be used without importing numba — a prerequisite for the JAX/diffrax
# migration, where these legacy backends are being removed entirely.
_LAZY_ATTRS = {
    "ParcelModel": "pyrcel.parcel",
    "run_model": "pyrcel.driver",
    "iterate_runs": "pyrcel.driver",
}


def __getattr__(name):
    module_path = _LAZY_ATTRS.get(name)
    if module_path is not None:
        import importlib

        module = importlib.import_module(module_path)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + list(_LAZY_ATTRS))
