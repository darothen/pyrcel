"""Activation parameterizations — shim for PR #33.

This module re-exports the NumPy/SciPy implementations from :mod:`pyrcel.legacy.activation`
while the JAX-native activation subpackage (``pyrcel/activation/``) is being built in #36.
After #36 merges, this file will be replaced by ``pyrcel/activation/__init__.py``.
"""

from .legacy.activation import (  # noqa: F401
    arg2000,
    binned_activation,
    lognormal_activation,
    mbn2014,
    ming2006,
    multi_mode_activation,
    shipwayabel2010,
)
