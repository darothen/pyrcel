"""Activation parameterizations.

JAX-native (differentiable) implementations live in this package.
Legacy NumPy/SciPy implementations are preserved in
:mod:`pyrcel.legacy.activation` and re-exported here for backward
compatibility.

JAX-native
----------
:func:`arg2000`
    Abdul-Razzak & Ghan (2000) — closed-form, fully differentiable.
:class:`ARG2000`
    Callable class wrapper around :func:`arg2000`.
:class:`ActivationScheme`
    Abstract base class for all activation schemes.

Legacy (NumPy/SciPy)
--------------------
The following names are re-exported from :mod:`pyrcel.legacy.activation`
so that code importing from ``pyrcel.activation`` continues to work.
``mbn2014``, ``binned_activation``, ``lognormal_activation``,
``multi_mode_activation``, ``ming2006``, ``shipwayabel2010``.
"""

# Legacy re-exports — keep the original names accessible.
from ..legacy.activation import (  # noqa: F401
    binned_activation,
    lognormal_activation,
    mbn2014,
    ming2006,
    multi_mode_activation,
    shipwayabel2010,
)
from ._arg2000 import ARG2000, arg2000  # noqa: F401
from ._scheme import ActivationScheme  # noqa: F401
