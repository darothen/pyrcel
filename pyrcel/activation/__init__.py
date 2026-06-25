"""Activation parameterizations.

JAX-native (differentiable) implementations live in this package.
Legacy NumPy/SciPy implementations are preserved in
`pyrcel.legacy.activation`; two names (`ming2006`, `shipwayabel2010`)
are re-exported here for backward compatibility.
"""

from ._arg2000 import ARG2000, arg2000  # noqa: F401
from ._common import (  # noqa: F401
    binned_activation,
    lognormal_activation,
    multi_mode_activation,
)
from ._mbn2014 import MBN2014, mbn2014  # noqa: F401
from ._scheme import ActivationScheme  # noqa: F401

# Legacy-only parameterizations re-exported for backward compatibility.
from ..legacy.activation import ming2006, shipwayabel2010  # noqa: F401
