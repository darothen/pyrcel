"""Shared type aliases for pyrcel's public and internal APIs."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

# Accepts a Python scalar or a NumPy array of any floating dtype.
# Used for distribution pdf/cdf inputs and outputs that may be either.
FloatND = float | NDArray[np.floating[Any]]
