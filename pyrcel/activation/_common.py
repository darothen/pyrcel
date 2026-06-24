"""Shared utility functions for JAX-native activation parameterizations."""

from __future__ import annotations

import jax

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402
from jax import Array  # noqa: E402
from jax.scipy.special import erfc  # noqa: E402
from jax.typing import ArrayLike  # noqa: E402


def _lognormal_act(
    smax: ArrayLike, sigma: ArrayLike, N: ArrayLike, sgi: ArrayLike
) -> tuple[Array, Array]:
    """Activated fraction of one lognormal mode.

    Parameters
    ----------
    smax : float
        Maximum parcel supersaturation (decimal).
    sigma : float
        Geometric standard deviation.
    N : float
        Total number concentration (any unit; returned N_act is in the same unit).
    sgi : float
        Modal critical supersaturation (decimal).

    Returns
    -------
    N_act : float
        Activated number concentration.
    act_frac : float
        Activated fraction (0–1).
    """
    ui = 2.0 * jnp.log(sgi / smax) / (3.0 * jnp.sqrt(2.0) * jnp.log(sigma))
    N_act = 0.5 * N * erfc(ui)
    return N_act, N_act / N
