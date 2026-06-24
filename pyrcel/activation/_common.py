"""Shared utility functions for JAX-native activation parameterizations."""

from __future__ import annotations

import jax

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402
from jax import Array  # noqa: E402
from jax.scipy.special import erfc  # noqa: E402
from jax.typing import ArrayLike  # noqa: E402

from .. import constants as c  # noqa: E402
from ..thermo import sigma_w  # noqa: E402


def _lognormal_act(
    smax: ArrayLike, sigma: ArrayLike, N: ArrayLike, sgi: ArrayLike
) -> tuple[Array, Array]:
    r"""Activated fraction of one lognormal mode.

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

    Notes
    -----
    Integrates the lognormal number distribution above the critical dry size
    corresponding to $S_\mathrm{max}$:

    $$N_\mathrm{act} = \frac{N}{2}\,\mathrm{erfc}(u_i), \qquad
      u_i = \frac{2\ln(s_{g,i}/S_\mathrm{max})}{3\sqrt{2}\ln\sigma}$$

    where $s_{g,i}$ is the modal critical supersaturation.
    """
    ui = 2.0 * jnp.log(sgi / smax) / (3.0 * jnp.sqrt(2.0) * jnp.log(sigma))
    N_act = 0.5 * N * erfc(ui)
    return N_act, N_act / N


def _kohler_crit_approx(T: ArrayLike, r_dry: ArrayLike, kappa: ArrayLike) -> tuple[Array, Array]:
    r"""Approximate critical radius and supersaturation (JAX-traceable, kappa > 0).

    Vectorises over ``r_dry`` / ``kappa`` for scalar ``T``.

    Returns
    -------
    r_crit : Array
        Critical wet radius (m).
    s_crit : Array
        Critical supersaturation (decimal).

    Notes
    -----
    From the approximate κ-Köhler equation $S_\mathrm{eq} \approx A/r -
    \kappa r_d^3/r^3$, setting $dS_\mathrm{eq}/dr = 0$ gives:

    $$r_\mathrm{crit} = \sqrt{\frac{3\kappa r_d^3}{A}}, \qquad
      s_\mathrm{crit} = \sqrt{\frac{4A^3}{27\kappa r_d^3}}$$

    where $A = 2 M_w \sigma_w(T) / (\rho_w R T)$ is the Kelvin parameter.
    """
    A = (2.0 * c.Mw * sigma_w(T)) / (c.rho_w * c.R * T)
    r_crit = jnp.sqrt((3.0 * kappa * r_dry**3) / A)
    s_crit = jnp.sqrt((4.0 * A**3) / (27.0 * kappa * r_dry**3))
    return r_crit, s_crit
