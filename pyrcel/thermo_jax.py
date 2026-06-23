"""JAX implementations of the aerosol/atmospheric thermodynamics functions.

This is the v2 (JAX + diffrax) counterpart to :mod:`pyrcel.thermo`. Every function
here is a faithful, elementwise translation of the corresponding NumPy function and
is intended to agree with it to within floating-point round-off (see
``tests/test_thermo_jax.py``). Unlike the NumPy version these functions are
``jax``-traceable: they accept scalars or arrays, broadcast naturally, and are
differentiable via ``jax.grad``/``jax.jacfwd``.

float64 is **mandatory** for this model (radii ~1e-8 m and the ``r**3 - r_dry**3``
cancellation in :func:`Seq`). We enable it at import time so the module is correct
regardless of how it is first used; this is idempotent and matches the design doc
(§6.1).
"""

from __future__ import annotations

import math

import jax

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402

from . import constants as c  # noqa: E402

PI = math.pi


def sigma_w(T):
    """Surface tension of water for a given temperature.

    Parameters
    ----------
    T : float
        Ambient temperature, K.

    Returns
    -------
    float
        Surface tension, J/m^2.

    See Also
    --------
    pyrcel.thermo.sigma_w
    """
    return 0.0761 - 1.55e-4 * (T - 273.15)


def dv_cont(T, P):
    """Diffusivity of water vapor in air, neglecting non-continuum effects.

    Parameters
    ----------
    T : float
        Ambient temperature, K.
    P : float
        Ambient pressure, Pa.

    Returns
    -------
    float
        Water vapor diffusivity, m^2/s.

    See Also
    --------
    pyrcel.thermo.dv_cont
    """
    P_atm = P * 1.01325e-5  # Pa -> atm
    return 1e-4 * (0.211 / P_atm) * ((T / 273.0) ** 1.94)


def dv(T, r, P, accom=c.ac):
    """Diffusivity of water vapor in air, modified for non-continuum effects.

    Parameters
    ----------
    T : float
        Ambient temperature, K.
    r : array or float
        Droplet/particle radius, m.
    P : float
        Ambient pressure, Pa.
    accom : float, optional
        Condensation coefficient (default :data:`pyrcel.constants.ac`).

    Returns
    -------
    array or float
        Non-continuum-corrected water vapor diffusivity, m^2/s.

    See Also
    --------
    pyrcel.thermo.dv
    """
    dv_t = dv_cont(T, P)
    denom = 1.0 + (dv_t / (accom * r)) * jnp.sqrt((2.0 * PI * c.Mw) / (c.R * T))
    return dv_t / denom


def es(T_c):
    """Saturation vapor pressure over water for a given temperature.

    Parameters
    ----------
    T_c : float
        Ambient temperature, degrees C.

    Returns
    -------
    float
        Saturation vapor pressure, Pa.

    See Also
    --------
    pyrcel.thermo.es
    """
    return 611.2 * jnp.exp(17.67 * T_c / (T_c + 243.5))


def ka_cont(T):
    """Thermal conductivity of air, neglecting non-continuum effects.

    Parameters
    ----------
    T : float
        Ambient temperature, K.

    Returns
    -------
    float
        Thermal conductivity, J/m/s/K.

    See Also
    --------
    pyrcel.thermo.ka_cont
    """
    return 1e-3 * (4.39 + 0.071 * T)


def ka(T, rho, r):
    """Thermal conductivity of air, modified for non-continuum effects.

    Parameters
    ----------
    T : float
        Ambient temperature, K.
    rho : float
        Ambient air density, kg/m^3.
    r : array or float
        Droplet/particle radius, m.

    Returns
    -------
    array or float
        Non-continuum-corrected thermal conductivity, J/m/s/K.

    See Also
    --------
    pyrcel.thermo.ka
    """
    ka_t = ka_cont(T)
    denom = 1.0 + (ka_t / (c.at * r * rho * c.Cp)) * jnp.sqrt(
        (2.0 * PI * c.Ma) / (c.R * T)
    )
    return ka_t / denom


def rho_air(T, P, RH=1.0):
    """Density of moist air for a given temperature, pressure, and relative humidity.

    Parameters
    ----------
    T : float
        Ambient temperature, K.
    P : float
        Ambient pressure, Pa.
    RH : float, optional
        Relative humidity, decimal (default 1.0).

    Returns
    -------
    float
        Air density, kg/m^3.

    See Also
    --------
    pyrcel.thermo.rho_air
    """
    qsat = RH * 0.622 * (es(T - 273.15) / P)
    Tv = T * (1.0 + 0.61 * qsat)
    return P / c.Rd / Tv


def Seq(r, r_dry, T, kappa):
    """κ-Köhler equilibrium supersaturation over an aerosol particle.

    Faithful to the numba derivative's :func:`Seq`, including the ``kappa > 0``
    guard (for ``kappa <= 0`` the curvature/solute term ``B`` collapses to 1).

    Parameters
    ----------
    r : array or float
        Droplet radius, m.
    r_dry : array or float
        Dry particle radius, m.
    T : float
        Ambient temperature, K.
    kappa : array or float
        Particle hygroscopicity parameter.

    Returns
    -------
    array or float
        Equilibrium supersaturation :math:`S_{eq}`.

    See Also
    --------
    pyrcel.thermo.Seq
    """
    A = (2.0 * c.Mw * sigma_w(T)) / (c.R * T * c.rho_w * r)
    r3 = r**3
    rd3 = r_dry**3
    B_full = (r3 - rd3) / (r3 - rd3 * (1.0 - kappa))
    B = jnp.where(kappa > 0.0, B_full, 1.0)
    return jnp.exp(A) * B - 1.0


def Seq_approx(r, r_dry, T, kappa):
    """Approximate κ-Köhler equilibrium supersaturation over an aerosol particle.

    Parameters
    ----------
    r : array or float
        Droplet radius, m.
    r_dry : array or float
        Dry particle radius, m.
    T : float
        Ambient temperature, K.
    kappa : array or float
        Particle hygroscopicity parameter.

    Returns
    -------
    array or float
        Approximate equilibrium supersaturation.

    See Also
    --------
    pyrcel.thermo.Seq_approx
    """
    A = (2.0 * c.Mw * sigma_w(T)) / (c.R * T * c.rho_w * r)
    return A - kappa * (r_dry**3) / (r**3)
