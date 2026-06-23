"""JAX implementations of the aerosol/atmospheric thermodynamics functions.

This is the v2 (JAX + diffrax) primary thermodynamics module. Every function
here is a faithful, elementwise translation of the corresponding NumPy reference
in :mod:`pyrcel.legacy.thermo` and is intended to agree with it to within
floating-point round-off (see ``tests/test_thermo.py``). Unlike the NumPy
version these functions are ``jax``-traceable: they accept scalars or arrays,
broadcast naturally, and are differentiable via ``jax.grad``/``jax.jacfwd``.

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
    pyrcel.legacy.thermo.sigma_w
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
    pyrcel.legacy.thermo.dv_cont
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
    pyrcel.legacy.thermo.dv
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
    pyrcel.legacy.thermo.es
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
    pyrcel.legacy.thermo.ka_cont
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
    pyrcel.legacy.thermo.ka
    """
    ka_t = ka_cont(T)
    denom = 1.0 + (ka_t / (c.at * r * rho * c.Cp)) * jnp.sqrt((2.0 * PI * c.Ma) / (c.R * T))
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
    pyrcel.legacy.thermo.rho_air
    """
    qsat = RH * 0.622 * (es(T - 273.15) / P)
    Tv = T * (1.0 + 0.61 * qsat)
    return P / c.Rd / Tv


def Seq(r, r_dry, T, kappa):
    """κ-Köhler equilibrium supersaturation over an aerosol particle.

    Two numerical stability improvements over the naïve formulation:

    1. **Difference-of-cubes**: ``r³ - r_dry³`` is rewritten as
       ``(r - r_dry)(r² + r·r_dry + r_dry²)`` to avoid catastrophic
       cancellation when ``r ≈ r_dry`` (near the dry particle limit).
    2. **Compensated final step**: ``exp(A)·B - 1`` is replaced by
       ``B·expm1(A) + (B - 1)`` where ``B - 1 = -κ·r_dry³ / denom`` is
       computed without cancellation, so the result is accurate even when
       ``A ≪ 1`` (large droplets) and ``B ≈ 1`` (dilute solution).

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
    pyrcel.legacy.thermo.Seq
    """
    A = (2.0 * c.Mw * sigma_w(T)) / (c.R * T * c.rho_w * r)

    # r³ - r_dry³ = (r - r_dry)(r² + r·r_dry + r_dry²); accurate near r ≈ r_dry.
    delta = r - r_dry
    sum_sq = r**2 + r * r_dry + r_dry**2
    num = delta * sum_sq  # r³ - r_dry³
    krd3 = kappa * r_dry**3
    denom = num + krd3  # r³ - r_dry³·(1 - κ)
    B_full = num / denom
    Bm1_full = -krd3 / denom  # B - 1, avoids cancellation when B ≈ 1

    B = jnp.where(kappa > 0.0, B_full, 1.0)
    Bm1 = jnp.where(kappa > 0.0, Bm1_full, 0.0)

    # exp(A)·B - 1 = B·expm1(A) + (B - 1); accurate when A ≪ 1.
    return B * jnp.expm1(A) + Bm1


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
    pyrcel.legacy.thermo.Seq_approx
    """
    A = (2.0 * c.Mw * sigma_w(T)) / (c.R * T * c.rho_w * r)
    return A - kappa * (r_dry**3) / (r**3)
