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
from jax import Array  # noqa: E402
from jax.typing import ArrayLike  # noqa: E402

from . import constants as c  # noqa: E402

PI = math.pi


def sigma_w(T: ArrayLike) -> Array:
    """Surface tension of water for a given temperature.

    Parameters
    ----------
    T : float
        Ambient temperature, K.

    Returns
    -------
    float
        Surface tension, J/m².

    Notes
    -----
    Linear fit to experimental data:

    $$\sigma_w(T) = 0.0761 - 1.55 \times 10^{-4}(T - 273.15)$$

    See Also
    --------
    pyrcel.legacy.thermo.sigma_w
    """
    return 0.0761 - 1.55e-4 * (T - 273.15)


def dv_cont(T: ArrayLike, P: ArrayLike) -> Array:
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
        Water vapor diffusivity, m²/s.

    Notes
    -----
    $$D_v(T, P) = 10^{-4} \cdot \frac{0.211}{P_\mathrm{atm}}
      \left(\frac{T}{273}\right)^{1.94}$$

    where $P_\mathrm{atm}$ is pressure in atmospheres.

    See Also
    --------
    pyrcel.legacy.thermo.dv_cont
    """
    P_atm = P * 1.01325e-5  # Pa -> atm
    return 1e-4 * (0.211 / P_atm) * ((T / 273.0) ** 1.94)


def dv(T: ArrayLike, r: ArrayLike, P: ArrayLike, accom: ArrayLike = c.ac) -> Array:
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
        Non-continuum-corrected water vapor diffusivity, m²/s.

    Notes
    -----
    $$D_v'(T, r, P) = \frac{D_v(T,P)}{1 + \dfrac{D_v(T,P)}{\alpha_c\, r}
      \sqrt{\dfrac{2\pi M_w}{R T}}}$$

    where $\alpha_c$ is the condensation (accommodation) coefficient.

    See Also
    --------
    pyrcel.legacy.thermo.dv
    """
    dv_t = dv_cont(T, P)
    denom = 1.0 + (dv_t / (accom * r)) * jnp.sqrt((2.0 * PI * c.Mw) / (c.R * T))
    return dv_t / denom


def es(T_c: ArrayLike) -> Array:
    """Saturation vapor pressure over water for a given temperature.

    Parameters
    ----------
    T_c : float
        Ambient temperature, °C.

    Returns
    -------
    float
        Saturation vapor pressure, Pa.

    Notes
    -----
    Magnus formula:

    $$e_s(T_c) = 611.2 \exp\!\left(\frac{17.67\, T_c}{T_c + 243.5}\right)$$

    See Also
    --------
    pyrcel.legacy.thermo.es
    """
    return 611.2 * jnp.exp(17.67 * T_c / (T_c + 243.5))


def ka_cont(T: ArrayLike) -> Array:
    """Thermal conductivity of air, neglecting non-continuum effects.

    Parameters
    ----------
    T : float
        Ambient temperature, K.

    Returns
    -------
    float
        Thermal conductivity, J m⁻¹ s⁻¹ K⁻¹.

    Notes
    -----
    $$k_a(T) = 10^{-3}(4.39 + 0.071\, T)$$

    See Also
    --------
    pyrcel.legacy.thermo.ka_cont
    """
    return 1e-3 * (4.39 + 0.071 * T)


def ka(T: ArrayLike, rho: ArrayLike, r: ArrayLike) -> Array:
    """Thermal conductivity of air, modified for non-continuum effects.

    Parameters
    ----------
    T : float
        Ambient temperature, K.
    rho : float
        Ambient air density, kg/m³.
    r : array or float
        Droplet/particle radius, m.

    Returns
    -------
    array or float
        Non-continuum-corrected thermal conductivity, J m⁻¹ s⁻¹ K⁻¹.

    Notes
    -----
    $$k_a'(T, \rho, r) = \frac{k_a(T)}{1 + \dfrac{k_a(T)}{\alpha_t\, r\, \rho\, C_p}
      \sqrt{\dfrac{2\pi M_a}{R T}}}$$

    where $\alpha_t$ is the thermal accommodation coefficient.

    See Also
    --------
    pyrcel.legacy.thermo.ka
    """
    ka_t = ka_cont(T)
    denom = 1.0 + (ka_t / (c.at * r * rho * c.Cp)) * jnp.sqrt((2.0 * PI * c.Ma) / (c.R * T))
    return ka_t / denom


def rho_air(T: ArrayLike, P: ArrayLike, RH: ArrayLike = 1.0) -> Array:
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
        Air density, kg/m³.

    Notes
    -----
    Uses the virtual temperature $T_v = T(1 + 0.61\, q_\mathrm{sat})$ with

    $$q_\mathrm{sat} = \mathrm{RH} \cdot 0.622 \cdot \frac{e_s(T)}{P}$$

    so that

    $$\rho = \frac{P}{R_d\, T_v}$$

    See Also
    --------
    pyrcel.legacy.thermo.rho_air
    """
    qsat = RH * 0.622 * (es(T - 273.15) / P)
    Tv = T * (1.0 + 0.61 * qsat)
    return P / c.Rd / Tv


def Seq(r: ArrayLike, r_dry: ArrayLike, T: ArrayLike, kappa: ArrayLike) -> Array:
    """κ-Köhler equilibrium supersaturation over an aerosol particle.

    Two numerical stability improvements over the naïve formulation:

    1. **Difference-of-cubes**: ``r³ - r_dry³`` is rewritten as
       ``(r - r_dry)(r² + r·r_dry + r_dry²)`` to avoid catastrophic
       cancellation when ``r ≈ r_dry`` (near the dry particle limit).
    2. **Compensated final step**: ``exp(A)·B - 1`` is replaced by
       ``B·expm1(A) + (B - 1)`` where ``B - 1 = -κ·r_dry³ / denom`` is
       computed without cancellation, so the result is accurate even when
       ``A ≪ 1`` (large droplets) and ``B ≈ 1`` (dilute solution).

    For ``kappa <= 0`` the solute term is suppressed (``B = 1``,
    ``B - 1 = 0``), so the result reduces to the pure-curvature Kelvin
    term ``expm1(A)``.

    Parameters
    ----------
    r : array or float
        Droplet radius, m.
    r_dry : array or float
        Dry particle radius, m.
    T : float
        Ambient temperature, K.
    kappa : array or float
        Particle hygroscopicity parameter.  For ``kappa <= 0`` the solute
        term is zeroed out and only the Kelvin curvature term remains.

    Returns
    -------
    array or float
        Equilibrium supersaturation $S_\mathrm{eq}$.

    Notes
    -----
    The full κ-Köhler equation [Petters2007]_:

    $$S_\mathrm{eq}(r) = \exp\!\left(\frac{A}{r}\right)
      \frac{r^3 - r_d^3}{r^3 - r_d^3(1 - \kappa)} - 1$$

    where the Kelvin parameter is

    $$A = \frac{2 M_w \sigma_w(T)}{R T \rho_w}$$

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


def Seq_approx(r: ArrayLike, r_dry: ArrayLike, T: ArrayLike, kappa: ArrayLike) -> Array:
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

    Notes
    -----
    First-order expansion of :func:`Seq` valid when $A/r \ll 1$ and
    $\kappa r_d^3 / r^3 \ll 1$:

    $$S_\mathrm{eq}(r) \approx \frac{A}{r} - \kappa \frac{r_d^3}{r^3}$$

    Used for the critical-radius approximation in :mod:`pyrcel.equilibrate`.

    See Also
    --------
    pyrcel.legacy.thermo.Seq_approx
    """
    A = (2.0 * c.Mw * sigma_w(T)) / (c.R * T * c.rho_w * r)
    return A - kappa * (r_dry**3) / (r**3)
