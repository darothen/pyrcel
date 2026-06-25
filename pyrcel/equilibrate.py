"""JAX/optimistix equilibration: initial wet radii and parcel state vector.

This is the v2 counterpart to ``ParcelModel._setup_run`` in ``pyrcel/parcel.py``.
For each aerosol bin it finds the equilibrium wet radius ``r0`` -- the root of

    Seq(r, r_dry, T0, kappa) - S0 = 0

on the rising (stable) branch of the Köhler curve -- and assembles the initial
state vector ``y0 = [z, P, T, wv, wc, wi, S, r0_0 ... r0_{nr-1}]``.

Differences from ``master`` (and why they don't matter for fidelity):

* ``scipy.optimize.bisect`` -> [optimistix.Bisection][] (same bracketing
  method, kept in the JAX ecosystem; ``vmap``-ed over bins).
* The upper bracket is the **analytic** approximate Köhler critical radius rather
  than ``scipy.optimize.fminbound`` on ``-Seq``. The root ``r0`` is independent of
  the bracket as long as the bracket is valid (``Seq(r_dry) < S0 < Seq(r_crit)``),
  which holds for all sub-critical ``S0`` -- so ``r0`` still matches ``master`` to
  solver tolerance (see ``tests/test_equilibration.py``, design §7.2b).

Per the locked decision, gradients are required through the *integration* w.r.t.
``y0``, not through this solve; equilibration may run eagerly. (As a bonus,
``optimistix`` root-finds are differentiable via the implicit function theorem.)
"""

from __future__ import annotations

from typing import Any

import jax

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402
import optimistix as optx  # noqa: E402
from jax import Array  # noqa: E402
from jax.typing import ArrayLike  # noqa: E402

from . import constants as c  # noqa: E402
from .thermo import Seq, es, rho_air, sigma_w  # noqa: E402

#: Default root-find tolerances. ``master`` uses ``bisect(xtol=1e-30)`` (machine
#: precision); these reproduce ``r0`` to ~1e-13 relative, well inside §7.2b's 1e-10.
_RTOL = 1e-13
_ATOL = 1e-30
_MAX_STEPS = 200


def kohler_crit_approx(T: ArrayLike, r_dry: ArrayLike, kappa: ArrayLike) -> tuple[Array, Array]:
    """Analytic approximate Köhler critical radius and supersaturation.

    Mirrors `pyrcel.legacy.thermo.kohler_crit` with ``approx=True``. This
    approximation can return ``r_crit < r_dry`` for very small, low-κ particles,
    so it is **not** used to bracket the equilibrium root — `kohler_crit`
    is. Kept for reference and the analytic ``s_crit``.

    Parameters
    ----------
    T : float
        Ambient temperature, K.
    r_dry : float
        Dry particle radius, m.
    kappa : float
        Particle hygroscopicity parameter.

    Returns
    -------
    r_crit : float
        Approximate critical wet radius, m.
    s_crit : float
        Approximate critical supersaturation.

    See Also
    --------
    kohler_crit : Exact numerical Köhler critical radius.
    pyrcel.legacy.thermo.kohler_crit : NumPy equivalent with ``approx=True``.
    """
    A = (2.0 * c.Mw * sigma_w(T)) / (c.R * T * c.rho_w)
    r_crit = jnp.sqrt((3.0 * kappa * (r_dry**3)) / A)
    s_crit = jnp.sqrt((4.0 * (A**3)) / (27.0 * kappa * (r_dry**3)))
    return r_crit, s_crit


def kohler_crit(
    T: ArrayLike,
    r_dry: ArrayLike,
    kappa: ArrayLike,
    *,
    rtol: float = _RTOL,
    atol: float = _ATOL,
    max_steps: int = _MAX_STEPS,
) -> Array:
    """Exact Köhler critical (peak) radius via a root-find on ``dSeq/dr = 0``.

    The JAX analog of ``pyrcel.legacy.thermo.kohler_crit`` (which uses
    ``scipy.optimize.fminbound`` on ``-Seq``). ``Seq`` is unimodal for ``kappa > 0``,
    so its derivative -- obtained analytically with ``jax.grad`` -- has a single
    zero (the peak) on ``[r_dry, r_dry * 1e4]``. Always returns ``r_crit > r_dry``,
    making it a valid upper bracket for the equilibrium solve at any particle size.

    Parameters
    ----------
    T : float
        Ambient temperature, K.
    r_dry : float
        Dry particle radius, m.
    kappa : float
        Particle hygroscopicity parameter.
    rtol, atol : float, optional
        Root-find tolerances.
    max_steps : int, optional
        Maximum bisection iterations.

    Returns
    -------
    float
        Critical wet radius, m. Always satisfies ``r_crit > r_dry``.

    See Also
    --------
    kohler_crit_approx : Analytic approximation.
    pyrcel.legacy.thermo.kohler_crit : NumPy equivalent using ``scipy.optimize.fminbound``.
    """
    solver = optx.Bisection(rtol=rtol, atol=atol, expand_if_necessary=True)  # pyrefly: ignore[missing-argument]

    def dseq_dr(r, args):
        rd, kap = args
        return jax.grad(lambda rr: Seq(rr, rd, T, kap))(r)

    lower = r_dry
    upper = r_dry * 1e4
    guess = jnp.sqrt(lower * upper)
    sol = optx.root_find(
        dseq_dr,
        solver,
        guess,
        args=(r_dry, kappa),
        options=dict(lower=lower, upper=upper),
        max_steps=max_steps,
        throw=False,
    )
    return sol.value


def equilibrate_radii(
    T0: float,
    S0: float,
    r_drys: ArrayLike,
    kappas: ArrayLike,
    *,
    rtol: float = _RTOL,
    atol: float = _ATOL,
    max_steps: int = _MAX_STEPS,
) -> Array:
    """Equilibrium wet radii for every aerosol bin (vectorized over bins).

    Parameters
    ----------
    T0 : float
        Initial temperature, K.
    S0 : float
        Initial supersaturation (0 == 100% RH); sub-critical (typically < 0).
    r_drys, kappas : array, shape ``(nr,)``
        Dry radii (m) and hygroscopicities.

    Returns
    -------
    array, shape ``(nr,)``
        Equilibrium wet radii ``r0``, each in ``(r_dry, r_crit)``.
    """
    r_drys = jnp.asarray(r_drys)
    kappas = jnp.asarray(kappas)
    solver = optx.Bisection(rtol=rtol, atol=atol, expand_if_necessary=True)  # pyrefly: ignore[missing-argument]

    def residual(r, args):
        r_dry, kappa = args
        return Seq(r, r_dry, T0, kappa) - S0

    def solve_one(r_dry, kappa):
        r_crit = kohler_crit(T0, r_dry, kappa, rtol=rtol, atol=atol, max_steps=max_steps)
        lower = r_dry
        upper = r_crit
        guess = 0.5 * (lower + upper)
        sol = optx.root_find(
            residual,
            solver,
            guess,
            args=(r_dry, kappa),
            options=dict(lower=lower, upper=upper),
            max_steps=max_steps,
            throw=False,
        )
        return sol.value

    return jax.vmap(solve_one)(r_drys, kappas)


def equilibrate_initial_state(
    T0: float,
    S0: float,
    P0: float,
    r_drys: ArrayLike,
    kappas: ArrayLike,
    Nis: ArrayLike,
    **kwargs: Any,
) -> Array:
    """Assemble the equilibrated initial state vector ``y0``.

    Faithful to ``ParcelModel._setup_run``: water-vapor mixing ratio from the
    ambient RH, equilibrium wet radii, and the droplet liquid-water content they
    imply (normalized by the dry-air density).

    Parameters
    ----------
    T0 : float
        Initial temperature, K.
    S0 : float
        Initial supersaturation (0 == 100% RH); typically sub-critical (< 0).
    P0 : float
        Initial pressure, Pa.
    r_drys : array, shape ``(nr,)``
        Dry radii, m.
    kappas : array, shape ``(nr,)``
        Hygroscopicities.
    Nis : array, shape ``(nr,)``
        Number concentrations, m^-3.
    **kwargs
        Forwarded to [equilibrate_radii][pyrcel.equilibrate.equilibrate_radii] (``rtol``, ``atol``,
        ``max_steps``).

    Returns
    -------
    array, shape ``(7 + nr,)``
        ``[z=0, P0, T0, wv0, wc0, wi0=0, S0, r0_0 ... r0_{nr-1}]``.
    """
    r_drys = jnp.asarray(r_drys)
    kappas = jnp.asarray(kappas)
    Nis = jnp.asarray(Nis)

    es0 = es(T0 - 273.15)
    wv0 = (S0 + 1.0) * (c.epsilon * es0 / (P0 - es0))

    r0s = equilibrate_radii(T0, S0, r_drys, kappas, **kwargs)

    water_vol = (4.0 * jnp.pi / 3.0) * c.rho_w * Nis * (r0s**3 - r_drys**3)
    wc0 = jnp.sum(water_vol) / rho_air(T0, P0, 0.0)
    wi0 = 0.0

    bulk = jnp.array([0.0, P0, T0, wv0, wc0, wi0, S0])
    return jnp.concatenate([bulk, r0s])
