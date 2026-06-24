"""JAX-native Morales Betancourt & Nenes (2014) activation parameterization.

The scheme finds the maximum parcel supersaturation by solving the implicit
equation F(smax, θ) = 0 via bisection, then computes per-mode activated
fractions from lognormal statistics.  The bisection uses
:func:`jax.lax.while_loop` so it is JIT-compilable but not differentiable
through the loop itself.  Differentiability w.r.t. the physical inputs is
recovered by the **implicit function theorem** (IFT), implemented via
:func:`jax.custom_vjp`:

    ∂smax/∂θ  =  −(∂F/∂θ) / (∂F/∂smax)   at  F(smax, θ) = 0

Both partial derivatives are evaluated by applying :func:`jax.grad` to the
residual :func:`_mbn_residual`, which is a plain differentiable JAX function
(no while_loop).

References
----------
.. [MBN2014] Morales Betancourt, R. and Nenes, A.: Droplet activation
   parameterization: the population splitting concept revisited,
   Geosci. Model Dev. Discuss., 7, 2903–2932,
   doi:10.5194/gmdd-7-2903-2014, 2014.
"""

from __future__ import annotations

import jax

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402
from jax import Array  # noqa: E402
from jax.typing import ArrayLike  # noqa: E402

from .. import constants as c  # noqa: E402
from ..thermo import ka_cont  # noqa: E402
from ._arg2000 import _lognormal_act  # noqa: E402

# Bisection search bounds and convergence settings (match legacy defaults).
_XMIN: float = 1e-5
_XMAX: float = 0.5
_TOL: float = 1e-6
_MAX_ITERS: int = 200


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _vpres_jax(T: ArrayLike) -> Array:
    """Saturated water vapour pressure (mb) — polynomial fit.

    Uses the same coefficients and reference temperature (273 K, not 273.15 K)
    as :func:`pyrcel.legacy.activation._vpres` to preserve numerical agreement
    with the legacy implementation.
    """
    A = jnp.array(
        [
            6.107799610e0,
            4.436518521e-1,
            1.428945805e-2,
            2.650648471e-4,
            3.031240396e-6,
            2.034080948e-8,
            6.136820929e-11,
        ]
    )
    Tc = T - 273.0  # legacy uses 273, not 273.15
    vp = A[6]
    for ai in [A[5], A[4], A[3], A[2], A[1]]:
        vp = vp * Tc + ai
    return vp * Tc + A[0]


def _erfp_jax(x: ArrayLike) -> Array:
    """Polynomial erf approximation — JAX-traceable analogue of legacy _erfp.

    Used *only* inside the MBN2014 condensation integrals to match the
    published scheme exactly.  The final lognormal activation step uses
    :func:`jax.scipy.special.erfc` (accurate).
    """
    AA = jnp.array([0.278393, 0.230389, 0.000972, 0.078108])
    y = jnp.abs(x)
    axx = 1.0 + y * (AA[0] + y * (AA[1] + y * (AA[2] + y * AA[3])))
    axx = axx**4
    axx = 1.0 - 1.0 / axx
    return jnp.where(x <= 0.0, -axx, axx)


def _sintegral_jax(
    smax: Array,
    sgis: Array,
    Ns_m: Array,
    sigmas: Array,
    A: Array,
    beta2: Array,
    alpha: Array,
    V: Array,
) -> tuple[Array, Array]:
    """Vectorised condensation integrals I₁ and I₂ from MBN2014 (App. B).

    Parameters
    ----------
    smax : float
        Trial maximum supersaturation (decimal).
    sgis : array, shape (n,)
        Modal critical supersaturations (decimal).
    Ns_m : array, shape (n,)
        Modal number concentrations (m⁻³).
    sigmas : array, shape (n,)
        Geometric standard deviations.
    A : float
        Köhler surface-tension coefficient (m).
    beta2 : float
        Inverse condensational growth coefficient (s/m²).
    alpha : float
        Thermodynamic alpha coefficient (1/m).
    V : float
        Updraft speed (m/s).

    Returns
    -------
    I1, I2 : float
        Summed condensation integrals (SI units).
    """
    sqrtwo = jnp.sqrt(2.0)
    log_sigmas = jnp.log(sigmas)
    d_eq = A * 2.0 / sgis / 3.0 / jnp.sqrt(3.0)

    # --- population-splitting supersaturation ssplt2 (and ssplt1) -----------
    zeta_c = ((16.0 / 9.0) * alpha * V * beta2 * A**2) ** 0.25
    delta = 1.0 - (zeta_c / smax) ** 4.0
    critical = delta <= 0.0

    # Critical branch: ratio formula from MBN2014 eq. (B12)
    ratio_c = (1.0 / sqrtwo) + (2e7 / 3.0) * A * (smax**-0.3824 - zeta_c**-0.3824)
    ssplt2_c = smax * jnp.minimum(ratio_c, 1.0)

    # Non-critical branch: roots of the quartic, eq. (B10)
    safe_delta = jnp.maximum(delta, 0.0)
    sqrt_delta = jnp.sqrt(safe_delta)
    ssplt2_nc = jnp.sqrt(0.5 * (1.0 + sqrt_delta)) * smax
    # ssplt1 could be 0 when delta → 0; guard the log below
    safe_ssplt1_sq = jnp.maximum(0.5 * (1.0 - sqrt_delta), 0.0)
    ssplt1_nc = jnp.sqrt(safe_ssplt1_sq) * smax

    ssplt2 = jnp.where(critical, ssplt2_c, ssplt2_nc)
    ssplt1 = ssplt1_nc  # guarded; only selected in non-critical branch

    # --- integ2: same formula in both branches (eq. B5) ---------------------
    log_sgi_sp2 = jnp.log(sgis / ssplt2)
    log_sgi_smax = jnp.log(sgis / smax)
    log_factor = 3.0 * log_sigmas / (2.0 * sqrtwo)
    u_sp2 = 2.0 * log_sgi_sp2 / (3.0 * sqrtwo * log_sigmas)
    u_smax = 2.0 * log_sgi_smax / (3.0 * sqrtwo * log_sigmas)
    exp98 = jnp.exp(9.0 / 8.0 * log_sigmas**2)
    integ2 = (exp98 * Ns_m / sgis) * (
        _erfp_jax(u_sp2 - log_factor) - _erfp_jax(u_smax - log_factor)
    )

    # --- integ1 + I_extra: differs between branches --------------------------
    # Both branches are evaluated (jnp.where); guarded quantities avoid NaN.

    # Critical branch (eq. B13)
    u_sp_plus_c = sqrtwo * log_sgi_sp2 / 3.0 / log_sigmas
    I_extra_c = (
        Ns_m
        * d_eq
        * exp98
        * (1.0 - _erfp_jax(u_sp_plus_c - log_factor))
        * jnp.sqrt(beta2 * alpha * V)
    )
    integ1_c = jnp.zeros_like(Ns_m)

    # Non-critical branch (eqs. B3–B4, B6)
    g_i = jnp.exp(4.5 * log_sigmas**2)
    sg_ratio_sq = (sgis / smax) ** 2
    safe_ssplt1 = jnp.maximum(ssplt1, 1e-30 * smax)
    log_sgi_sp1 = jnp.log(sgis / safe_ssplt1)
    u_sp1 = 2.0 * log_sgi_sp1 / (3.0 * sqrtwo * log_sigmas)
    shift = 3.0 * log_sigmas / sqrtwo

    def _i1_part(u_sp: Array) -> Array:
        return (
            Ns_m
            * smax
            * ((1.0 - _erfp_jax(u_sp)) - 0.5 * sg_ratio_sq * g_i * (1.0 - _erfp_jax(u_sp + shift)))
        )

    integ1_nc = _i1_part(u_sp2) - _i1_part(u_sp1)
    u_sp1_inertial = sqrtwo * log_sgi_sp1 / 3.0 / log_sigmas
    I_extra_nc = (
        Ns_m
        * d_eq
        * exp98
        * (1.0 - _erfp_jax(u_sp1_inertial - log_factor))
        * jnp.sqrt(beta2 * alpha * V)
    )

    integ1 = jnp.where(critical, integ1_c, integ1_nc)
    I_extra = jnp.where(critical, I_extra_c, I_extra_nc)

    return jnp.sum(integ1 + I_extra), jnp.sum(integ2)


def _mbn_residual(
    smax: Array,
    V: Array,
    T: Array,
    P: Array,
    mus_m: Array,
    sigmas: Array,
    Ns_m: Array,
    kappas: Array,
    accom: Array,
) -> Array:
    """Implicit equation F(smax, θ) = 0 whose root defines smax.

    All operations are standard JAX, making this fully differentiable w.r.t.
    every argument (used by the IFT backward pass).

    Parameters
    ----------
    smax : float
        Candidate supersaturation.
    V, T, P : float
        Updraft speed (m/s), temperature (K), pressure (Pa).
    mus_m : array, shape (n,)
        Median dry *radii* in **metres** (not μm; conversion done by caller).
    sigmas : array, shape (n,)
        Geometric standard deviations.
    Ns_m : array, shape (n,)
        Number concentrations in **m⁻³**.
    kappas : array, shape (n,)
        Hygroscopicity parameters.
    accom : float
        Condensation accommodation coefficient.

    Returns
    -------
    float
        Residual F(smax, θ).  Zero at the physical solution.
    """
    rho_a = P * c.Ma / c.R / T

    # Surface tension (273 K reference to match legacy implementation exactly)
    surt = 0.0761 - 1.55e-4 * (T - 273.0)
    A = 4.0 * c.Mw * surt / c.R / T / c.rho_w

    # Critical supersaturation per mode (Köhler approximation, eq. 5 of MBN2014)
    dpgs = 2.0 * mus_m  # diameter in m
    sgis = jnp.exp(jnp.sqrt(4.0 * A**3 / 27.0 / kappas / dpgs**3)) - 1.0

    # Saturated vapour pressure (mb → Pa)
    wv_pres_sat = _vpres_jax(T) * (1e5 / 1e3)

    # Thermodynamic alpha / beta1 (same as parcel ODE alpha / gamma_p)
    alpha = c.g * c.Mw * c.L / c.Cp / c.R / T**2 - c.g * c.Ma / c.R / T
    beta1 = P * c.Ma / wv_pres_sat / c.Mw + c.Mw * c.L**2 / c.Cp / c.R / T**2

    # Nenes diffusivity correction (continuum form matched to legacy)
    dv_orig = 1e-4 * (0.211 / (P / 1.013e5) * (T / 273.0) ** 1.94)
    dv_big = 5.0e-6
    dv_low = 1e-6 * 0.207683 * accom**-0.33048
    coef = jnp.sqrt(2.0 * jnp.pi * c.Mw / (c.R * T))
    dv_lam = 2.0 * dv_orig / accom * coef  # mean-free-path correction length
    dv_ave = (dv_orig / (dv_big - dv_low)) * (
        (dv_big - dv_low) - dv_lam * jnp.log((dv_big + dv_lam) / (dv_low + dv_lam))
    )

    # beta2 = 1/G (inverse condensational growth coefficient)
    beta2 = c.R * T * c.rho_w / wv_pres_sat / dv_ave / c.Mw / 4.0 + (
        c.L * c.rho_w / 4.0 / ka_cont(T) / T * (c.L * c.Mw / c.R / T - 1.0)
    )

    beta = 0.5 * jnp.pi * beta1 * c.rho_w / beta2 / alpha / V / rho_a
    cf1 = 0.5 * jnp.sqrt(1.0 / beta2 / (alpha * V))
    cf2 = A / 3.0

    I1, I2 = _sintegral_jax(smax, sgis, Ns_m, sigmas, A, beta2, alpha, V)
    return (I1 * cf1 + I2 * cf2) * beta * smax - 1.0


# ---------------------------------------------------------------------------
# while_loop bisection — JIT-compilable but not differentiable through loop
# ---------------------------------------------------------------------------


def _bisect(F, x_lo: Array, x_hi: Array) -> Array:
    """Bisect F on [x_lo, x_hi] using jax.lax.while_loop."""
    y_lo = F(x_lo)
    y_hi = F(x_hi)  # noqa: F841  # evaluated to establish bracket validity

    def cond(carry):
        x1, _, x2, _, it = carry
        return (jnp.abs(x2 - x1) > _TOL * jnp.abs(x1)) & (it < _MAX_ITERS)

    def body(carry):
        x1, y1, x2, y2, it = carry
        xm = 0.5 * (x1 + x2)
        ym = F(xm)
        same_sign = (y1 * ym) > 0.0
        return (
            jnp.where(same_sign, xm, x1),
            jnp.where(same_sign, ym, y1),
            jnp.where(same_sign, x2, xm),
            jnp.where(same_sign, y2, ym),
            it + 1,
        )

    x1, _, x2, _, _ = jax.lax.while_loop(cond, body, (x_lo, y_lo, x_hi, y_hi, 0))
    return 0.5 * (x1 + x2)


# ---------------------------------------------------------------------------
# custom_vjp wrapper — IFT gradient over the bisection
# ---------------------------------------------------------------------------


@jax.custom_vjp
def _smax_bisect(
    V: Array,
    T: Array,
    P: Array,
    mus_m: Array,
    sigmas: Array,
    Ns_m: Array,
    kappas: Array,
    accom: Array,
) -> Array:
    """Find smax via while_loop bisection (forward only; IFT provides the VJP)."""

    def F(s: Array) -> Array:
        return _mbn_residual(s, V, T, P, mus_m, sigmas, Ns_m, kappas, accom)

    return _bisect(F, jnp.asarray(_XMIN), jnp.asarray(_XMAX))


def _smax_bisect_fwd(
    V: Array,
    T: Array,
    P: Array,
    mus_m: Array,
    sigmas: Array,
    Ns_m: Array,
    kappas: Array,
    accom: Array,
) -> tuple[Array, tuple]:
    smax = _smax_bisect(V, T, P, mus_m, sigmas, Ns_m, kappas, accom)
    return smax, (smax, V, T, P, mus_m, sigmas, Ns_m, kappas, accom)


def _smax_bisect_bwd(res: tuple, g: Array) -> tuple:
    smax, V, T, P, mus_m, sigmas, Ns_m, kappas, accom = res
    # Evaluate ∂F/∂smax and ∂F/∂θ at the solution via automatic differentiation
    # of the differentiable residual (no while_loop involved here).
    all_grads = jax.grad(_mbn_residual, argnums=range(9))(
        smax, V, T, P, mus_m, sigmas, Ns_m, kappas, accom
    )
    dF_dsmax = all_grads[0]
    # IFT: ∂smax/∂θᵢ = −(∂F/∂θᵢ) / (∂F/∂smax)
    scale = -g / dF_dsmax
    return tuple(
        scale * gi for gi in all_grads[1:]
    )  # (dV, dT, dP, dmus, dsigmas, dNs, dkappas, daccom)


_smax_bisect.defvjp(_smax_bisect_fwd, _smax_bisect_bwd)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def mbn2014(
    V: ArrayLike,
    T: ArrayLike,
    P: ArrayLike,
    mus: ArrayLike,
    sigmas: ArrayLike,
    Ns: ArrayLike,
    kappas: ArrayLike,
    accom: float = c.ac,
) -> tuple[Array, Array, Array]:
    """JAX-native Morales Betancourt & Nenes (2014) activation parameterization.

    A faithful JAX re-implementation of :func:`pyrcel.legacy.activation.mbn2014`.
    The maximum parcel supersaturation is found by bisecting the implicit
    equation derived in [MBN2014]_ using :func:`jax.lax.while_loop`.
    Gradients w.r.t. all physical inputs are available via :func:`jax.grad`
    through the **implicit function theorem** (IFT).

    Parameters
    ----------
    V : float
        Updraft speed (m/s).
    T : float
        Parcel temperature (K).
    P : float
        Parcel pressure (Pa).
    mus : array-like, shape (n_modes,)
        Median dry radii of each lognormal mode (μm).
    sigmas : array-like, shape (n_modes,)
        Geometric standard deviations (> 1).
    Ns : array-like, shape (n_modes,)
        Total number concentrations (cm⁻³).
    kappas : array-like, shape (n_modes,)
        Hygroscopicity parameters.
    accom : float, optional
        Condensation accommodation coefficient (default :data:`pyrcel.constants.ac`).

    Returns
    -------
    smax : float
        Maximum parcel supersaturation (decimal).
    N_acts : jnp.ndarray, shape (n_modes,)
        Activated number concentration per mode (cm⁻³).
    act_fracs : jnp.ndarray, shape (n_modes,)
        Activated number fraction per mode.

    References
    ----------
    .. [MBN2014] Morales Betancourt, R. and Nenes, A.: Droplet activation
       parameterization: the population splitting concept revisited,
       Geosci. Model Dev. Discuss., 7, 2903–2932,
       doi:10.5194/gmdd-7-2903-2014, 2014.

    See Also
    --------
    pyrcel.legacy.activation.mbn2014 : NumPy reference implementation.
    """
    mus = jnp.asarray(mus, dtype=float)
    sigmas = jnp.asarray(sigmas, dtype=float)
    Ns = jnp.asarray(Ns, dtype=float)
    kappas = jnp.asarray(kappas, dtype=float)
    V = jnp.asarray(V, dtype=float)
    T = jnp.asarray(T, dtype=float)
    P = jnp.asarray(P, dtype=float)
    accom = jnp.asarray(accom, dtype=float)

    mus_m = mus * 1e-6  # μm → m
    Ns_m = Ns * 1e6  # cm⁻³ → m⁻³

    smax = _smax_bisect(V, T, P, mus_m, sigmas, Ns_m, kappas, accom)

    # Recompute sgis at converged smax for the lognormal activation step
    surt = 0.0761 - 1.55e-4 * (T - 273.0)
    A = 4.0 * c.Mw * surt / c.R / T / c.rho_w
    dpgs = 2.0 * mus_m
    sgis = jnp.exp(jnp.sqrt(4.0 * A**3 / 27.0 / kappas / dpgs**3)) - 1.0

    N_acts_list = []
    act_fracs_list = []
    for i in range(len(mus)):
        N_act, act_frac = _lognormal_act(smax, sigmas[i], Ns[i], sgis[i])
        N_acts_list.append(N_act)
        act_fracs_list.append(act_frac)

    return smax, jnp.stack(N_acts_list), jnp.stack(act_fracs_list)


class MBN2014:
    """Morales Betancourt & Nenes (2014) activation scheme.

    A thin callable wrapper around :func:`mbn2014` satisfying the
    :class:`~pyrcel.activation.ActivationScheme` interface.

    Parameters
    ----------
    accom : float, optional
        Condensation accommodation coefficient forwarded to every call.
    """

    def __init__(self, accom: float = c.ac) -> None:
        self.accom = accom

    def __call__(
        self,
        V: ArrayLike,
        T: ArrayLike,
        P: ArrayLike,
        mus: ArrayLike,
        sigmas: ArrayLike,
        Ns: ArrayLike,
        kappas: ArrayLike,
        accom: float | None = None,
    ) -> tuple[Array, Array, Array]:
        return mbn2014(
            V, T, P, mus, sigmas, Ns, kappas, accom=self.accom if accom is None else accom
        )

    def __repr__(self) -> str:
        return f"MBN2014(accom={self.accom})"
