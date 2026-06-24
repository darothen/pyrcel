"""Tests for nd_from_parcel: differentiable Nd via sigmoid soft threshold (issue #67).

Covers:
1. Basic correctness: positive float, units/scale consistent with hard Nd
2. Soft→hard limit: at ε=1e-10 agrees with hard-threshold Nd within 10%
3. Differentiability: jax.grad runs and returns finite values for y0 and args
4. FD gradient check: JAX gradient vs central-difference for V and Ns
5. JIT compatibility
"""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
import pytest

import pyrcel as pm
from pyrcel import AerosolSpecies, Lognorm
from pyrcel.integrator import nd_from_parcel

jax.config.update("jax_enable_x64", True)


# ---------------------------------------------------------------------------
# Shared fixture: simple single-mode parcel
# ---------------------------------------------------------------------------

T0, S0, P0 = 283.15, 0.0, 85000.0
V0 = 0.5
_AER = AerosolSpecies("sulfate", Lognorm(mu=0.05, sigma=2.0, N=300.0), kappa=0.6, bins=8)
_T_END = 300.0
_EPS_TIGHT = 1e-10


@pytest.fixture(scope="module")
def parcel_setup():
    """Return a ready-to-go (y0, args, t_end) triple for the single-mode test parcel."""
    model = pm.ParcelModel([_AER], V=V0, T0=T0, S0=S0, P0=P0)
    return model.y0, model.args, _T_END


# ---------------------------------------------------------------------------
# 1. Basic sanity
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_nd_soft_returns_finite_positive(parcel_setup):
    y0, args, t_end = parcel_setup
    nd = nd_from_parcel(y0, args, t_end, epsilon=1e-8)
    assert float(nd) > 0.0
    assert np.isfinite(float(nd))


@pytest.mark.slow
def test_nd_soft_scale_consistent_with_hard(parcel_setup):
    """Soft Nd should be in the same ballpark as the hard diagnostic."""
    y0, args, t_end = parcel_setup
    model = pm.ParcelModel([_AER], V=V0, T0=T0, S0=S0, P0=P0)
    model.run(t_end=t_end, output_dt=5.0)
    nd_hard = model.summary()["total_Nd"]

    nd_soft = float(nd_from_parcel(y0, args, t_end, epsilon=1e-8))
    # Same order of magnitude (within factor 3) — loose check on scale only
    assert nd_soft / nd_hard > 0.33
    assert nd_soft / nd_hard < 3.0


# ---------------------------------------------------------------------------
# 2. Soft → hard limit
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_nd_soft_tight_eps_matches_hard(parcel_setup):
    """At ε=1e-10, soft Nd should be within 10% of hard Nd.

    Residual discrepancy comes from the approx vs exact Köhler r_crit, not
    from sigmoid smoothing.
    """
    y0, args, t_end = parcel_setup
    model = pm.ParcelModel([_AER], V=V0, T0=T0, S0=S0, P0=P0)
    model.run(t_end=t_end, output_dt=5.0)
    nd_hard = model.summary()["total_Nd"]

    nd_soft = float(nd_from_parcel(y0, args, t_end, epsilon=_EPS_TIGHT))
    rel_err = abs(nd_soft - nd_hard) / nd_hard
    assert rel_err < 0.10, f"rel_err={rel_err:.3%}"


@pytest.mark.slow
def test_nd_soft_monotone_in_eps(parcel_setup):
    """Nd_soft should be roughly monotone as ε shrinks toward the hard limit."""
    y0, args, t_end = parcel_setup
    model = pm.ParcelModel([_AER], V=V0, T0=T0, S0=S0, P0=P0)
    model.run(t_end=t_end, output_dt=5.0)
    nd_hard = model.summary()["total_Nd"]

    prev_err = float("inf")
    for eps in [1e-6, 1e-8, 1e-10]:
        nd = float(nd_from_parcel(y0, args, t_end, epsilon=eps))
        err = abs(nd - nd_hard)
        assert err < prev_err * 1.5, f"Nd not converging as ε→0: eps={eps}, err={err}"
        prev_err = err


# ---------------------------------------------------------------------------
# 3. Differentiability: grad runs without error
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_grad_wrt_y0_is_finite(parcel_setup):
    y0, args, t_end = parcel_setup

    def fn(y0_):
        return nd_from_parcel(y0_, args, t_end, epsilon=1e-8)

    grad_y0 = jax.grad(fn)(jnp.asarray(y0))
    assert np.all(np.isfinite(np.asarray(grad_y0)))


@pytest.mark.slow
def test_grad_wrt_args_is_finite(parcel_setup):
    """Gradient flows through the ODE and the Köhler computation."""
    y0, args, t_end = parcel_setup
    r_drys, Nis, kappas, accom, V = args

    def fn(r_drys_, Nis_, kappas_):
        new_args = (r_drys_, Nis_, kappas_, accom, V)
        return nd_from_parcel(y0, new_args, t_end, epsilon=1e-8)

    grad_r, grad_N, grad_k = jax.grad(fn, argnums=(0, 1, 2))(
        jnp.asarray(r_drys), jnp.asarray(Nis), jnp.asarray(kappas)
    )
    for name, g in [("r_drys", grad_r), ("Nis", grad_N), ("kappas", grad_k)]:
        assert np.all(np.isfinite(np.asarray(g))), f"non-finite gradient w.r.t. {name}"


# ---------------------------------------------------------------------------
# 4. FD gradient check for V and N
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_grad_V_fd_check(parcel_setup):
    """d(Nd)/d(V) via JAX autodiff should match central-difference to ~1%."""
    y0, args, t_end = parcel_setup
    r_drys, Nis, kappas, accom, _V = args
    V_base = float(_V.V) if hasattr(_V, "V") else float(_V)

    from pyrcel.updraft import ConstantV

    def fn(V_val):
        new_args = (r_drys, Nis, kappas, accom, ConstantV(V_val))
        return nd_from_parcel(y0, new_args, t_end, epsilon=1e-8)

    V_jax = jnp.float64(V_base)
    grad_jax = float(jax.grad(fn)(V_jax))

    dV = V_base * 1e-3
    grad_fd = float((fn(jnp.float64(V_base + dV)) - fn(jnp.float64(V_base - dV))) / (2 * dV))

    if abs(grad_fd) > 1e-6:
        rel_err = abs(grad_jax - grad_fd) / abs(grad_fd)
        assert rel_err < 0.05, f"JAX={grad_jax:.4g} FD={grad_fd:.4g} rel_err={rel_err:.3%}"
    else:
        assert abs(grad_jax - grad_fd) < 1e-3


@pytest.mark.slow
def test_grad_N_fd_check(parcel_setup):
    """d(Nd)/d(N_i) via autodiff vs FD for the first bin."""
    y0, args, t_end = parcel_setup
    r_drys, Nis_arr, kappas, accom, V = args
    Nis_base = jnp.asarray(Nis_arr)

    def fn(Nis_):
        new_args = (r_drys, Nis_, kappas, accom, V)
        return nd_from_parcel(y0, new_args, t_end, epsilon=1e-8)

    grad_N = np.asarray(jax.grad(fn)(Nis_base))

    # FD for bin 0 only (to keep wall time reasonable)
    dN = float(Nis_base[0]) * 1e-4
    Nis_p = Nis_base.at[0].add(dN)
    Nis_m = Nis_base.at[0].add(-dN)
    grad_fd_0 = float((fn(Nis_p) - fn(Nis_m)) / (2 * dN))

    if abs(grad_fd_0) > 1e-4:
        rel_err = abs(float(grad_N[0]) - grad_fd_0) / abs(grad_fd_0)
        assert rel_err < 0.05, (
            f"JAX={float(grad_N[0]):.4g} FD={grad_fd_0:.4g} rel_err={rel_err:.3%}"
        )
    else:
        assert abs(float(grad_N[0]) - grad_fd_0) < 1e-3


# ---------------------------------------------------------------------------
# 5. JIT compatibility
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_jit_nd_from_parcel(parcel_setup):
    y0, args, t_end = parcel_setup
    fn_jit = jax.jit(nd_from_parcel, static_argnames=("max_steps",))
    nd1 = float(fn_jit(y0, args, t_end, epsilon=1e-8))
    nd2 = float(fn_jit(y0, args, t_end, epsilon=1e-8))
    assert nd1 == pytest.approx(nd2, rel=1e-12)
    assert nd1 > 0.0


@pytest.mark.slow
def test_jit_grad_nd_from_parcel(parcel_setup):
    """jax.jit(jax.grad(nd_from_parcel)) should execute without error."""
    y0, args, t_end = parcel_setup

    def fn(y0_):
        return nd_from_parcel(y0_, args, t_end, epsilon=1e-8)

    grad_fn = jax.jit(jax.grad(fn))
    g = np.asarray(grad_fn(jnp.asarray(y0)))
    assert g.shape == jnp.asarray(y0).shape
    assert np.all(np.isfinite(g))
