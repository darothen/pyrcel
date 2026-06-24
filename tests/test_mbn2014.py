"""Tests for the JAX-native MBN2014 activation parameterization (issue #56).

Three layers:
1. **Cross-check** — JAX :func:`mbn2014` vs legacy NumPy :func:`mbn2014` over
   a grid of scenarios; smax agreement within 1e-5 relative (bisection uses
   the same tolerance on both sides but floating-point order differs slightly).
2. **Gradient check** — :func:`jax.grad` of ``smax`` w.r.t. V, T, P, Ns, kappas
   agrees with central finite differences to ~1e-2 relative (IFT gradient at
   a converged root; FD is limited by bisection tolerance, not precision).
3. **Physical invariants** — smax increases with V; activated fraction is in [0, 1].
4. **Class interface** — :class:`MBN2014` delegates correctly.
"""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
import pytest

jax.config.update("jax_enable_x64", True)

from pyrcel.activation._mbn2014 import mbn2014 as mbn2014_jax  # noqa: E402
from pyrcel.legacy.activation import mbn2014 as mbn2014_legacy  # noqa: E402

# ---------------------------------------------------------------------------
# Shared scenarios
# ---------------------------------------------------------------------------

SINGLE = dict(mus=[0.05], sigmas=[2.0], Ns=[300.0], kappas=[0.6])
TWO_MODE = dict(mus=[0.02, 0.1], sigmas=[1.5, 2.0], Ns=[500.0, 200.0], kappas=[0.6, 0.2])
MARINE = dict(
    mus=[0.07 / 2.0, 0.62 / 2.0],
    sigmas=[2.0, 2.7],
    Ns=[60.0, 3.1],
    kappas=[1.12, 1.12],
)

UPDRAFT_CONDITIONS = [
    (0.25, 283.15, 90000.0),
    (0.50, 283.15, 85000.0),
    (1.00, 280.00, 80000.0),
    (2.00, 275.00, 75000.0),
]


# ---------------------------------------------------------------------------
# 1. Cross-check: JAX vs legacy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("aerosol", [SINGLE, TWO_MODE, MARINE])
@pytest.mark.parametrize("V,T,P", UPDRAFT_CONDITIONS)
def test_mbn2014_smax_vs_legacy(aerosol, V, T, P):
    """JAX and legacy mbn2014 return the same smax to within 1e-5 relative."""
    smax_jax, _, _ = mbn2014_jax(V, T, P, **aerosol)
    smax_leg, _, _ = mbn2014_legacy(V, T, P, **aerosol)
    assert float(smax_jax) == pytest.approx(float(smax_leg), rel=1e-5)


@pytest.mark.parametrize("aerosol", [SINGLE, TWO_MODE])
@pytest.mark.parametrize("V,T,P", [(0.5, 283.15, 85000.0)])
def test_mbn2014_act_fracs_vs_legacy(aerosol, V, T, P):
    """Activated fractions agree within 1e-4 relative."""
    _, _, fracs_jax = mbn2014_jax(V, T, P, **aerosol)
    _, _, fracs_leg = mbn2014_legacy(V, T, P, **aerosol)
    np.testing.assert_allclose(np.asarray(fracs_jax), np.asarray(fracs_leg), rtol=1e-4)


def test_mbn2014_accom_nonunity_vs_legacy():
    """Non-unity accommodation correction agrees with legacy."""
    V, T, P = 0.5, 283.15, 85000.0
    smax_jax, _, _ = mbn2014_jax(V, T, P, accom=0.042, **SINGLE)
    smax_leg, _, _ = mbn2014_legacy(V, T, P, accom=0.042, **SINGLE)
    assert float(smax_jax) == pytest.approx(float(smax_leg), rel=1e-5)


# ---------------------------------------------------------------------------
# 2. Gradient check: IFT jax.grad vs central finite differences
# ---------------------------------------------------------------------------

_BASE = dict(V=0.5, T=283.15, P=85000.0, mus=[0.05], sigmas=[2.0], Ns=[300.0], kappas=[0.6])


def _smax(**kw):
    smax, _, _ = mbn2014_jax(**kw)
    return smax


def _fd_grad(param, delta, **base):
    kw_p = dict(base)
    kw_m = dict(base)
    kw_p[param] = base[param] + delta
    kw_m[param] = base[param] - delta
    return (float(_smax(**kw_p)) - float(_smax(**kw_m))) / (2.0 * delta)


@pytest.mark.parametrize(
    "param, delta",
    [
        ("V", 1e-3),
        ("T", 1e-2),
        ("P", 10.0),
    ],
)
def test_mbn2014_scalar_gradient(param, delta):
    """IFT jax.grad w.r.t. scalar inputs matches central FD to 1e-2 relative.

    The looser tolerance (vs ARG2000's 1e-4) reflects the bisection tolerance
    floor: FD perturbs smax by ~tol, so the FD estimate itself carries ~1e-3
    relative error relative to the true gradient.
    """
    base = dict(_BASE)

    def fn(x):
        kw = dict(base)
        kw[param] = x
        return _smax(**kw)

    grad_jax = float(jax.grad(fn)(jnp.float64(base[param])))
    grad_fd = _fd_grad(param, delta, **base)
    assert grad_jax == pytest.approx(grad_fd, rel=1e-2)


def test_mbn2014_gradient_wrt_N():
    """IFT jax.grad of smax w.r.t. Ns[0] agrees with central FD."""
    Ns = jnp.array([300.0])
    V, T, P = 0.5, 283.15, 85000.0
    mus = jnp.array([0.05])
    sigmas = jnp.array([2.0])
    kappas = jnp.array([0.6])

    def fn(Ns):
        smax, _, _ = mbn2014_jax(V, T, P, mus, sigmas, Ns, kappas)
        return smax

    grad_jax = float(jax.grad(fn)(Ns)[0])
    dN = 1.0
    grad_fd = (float(fn(Ns + dN)) - float(fn(Ns - dN))) / (2.0 * dN)
    assert grad_jax == pytest.approx(grad_fd, rel=1e-2)


def test_mbn2014_gradient_wrt_kappa():
    """IFT jax.grad of smax w.r.t. kappas[0] agrees with central FD."""
    kappas = jnp.array([0.6])
    V, T, P = 0.5, 283.15, 85000.0
    mus = jnp.array([0.05])
    sigmas = jnp.array([2.0])
    Ns = jnp.array([300.0])

    def fn(kappas):
        smax, _, _ = mbn2014_jax(V, T, P, mus, sigmas, Ns, kappas)
        return smax

    grad_jax = float(jax.grad(fn)(kappas)[0])
    dk = 1e-3
    grad_fd = (float(fn(kappas + dk)) - float(fn(kappas - dk))) / (2.0 * dk)
    assert grad_jax == pytest.approx(grad_fd, rel=1e-2)


# ---------------------------------------------------------------------------
# 3. Physical invariants
# ---------------------------------------------------------------------------


def test_mbn2014_smax_increases_with_updraft():
    """Higher updraft → higher smax."""
    Vs = [0.1, 0.3, 0.5, 1.0, 2.0]
    smaxes = [float(mbn2014_jax(V, 283.15, 85000.0, **SINGLE)[0]) for V in Vs]
    assert all(s1 < s2 for s1, s2 in zip(smaxes, smaxes[1:]))


def test_mbn2014_act_fracs_in_unit_interval():
    """Activated fractions must lie in [0, 1] for all test scenarios."""
    for aerosol in [SINGLE, TWO_MODE, MARINE]:
        for V, T, P in UPDRAFT_CONDITIONS:
            _, _, fracs = mbn2014_jax(V, T, P, **aerosol)
            assert jnp.all(fracs >= 0.0) and jnp.all(fracs <= 1.0)


def test_mbn2014_smax_positive():
    """smax must be strictly positive."""
    smax, _, _ = mbn2014_jax(0.5, 283.15, 85000.0, **SINGLE)
    assert float(smax) > 0.0


# ---------------------------------------------------------------------------
# 4. MBN2014 class interface
# ---------------------------------------------------------------------------


def test_mbn2014_class_matches_function():
    """MBN2014() class gives identical results to the mbn2014() function."""
    from pyrcel.activation import MBN2014

    scheme = MBN2014()
    V, T, P = 0.5, 283.15, 85000.0
    smax_cls, N_cls, f_cls = scheme(V, T, P, **SINGLE)
    smax_fn, N_fn, f_fn = mbn2014_jax(V, T, P, **SINGLE)
    assert float(smax_cls) == pytest.approx(float(smax_fn), rel=1e-12)


def test_mbn2014_class_call_time_accom_overrides_instance():
    """Explicit accom at call time overrides the instance default."""
    from pyrcel.activation import MBN2014

    scheme = MBN2014(accom=1.0)
    V, T, P = 0.5, 283.15, 85000.0
    smax_default, _, _ = scheme(V, T, P, **SINGLE)
    smax_override, _, _ = scheme(V, T, P, **SINGLE, accom=0.042)
    assert float(smax_override) != pytest.approx(float(smax_default), rel=1e-4)


def test_mbn2014_class_repr():
    from pyrcel.activation import MBN2014

    r = repr(MBN2014())
    assert "MBN2014" in r


def test_mbn2014_public_import():
    """mbn2014 and MBN2014 are accessible from pyrcel top-level."""
    import pyrcel

    assert callable(pyrcel.MBN2014)
    from pyrcel.activation import MBN2014, mbn2014  # noqa: F401
