"""Tests for pyrcel.activation JAX-native parameterizations.

Two layers:
1. **Cross-check** — JAX :func:`arg2000` vs legacy NumPy :func:`arg2000` over
   a grid of scenarios; agreement within float64 round-off.
2. **Gradient check** — :func:`jax.grad` of ``smax`` w.r.t. each input
   (V, T, P, mus, Ns) agrees with central finite differences to ~1e-4 relative.
3. **Physical invariants** — smax increases with V; activated fraction is in [0,1].
"""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
import pytest

jax.config.update("jax_enable_x64", True)

from pyrcel.activation import arg2000 as arg2000_jax  # noqa: E402
from pyrcel.legacy.activation import arg2000 as arg2000_legacy  # noqa: E402

# ---------------------------------------------------------------------------
# Test scenarios
# ---------------------------------------------------------------------------

SINGLE = dict(mus=[0.05], sigmas=[2.0], Ns=[300.0], kappas=[0.6])
TWO_MODE = dict(
    mus=[0.02, 0.1],
    sigmas=[1.5, 2.0],
    Ns=[500.0, 200.0],
    kappas=[0.6, 0.2],
)
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
# Cross-check: JAX vs legacy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("aerosol", [SINGLE, TWO_MODE, MARINE])
@pytest.mark.parametrize("V,T,P", UPDRAFT_CONDITIONS)
def test_arg2000_smax_vs_legacy(aerosol, V, T, P):
    """JAX and legacy arg2000 return the same smax to within 1e-6 relative."""
    smax_jax, _, _ = arg2000_jax(V, T, P, **aerosol)
    smax_leg, _, _ = arg2000_legacy(V, T, P, **aerosol)
    assert float(smax_jax) == pytest.approx(float(smax_leg), rel=1e-6)


@pytest.mark.parametrize("aerosol", [SINGLE, TWO_MODE])
@pytest.mark.parametrize("V,T,P", [(0.5, 283.15, 85000.0)])
def test_arg2000_act_fracs_vs_legacy(aerosol, V, T, P):
    """Activated fractions agree within 1e-5 relative."""
    _, N_acts_jax, fracs_jax = arg2000_jax(V, T, P, **aerosol)
    _, N_acts_leg, fracs_leg = arg2000_legacy(V, T, P, **aerosol)
    np.testing.assert_allclose(np.asarray(fracs_jax), np.asarray(fracs_leg), rtol=1e-5)
    np.testing.assert_allclose(np.asarray(N_acts_jax), np.asarray(N_acts_leg), rtol=1e-5)


def test_arg2000_accom_nonunity_vs_legacy():
    """Non-unity accommodation correction agrees with legacy."""
    kw = dict(**SINGLE, accom=0.042)
    V, T, P = 0.5, 283.15, 85000.0
    smax_jax, _, _ = arg2000_jax(V, T, P, **kw)
    smax_leg, _, _ = arg2000_legacy(V, T, P, **kw)
    assert float(smax_jax) == pytest.approx(float(smax_leg), rel=1e-6)


# ---------------------------------------------------------------------------
# Gradient check: jax.grad vs central finite differences
# ---------------------------------------------------------------------------

_BASE = dict(V=0.5, T=283.15, P=85000.0, mus=[0.05], sigmas=[2.0], Ns=[300.0], kappas=[0.6])


def _smax(**kw):
    smax, _, _ = arg2000_jax(**kw)
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
        ("V", 1e-4),
        ("T", 1e-3),
        ("P", 1.0),
    ],
)
def test_arg2000_scalar_gradient(param, delta):
    """jax.grad w.r.t. scalar inputs matches central FD to 1e-4 relative."""
    base = {**_BASE}

    def fn(x):
        kw = dict(base)
        kw[param] = x
        return _smax(**kw)

    grad_jax = float(jax.grad(fn)(float(base[param])))
    grad_fd = _fd_grad(param, delta, **base)
    assert grad_jax == pytest.approx(grad_fd, rel=1e-4)


def test_arg2000_gradient_wrt_N():
    """jax.grad of smax w.r.t. Ns[0] agrees with central FD."""
    Ns = jnp.array([300.0])
    sigmas = jnp.array([2.0])
    mus = jnp.array([0.05])
    kappas = jnp.array([0.6])
    V, T, P = 0.5, 283.15, 85000.0

    def fn(Ns):
        smax, _, _ = arg2000_jax(V, T, P, mus, sigmas, Ns, kappas)
        return smax

    grad_jax = float(jax.grad(fn)(Ns)[0])
    dN = 0.1
    grad_fd = (float(fn(Ns + dN)) - float(fn(Ns - dN))) / (2.0 * dN)
    assert grad_jax == pytest.approx(grad_fd, rel=1e-4)


def test_arg2000_gradient_wrt_kappa():
    """jax.grad of smax w.r.t. kappas[0] agrees with central FD."""
    kappas = jnp.array([0.6])
    mus = jnp.array([0.05])
    sigmas = jnp.array([2.0])
    Ns = jnp.array([300.0])
    V, T, P = 0.5, 283.15, 85000.0

    def fn(kappas):
        smax, _, _ = arg2000_jax(V, T, P, mus, sigmas, Ns, kappas)
        return smax

    grad_jax = float(jax.grad(fn)(kappas)[0])
    dk = 1e-4
    grad_fd = (float(fn(kappas + dk)) - float(fn(kappas - dk))) / (2.0 * dk)
    assert grad_jax == pytest.approx(grad_fd, rel=1e-3)


# ---------------------------------------------------------------------------
# Physical invariants
# ---------------------------------------------------------------------------


def test_smax_increases_with_updraft():
    """Higher updraft → higher smax."""
    Vs = [0.1, 0.3, 0.5, 1.0, 2.0]
    smaxes = [float(arg2000_jax(V, 283.15, 85000.0, **SINGLE)[0]) for V in Vs]
    assert all(s1 < s2 for s1, s2 in zip(smaxes, smaxes[1:]))


def test_act_fracs_in_unit_interval():
    """Activated fractions must lie in [0, 1] for all test scenarios."""
    for aerosol in [SINGLE, TWO_MODE, MARINE]:
        for V, T, P in UPDRAFT_CONDITIONS:
            _, _, fracs = arg2000_jax(V, T, P, **aerosol)
            assert jnp.all(fracs >= 0.0) and jnp.all(fracs <= 1.0)


def test_smax_positive():
    """smax must be strictly positive."""
    smax, _, _ = arg2000_jax(0.5, 283.15, 85000.0, **SINGLE)
    assert float(smax) > 0.0


def test_two_mode_smax_lt_either_single_mode():
    """Multi-mode competition: joint smax < either single-mode smax (more CCN compete)."""
    V, T, P = 0.5, 283.15, 85000.0
    smax_all, _, _ = arg2000_jax(V, T, P, **TWO_MODE)
    mode1 = dict(
        mus=[TWO_MODE["mus"][0]],
        sigmas=[TWO_MODE["sigmas"][0]],
        Ns=[TWO_MODE["Ns"][0]],
        kappas=[TWO_MODE["kappas"][0]],
    )
    mode2 = dict(
        mus=[TWO_MODE["mus"][1]],
        sigmas=[TWO_MODE["sigmas"][1]],
        Ns=[TWO_MODE["Ns"][1]],
        kappas=[TWO_MODE["kappas"][1]],
    )
    smax1, _, _ = arg2000_jax(V, T, P, **mode1)
    smax2, _, _ = arg2000_jax(V, T, P, **mode2)
    assert float(smax_all) < float(min(smax1, smax2))


# ---------------------------------------------------------------------------
# ARG2000 class interface
# ---------------------------------------------------------------------------


def test_arg2000_class_matches_function():
    """ARG2000() class gives identical results to the arg2000() function."""
    from pyrcel.activation import ARG2000

    scheme = ARG2000()
    V, T, P = 0.5, 283.15, 85000.0
    smax_cls, N_cls, f_cls = scheme(V, T, P, **SINGLE)
    smax_fn, N_fn, f_fn = arg2000_jax(V, T, P, **SINGLE)
    assert float(smax_cls) == pytest.approx(float(smax_fn), rel=1e-12)


def test_arg2000_class_repr():
    from pyrcel.activation import ARG2000

    r = repr(ARG2000())
    assert "ARG2000" in r
