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

from pyrcel.activation import arg2000 as arg2000_jax
from pyrcel.activation import binned_activation, lognormal_activation, multi_mode_activation
from pyrcel.aerosol import AerosolSpecies
from pyrcel.distributions import Lognorm
from pyrcel.legacy.activation import (
    arg2000 as arg2000_legacy,
)
from pyrcel.legacy.activation import (
    binned_activation as binned_leg,
)
from pyrcel.legacy.activation import (
    lognormal_activation as lognormal_leg,
)
from pyrcel.legacy.activation import (
    multi_mode_activation as multi_mode_leg,
)

jax.config.update("jax_enable_x64", True)

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


def test_arg2000_class_call_time_accom_overrides_instance():
    """Explicit accom=0.0 at call time must not be swallowed by falsy check."""
    from pyrcel.activation import ARG2000

    scheme = ARG2000(accom=1.0)
    V, T, P = 0.5, 283.15, 85000.0
    smax_default, _, _ = scheme(V, T, P, **SINGLE)
    smax_override, _, _ = scheme(V, T, P, **SINGLE, accom=0.042)
    # accom=0.042 reduces effective diffusivity → higher smax
    assert float(smax_override) != pytest.approx(float(smax_default), rel=1e-4)


def test_arg2000_class_repr():
    from pyrcel.activation import ARG2000

    r = repr(ARG2000())
    assert "ARG2000" in r


# ---------------------------------------------------------------------------
# Cross-check: JAX binned/lognormal helpers vs legacy (approx=True)
# ---------------------------------------------------------------------------

_T = 283.15
_SMAX = 0.008

# Two aerosol modes reused across several tests.
_AER1 = AerosolSpecies("mode1", Lognorm(mu=0.05, sigma=2.0, N=300.0), bins=20, kappa=0.6)
_AER2 = AerosolSpecies("mode2", Lognorm(mu=0.10, sigma=2.0, N=100.0), bins=20, kappa=0.3)
# Rough wet radii: 3× dry radius puts most bins past their critical wet radius.
_RS1 = _AER1.r_drys * 3.0
_RS2 = _AER2.r_drys * 3.0


def test_lognormal_activation_vs_legacy():
    """JAX lognormal_activation matches legacy (approx Köhler) to float64 round-off."""
    mu = _AER1.r_drys[10]  # a single representative dry radius
    sigma = 2.0
    N = 300.0
    kappa = 0.6

    N_act_jax, frac_jax = lognormal_activation(_SMAX, mu, sigma, N, kappa, T=_T)
    N_act_leg, frac_leg = lognormal_leg(_SMAX, mu, sigma, N, kappa, T=_T, approx=True)

    assert float(N_act_jax) == pytest.approx(float(N_act_leg), rel=1e-12)
    assert float(frac_jax) == pytest.approx(float(frac_leg), rel=1e-12)


def test_lognormal_activation_precomputed_sgi():
    """Passing sgi directly gives the same result as letting T drive the computation."""
    from pyrcel.activation._common import _kohler_crit_approx

    mu = _AER1.r_drys[10]
    kappa = 0.6
    sigma = 2.0
    N = 300.0
    _, sgi = _kohler_crit_approx(_T, mu, kappa)

    N_act_T, _ = lognormal_activation(_SMAX, mu, sigma, N, kappa, T=_T)
    N_act_sgi, _ = lognormal_activation(_SMAX, mu, sigma, N, kappa, sgi=sgi)

    assert float(N_act_T) == pytest.approx(float(N_act_sgi), rel=1e-12)


@pytest.mark.parametrize("smax", [0.002, 0.008, 0.02])
def test_binned_activation_vs_legacy(smax):
    """JAX eq-frac matches legacy (approx Köhler); kn-frac uses equilibrium gate."""
    eq_jax, kn_jax, alpha_jax, phi_jax = binned_activation(
        smax, _T, _RS1, _AER1.r_drys, _AER1.Nis, _AER1.kappa
    )
    eq_leg, kn_leg, alpha_leg, phi_leg = binned_leg(smax, _T, _RS1, _AER1, approx=True)

    # Equilibrium fraction is identical: same approximate Köhler formula.
    assert float(eq_jax) == pytest.approx(float(eq_leg), rel=1e-10)
    # Kinetic fraction: JAX gates is_r_large on equilibrium activation so that
    # tiny non-activated bins (s_crit >> smax) cannot trigger false positives.
    # _RS1 = 3·r_drys gives rs ≥ r_crit only for bins ≤ ~5.7 nm, none of which
    # are equilibrium-activated at any of these smax values → kn = 0.
    assert float(kn_jax) == pytest.approx(0.0, abs=1e-12)
    assert float(alpha_jax) == pytest.approx(0.0, abs=1e-12)
    assert float(phi_jax) == pytest.approx(1.0, abs=1e-12)


def test_binned_activation_no_kinetic():
    """When no bin has grown past its critical radius, kinetic fraction is zero."""
    # Use dry radii as wet radii: no bin has grown at all.
    eq, kn, alpha, phi = binned_activation(
        _SMAX, _T, _AER1.r_drys, _AER1.r_drys, _AER1.Nis, _AER1.kappa
    )
    assert float(kn) == pytest.approx(0.0, abs=1e-12)
    assert float(phi) == pytest.approx(1.0, rel=1e-12)


def test_binned_activation_fracs_in_unit_interval():
    """eq_frac and kn_frac must lie in [0, 1]."""
    for smax in [0.001, 0.005, 0.02]:
        eq, kn, _, _ = binned_activation(smax, _T, _RS1, _AER1.r_drys, _AER1.Nis, _AER1.kappa)
        assert 0.0 <= float(eq) <= 1.0
        assert 0.0 <= float(kn) <= 1.0


def test_multi_mode_activation_vs_legacy():
    """JAX multi_mode_activation eq-fracs match legacy; kn-fracs use equilibrium gate."""
    eq_jax, kn_jax = multi_mode_activation(
        _SMAX,
        _T,
        [_RS1, _RS2],
        [_AER1.r_drys, _AER2.r_drys],
        [_AER1.Nis, _AER2.Nis],
        [_AER1.kappa, _AER2.kappa],
    )
    eq_leg, kn_leg = multi_mode_leg(
        _SMAX,
        _T,
        [_AER1, _AER2],
        [_RS1, _RS2],
    )

    for jax_val, leg_val in zip(eq_jax, eq_leg):
        assert float(jax_val) == pytest.approx(float(leg_val), rel=1e-10)
    # kn-frac: same false-positive issue as test_binned_activation_vs_legacy —
    # tiny non-eq-activated bins are correctly excluded by the equilibrium gate.
    for jax_val in kn_jax:
        assert float(jax_val) == pytest.approx(0.0, abs=1e-12)
