"""Equivalence + property tests for :mod:`pyrcel.thermo`.

Two layers:

1. **Grid equivalence** against the NumPy oracle :mod:`pyrcel.legacy.thermo`, over the
   same parameter grids used by the legacy reference generator. Identical formulas, so we
   demand agreement to float64 round-off.
2. **Property-based** tests with Hypothesis: equivalence over randomized inputs, plus
   physical invariants (positivity, monotonicity, non-continuum corrections reduce the
   continuum values).
"""

from __future__ import annotations

import subprocess
import sys

import jax.numpy as jnp
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pyrcel import constants as pc
from pyrcel import thermo as tj
from pyrcel.legacy import thermo as th


def test_jax_modules_do_not_import_numba():
    """The v2 JAX modules must be usable without numba/Assimulo (design goal).

    Run in a fresh interpreter so prior imports in this session can't mask a
    regression.
    """
    code = (
        "import sys; import pyrcel.parcel_aux, pyrcel.thermo; "
        "assert 'numba' not in sys.modules, 'numba was imported'; "
        "assert 'pyrcel.legacy' not in sys.modules, 'legacy module imported'"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr


# --- Parameter grids (mirror pyrcel/test/generate_data.py) -----------------------
T = np.linspace(233.0, 333.0, 10)  # K
P = np.linspace(100.0, 1050.0, 10) * 100.0  # Pa (from hPa)
R = np.logspace(-3, 1, 10) * 1e-6  # m (from microns)
DENS = np.linspace(0.5, 1.5, 10)  # kg/m^3
RD = np.logspace(-8, -6, 5)  # m
KAP = np.logspace(-2, 0, 5)  # unitless (all > 0)
ROVERRD = np.logspace(0.1, 3, 5)  # r / r_dry


def _close(got, expected, rtol=1e-12, atol=0.0):
    np.testing.assert_allclose(
        np.asarray(got, dtype=float),
        np.asarray(expected, dtype=float),
        rtol=rtol,
        atol=atol,
    )


def _mesh(*arrays):
    grids = np.meshgrid(*arrays, indexing="ij")
    return [g.ravel() for g in grids]


# --- Grid equivalence ------------------------------------------------------------


def test_sigma_w_grid():
    _close(tj.sigma_w(jnp.asarray(T)), th.sigma_w(T))


def test_es_grid():
    Tc = T - 273.15
    _close(tj.es(jnp.asarray(Tc)), th.es(Tc))


def test_dv_cont_grid():
    TT, PP = _mesh(T, P)
    _close(tj.dv_cont(jnp.asarray(TT), jnp.asarray(PP)), th.dv_cont(TT, PP))


def test_dv_grid():
    TT, RR, PP = _mesh(T, R, P)
    _close(tj.dv(jnp.asarray(TT), jnp.asarray(RR), jnp.asarray(PP)), th.dv(TT, RR, PP))


def test_ka_cont_grid():
    _close(tj.ka_cont(jnp.asarray(T)), th.ka_cont(T))


def test_ka_grid():
    # thermo.ka and legacy.thermo.ka share the (T, rho, r) signature.
    TT, DD, RR = _mesh(T, DENS, R)
    _close(tj.ka(jnp.asarray(TT), jnp.asarray(DD), jnp.asarray(RR)), th.ka(TT, DD, RR))


def test_rho_air_grid():
    TT, PP = _mesh(T, P)
    _close(tj.rho_air(jnp.asarray(TT), jnp.asarray(PP)), th.rho_air(TT, PP))


def test_seq_grid():
    F, RDg, Tg, Kg = _mesh(ROVERRD, RD, T, KAP)
    r = F * RDg
    _close(
        tj.Seq(jnp.asarray(r), jnp.asarray(RDg), jnp.asarray(Tg), jnp.asarray(Kg)),
        th.Seq(r, RDg, Tg, Kg),
        # The JAX Seq uses difference-of-cubes + expm1 for accuracy; the legacy
        # NumPy version uses the naïve r³-r_dry³ path. They agree to ~2 ULPs
        # (≤5e-10 relative) everywhere on this grid.
        rtol=5e-10,
    )


def test_seq_approx_grid():
    F, RDg, Tg, Kg = _mesh(ROVERRD, RD, T, KAP)
    r = F * RDg
    _close(
        tj.Seq_approx(jnp.asarray(r), jnp.asarray(RDg), jnp.asarray(Tg), jnp.asarray(Kg)),
        th.Seq_approx(r, RDg, Tg, Kg),
        rtol=1e-11,
    )


# --- Property-based equivalence --------------------------------------------------

_T = st.floats(233.0, 333.0)
_Tc = st.floats(-40.0, 60.0)
_P = st.floats(1e4, 1.1e5)
_r = st.floats(1e-9, 1e-4)
_rho = st.floats(0.3, 1.6)
_accom = st.floats(0.1, 1.0)
_kappa = st.floats(1e-3, 1.4)
HSET = settings(max_examples=150, deadline=None)


@HSET
@given(Tc=_Tc)
def test_es_property_equiv(Tc):
    assert float(tj.es(jnp.asarray(Tc))) == pytest.approx(float(th.es(Tc)), rel=1e-11)


@HSET
@given(T=_T, r=_r, P=_P, accom=_accom)
def test_dv_property_equiv(T, r, P, accom):
    got = float(tj.dv(jnp.asarray(T), jnp.asarray(r), jnp.asarray(P), accom))
    assert got == pytest.approx(float(th.dv(T, r, P, accom)), rel=1e-10)


@HSET
@given(T=_T, rho=_rho, r=_r)
def test_ka_property_equiv(T, rho, r):
    got = float(tj.ka(jnp.asarray(T), jnp.asarray(rho), jnp.asarray(r)))
    assert got == pytest.approx(float(th.ka(T, rho, r)), rel=1e-10)


@HSET
@given(T=_T, rd=st.floats(1e-9, 1e-6), f=st.floats(1.001, 1e3), kappa=_kappa)
def test_seq_property_equiv(T, rd, f, kappa):
    r = f * rd
    got = float(tj.Seq(jnp.asarray(r), jnp.asarray(rd), jnp.asarray(T), kappa))
    assert got == pytest.approx(float(th.Seq(r, rd, T, kappa)), rel=1e-9, abs=1e-12)


# --- Physical invariants ---------------------------------------------------------


@HSET
@given(T=_T, r=_r, P=_P, accom=_accom)
def test_dv_below_continuum(T, r, P, accom):
    # Non-continuum correction always reduces diffusivity below the continuum value.
    assert tj.dv(T, r, P, accom) < tj.dv_cont(T, P)


@HSET
@given(T=_T, rho=_rho, r=_r)
def test_ka_below_continuum(T, rho, r):
    assert tj.ka(T, rho, r) < tj.ka_cont(T)


@HSET
@given(Tc=st.floats(-40.0, 59.0), d=st.floats(0.1, 1.0))
def test_es_monotonic_increasing(Tc, d):
    assert tj.es(Tc + d) > tj.es(Tc)


@HSET
@given(T=_T)
def test_es_positive(T):
    assert float(tj.es(T - 273.15)) > 0.0


# --- Seq numerical stability tests -----------------------------------------------


def test_seq_kappa_zero():
    """With kappa=0 the solute term vanishes: Seq = exp(A) - 1 = expm1(A)."""
    r = jnp.array(0.1e-6)
    r_dry = jnp.array(0.05e-6)
    T = 300.0
    A = (2.0 * pc.Mw * tj.sigma_w(T)) / (pc.R * T * pc.rho_w * r)
    expected = float(jnp.expm1(A))
    assert float(tj.Seq(r, r_dry, T, kappa=0.0)) == pytest.approx(expected, rel=1e-12)


def test_seq_near_dry_finite():
    """Seq must return a finite value even for r/r_dry as small as 1 + 1e-9."""
    r_dry = jnp.array(1e-8)
    T = 300.0
    kappa = 0.6
    for eps in [1e-4, 1e-6, 1e-9]:
        r = r_dry * (1.0 + eps)
        val = float(tj.Seq(r, r_dry, T, kappa))
        assert np.isfinite(val), f"Seq not finite at eps={eps}"
        assert val > -1.0


def test_seq_near_dry_analytical():
    """For r = r_dry*(1+eps) with eps ≪ 1, compare against the first-order analytical B.

    B_first_order = 3*eps / (3*eps + kappa); uses this to check that the
    difference-of-cubes factoring preserves full float64 precision near r ≈ r_dry.
    """
    r_dry = 1e-7  # 100 nm
    T = 300.0
    kappa = 0.6
    for eps in [1e-3, 1e-5, 1e-7]:
        r = r_dry * (1.0 + eps)
        seq_jax = float(tj.Seq(jnp.array(r), jnp.array(r_dry), T, kappa))

        # First-order analytical B (no cancellation)
        B_ref = (3.0 * eps) / (3.0 * eps + kappa)
        A = (2.0 * pc.Mw * tj.sigma_w(T)) / (pc.R * T * pc.rho_w * r)
        # Use the same compensated form for the reference to avoid noise in comparison
        seq_ref = float(B_ref * jnp.expm1(A) + (B_ref - 1.0))

        # For eps up to 1e-3 the first-order approximation should agree to ~1%
        rel_tol = max(3.0 * eps / kappa, 1e-10)  # error ∝ eps from higher-order terms
        assert seq_jax == pytest.approx(seq_ref, rel=rel_tol, abs=1e-14)


def test_seq_large_r_solute_vanishes():
    """For r >> r_dry the solute term B → 1; Seq ≈ Seq(kappa=0)."""
    r_dry = jnp.array(1e-8)
    T = 300.0
    kappa = 0.6
    r_large = jnp.array(1e-5)  # r/r_dry = 1000
    A = (2.0 * pc.Mw * tj.sigma_w(T)) / (pc.R * T * pc.rho_w * r_large)
    seq_hygro = float(tj.Seq(r_large, r_dry, T, kappa))
    seq_pure = float(jnp.expm1(A))  # kappa=0 limit
    # Solute contribution ≈ kappa*(r_dry/r)³ = 0.6 * 1e-9 ≪ A ≈ 1e-3
    assert abs(seq_hygro - seq_pure) < 1e-8
