"""Equilibration equivalence + property tests for :mod:`pyrcel.equilibrate_jax`.

Design doc §7.2b: the v2 ``optimistix`` equilibration must reproduce ``master``'s
frozen initial state ``y0`` (the ``scipy.bisect`` result) to a tight tolerance.

Layers:

1. **Equivalence** vs the frozen oracle ``y0`` for every scenario (radii, wv0, wc0,
   and the verbatim bulk slots).
2. **Physical correctness**: the solved radii satisfy ``Seq(r0) == S0`` and lie in
   the stable branch ``(r_dry, r_crit)``.
3. **Property-based** (Hypothesis): root validity, bracketing, and monotonicity of
   the equilibrium radius in ``S0`` over randomized inputs; plus ``kohler_crit``
   correctness (``dSeq/dr == 0`` at the peak).
"""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pyrcel.equilibrate_jax import (
    equilibrate_initial_state,
    equilibrate_radii,
    kohler_crit,
)
from pyrcel.thermo_jax import Seq

# Evidence-based tolerances: worst observed radius/wc rel-diff vs master ~1e-14.
R_RTOL = 1e-10
R_ATOL = 1e-15
WC_RTOL = 1e-9


# --- §7.2b: equivalence against frozen master y0 ---------------------------------

def test_y0_matches_master(oracle):
    T0 = float(oracle["T0"])
    S0 = float(oracle["S0"])
    P0 = float(oracle["P0"])
    y0_ref = np.asarray(oracle["y0"])
    y0_v2 = np.asarray(
        equilibrate_initial_state(T0, S0, P0, oracle["r_drys"], oracle["kappas"], oracle["Nis"])
    )

    # Equilibrium wet radii (the substantive output of the solve).
    np.testing.assert_allclose(y0_v2[7:], y0_ref[7:], rtol=R_RTOL, atol=R_ATOL)
    # Derived water content.
    np.testing.assert_allclose(y0_v2[3], y0_ref[3], rtol=1e-12)  # wv0 (analytic)
    np.testing.assert_allclose(y0_v2[4], y0_ref[4], rtol=WC_RTOL)  # wc0
    # Verbatim bulk slots: z, P, T, wi, S.
    np.testing.assert_array_equal(y0_v2[[0, 1, 2, 5, 6]], y0_ref[[0, 1, 2, 5, 6]])


def test_equilibrium_residual_is_zero(oracle):
    """Each solved radius must satisfy Seq(r0, r_dry, T0, kappa) == S0."""
    T0 = float(oracle["T0"])
    S0 = float(oracle["S0"])
    r0 = equilibrate_radii(T0, S0, oracle["r_drys"], oracle["kappas"])
    resid = np.asarray(Seq(r0, jnp.asarray(oracle["r_drys"]), T0, jnp.asarray(oracle["kappas"]))) - S0
    np.testing.assert_allclose(resid, 0.0, atol=1e-10)


def test_radii_in_stable_branch(oracle):
    T0 = float(oracle["T0"])
    S0 = float(oracle["S0"])
    r_drys = np.asarray(oracle["r_drys"])
    kappas = np.asarray(oracle["kappas"])
    r0 = np.asarray(equilibrate_radii(T0, S0, r_drys, kappas))
    r_crit = np.asarray(jax.vmap(lambda rd, k: kohler_crit(T0, rd, k))(jnp.asarray(r_drys), jnp.asarray(kappas)))
    assert np.all(r0 > r_drys), "equilibrium radius must exceed dry radius"
    assert np.all(r0 < r_crit), "equilibrium radius must be sub-critical (stable branch)"


# --- Property-based tests --------------------------------------------------------

_T = st.floats(250.0, 310.0)
_rd = st.floats(2e-9, 1e-6)
_kappa = st.floats(0.05, 1.4)
_S0 = st.floats(-0.2, -1e-4)
HSET = settings(max_examples=60, deadline=None)

# JIT scalar helpers so Hypothesis reuses one compiled solve across all examples
# (the nested bisections are expensive to re-trace eagerly per example).
_solve_one = jax.jit(
    lambda T0, S0, rd, kappa: equilibrate_radii(T0, S0, jnp.array([rd]), jnp.array([kappa]))[0]
)
_crit_one = jax.jit(lambda T0, rd, kappa: kohler_crit(T0, rd, kappa))


@HSET
@given(T0=_T, rd=_rd, kappa=_kappa, S0=_S0)
def test_equilibrium_root_valid(T0, rd, kappa, S0):
    r0 = float(_solve_one(T0, S0, rd, kappa))
    assert r0 > rd
    resid = float(Seq(r0, rd, T0, kappa)) - S0
    assert abs(resid) < 1e-8


@HSET
@given(T0=_T, rd=_rd, kappa=_kappa)
def test_kohler_crit_is_peak(T0, rd, kappa):
    r_crit = float(_crit_one(T0, rd, kappa))
    assert r_crit > rd
    dseq = float(jax.grad(lambda r: Seq(r, rd, T0, kappa))(r_crit))
    # Derivative vanishes at the peak; scale is tiny near r_crit so use a small abs tol.
    assert abs(dseq) < 1e-3


@HSET
@given(T0=_T, rd=_rd, kappa=_kappa, S0=_S0, dS=st.floats(1e-4, 1e-2))
def test_radius_monotonic_in_S0(T0, rd, kappa, S0, dS):
    # A higher ambient supersaturation yields a larger equilibrium radius.
    S_hi = min(S0 + dS, -1e-6)
    r_lo = float(_solve_one(T0, S0, rd, kappa))
    r_hi = float(_solve_one(T0, S_hi, rd, kappa))
    assert r_hi >= r_lo - 1e-18
