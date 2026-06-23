"""Cornerstone test: the JAX RHS must reproduce the frozen numba ``parcel_ode_sys``.

For every scenario we evaluate :func:`pyrcel.parcel_aux_jax.parcel_ode_sys` at the
exact ``(y, args)`` inputs the oracle used and compare against the frozen numba
``dy/dt``. This isolates the right-hand side from the ODE solver and catches any
transcription error immediately — and it runs with **no numba/Assimulo** present.

Tolerances: the port is elementwise identical up to float64 round-off (measured worst
absolute diff ~1e-14 across all fixtures), so we use ``rtol=1e-9, atol=1e-12``. ``atol``
sits comfortably above round-off yet far below any physical tendency, so genuine bugs
still fail loudly.

A second layer uses Hypothesis to assert structural identities and finiteness over
randomized states.
"""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from conftest import N_STATE_VARS, load_fixture
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra import numpy as hnp

from pyrcel.parcel_aux_jax import parcel_ode_sys, parcel_ode_sys_jit

RTOL = 1e-9
ATOL = 1e-12


def _args(oracle):
    return (
        jnp.asarray(oracle["r_drys"]),
        jnp.asarray(oracle["Nis"]),
        jnp.asarray(oracle["kappas"]),
        float(oracle["accom"]),
        float(oracle["V"]),
    )


def _eval_all(oracle):
    args = _args(oracle)
    Y = jnp.asarray(oracle["rhs_Y"])
    got = jax.vmap(lambda y: parcel_ode_sys(0.0, y, args))(Y)
    return np.asarray(got)


# --- Equivalence against frozen numba output -------------------------------------


def test_rhs_matches_numba(oracle):
    got = _eval_all(oracle)
    np.testing.assert_allclose(got, oracle["rhs_dYdt"], rtol=RTOL, atol=ATOL)


def test_jit_matches_eager(oracle):
    args = _args(oracle)
    y = jnp.asarray(oracle["rhs_Y"][0])
    got_jit = np.asarray(parcel_ode_sys_jit(0.0, y, args))
    got_eager = np.asarray(parcel_ode_sys(0.0, y, args))
    np.testing.assert_allclose(got_jit, got_eager, rtol=1e-12, atol=1e-15)


# --- Structural identities on the (physical + random) fixture states -------------


def test_structural_identities(oracle):
    V = float(oracle["V"])
    got = _eval_all(oracle)
    # dz/dt == V exactly
    np.testing.assert_allclose(got[:, 0], V, rtol=0, atol=0)
    # dwi/dt == 0 (no ice)
    np.testing.assert_array_equal(got[:, 5], 0.0)
    # mass balance: dwv/dt == -dwc/dt
    np.testing.assert_allclose(got[:, 3], -got[:, 4], rtol=1e-12, atol=0)


# --- Property-based: finiteness + identities over randomized states --------------

_F = load_fixture("simple_sulfate")
_NR = int(_F["nr"])
_RDRY = np.asarray(_F["r_drys"])
_KAP = np.asarray(_F["kappas"])
_NIS = np.asarray(_F["Nis"])


@settings(max_examples=60, deadline=None)
@given(
    T=st.floats(250.0, 305.0),
    P=st.floats(5e4, 1e5),
    S=st.floats(-0.05, 0.02),
    wv=st.floats(1e-3, 3e-2),
    log_growth=hnp.arrays(
        np.float64,
        (_NR,),
        elements=st.floats(0.0, 2.0),  # factor in [2, 101]
    ),
)
def test_rhs_finite_and_structural_random(T, P, S, wv, log_growth):
    rs = _RDRY * (1.0 + 10.0**log_growth)  # guarantees r > r_dry
    y = np.concatenate([[100.0, P, T, wv, 1e-4, 0.0, S], rs])
    args = (jnp.asarray(_RDRY), jnp.asarray(_NIS), jnp.asarray(_KAP), 1.0, 1.0)
    d = np.asarray(parcel_ode_sys(0.0, jnp.asarray(y), args))

    assert np.all(np.isfinite(d)), "non-finite RHS"
    assert d[0] == 1.0  # dz/dt == V
    assert d[5] == 0.0  # dwi/dt == 0
    np.testing.assert_allclose(d[3], -d[4], rtol=1e-10, atol=0)  # mass balance
    assert d.shape == (N_STATE_VARS + _NR,)


@settings(max_examples=40, deadline=None)
@given(V=st.floats(0.05, 10.0))
def test_dzdt_equals_updraft(V):
    y = _F["y0"]
    args = (jnp.asarray(_RDRY), jnp.asarray(_NIS), jnp.asarray(_KAP), 1.0, V)
    d = np.asarray(parcel_ode_sys(0.0, jnp.asarray(y), args))
    assert d[0] == pytest.approx(V, rel=0, abs=0)
