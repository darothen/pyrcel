"""Updraft abstraction + Equinox vector field (design §5.4, §6.7; Phase 5).

Checks that:

* :class:`ConstantV` / :class:`InterpolatedUpdraft` behave as plain callables of time,
* wrapping a scalar updraft in :class:`ConstantV` is *exactly* equivalent to passing the
  bare scalar (RHS and full integration), so the abstraction is zero-cost numerically,
* :class:`ParcelVectorField` round-trips through its ``args`` tuple, and
* a genuinely time-varying ``V(t)`` integrates and respects physical invariants.
"""

from __future__ import annotations

import jax.numpy as jnp
import numpy as np
import pytest

from pyrcel.integrator_diffrax import integrate_parcel, integrate_parcel_arrays
from pyrcel.parcel_aux_jax import ParcelVectorField, parcel_ode_sys
from pyrcel.updraft import ConstantV, InterpolatedUpdraft, as_updraft


def _args(oracle, V):
    return (
        jnp.asarray(oracle["r_drys"]),
        jnp.asarray(oracle["Nis"]),
        jnp.asarray(oracle["kappas"]),
        float(oracle["accom"]),
        V,
    )


# --- updraft callables -----------------------------------------------------------


def test_constantv_is_constant():
    V = ConstantV(2.5)
    assert float(V(0.0)) == 2.5
    assert float(V(123.0)) == 2.5


def test_interpolated_updraft():
    V = InterpolatedUpdraft(ts=[0.0, 10.0, 20.0], vs=[1.0, 3.0, 3.0])
    assert float(V(0.0)) == pytest.approx(1.0)
    assert float(V(5.0)) == pytest.approx(2.0)  # linear midpoint
    assert float(V(15.0)) == pytest.approx(3.0)
    # endpoints held constant outside the knot range
    assert float(V(-5.0)) == pytest.approx(1.0)
    assert float(V(99.0)) == pytest.approx(3.0)


def test_interpolated_constant_equals_constantv():
    Vc = ConstantV(1.7)
    Vi = InterpolatedUpdraft(ts=[0.0, 100.0], vs=[1.7, 1.7])
    for t in (0.0, 13.0, 88.0):
        assert float(Vi(t)) == pytest.approx(float(Vc(t)))


def test_as_updraft():
    assert isinstance(as_updraft(1.0), ConstantV)
    v = ConstantV(1.0)
    assert as_updraft(v) is v
    with pytest.raises(TypeError):
        as_updraft(lambda t: 1.0)


# --- zero-cost equivalence with the scalar path ----------------------------------


def test_constantv_matches_scalar_rhs(oracle):
    y0 = jnp.asarray(oracle["y0"])
    V = float(oracle["V"])
    f_scalar = parcel_ode_sys(0.0, y0, _args(oracle, V))
    f_const = parcel_ode_sys(0.0, y0, _args(oracle, ConstantV(V)))
    np.testing.assert_array_equal(np.asarray(f_scalar), np.asarray(f_const))


@pytest.mark.slow
def test_constantv_matches_scalar_integration(oracle):
    y0 = jnp.asarray(oracle["y0"])
    ts = jnp.asarray(oracle["traj_t"])
    V = float(oracle["V"])
    ys_scalar = integrate_parcel(y0, _args(oracle, V), ts).ys
    ys_const = integrate_parcel(y0, _args(oracle, ConstantV(V)), ts).ys
    np.testing.assert_array_equal(np.asarray(ys_scalar), np.asarray(ys_const))


# --- Equinox vector field --------------------------------------------------------


@pytest.mark.slow
def test_vectorfield_roundtrip_and_call(oracle):
    y0 = jnp.asarray(oracle["y0"])
    V = float(oracle["V"])
    args = _args(oracle, ConstantV(V))
    field = ParcelVectorField.from_args(args)
    # args property reproduces the tuple (arrays equal, V identical type)
    assert isinstance(field.V, ConstantV)
    np.testing.assert_array_equal(np.asarray(field.r_drys), np.asarray(args[0]))
    # __call__ matches the bare RHS
    f_field = field(0.0, y0)
    f_direct = parcel_ode_sys(0.0, y0, args)
    np.testing.assert_array_equal(np.asarray(f_field), np.asarray(f_direct))
    # and integrates identically via field.args
    ts = jnp.asarray(oracle["traj_t"])
    ys_field = integrate_parcel(y0, field.args, ts).ys
    ys_direct = integrate_parcel(y0, args, ts).ys
    np.testing.assert_array_equal(np.asarray(ys_field), np.asarray(ys_direct))


# --- time-varying V(t) -----------------------------------------------------------


@pytest.mark.slow
def test_time_varying_updraft_invariants():
    # Use the simple_sulfate fixture's initial state with a ramped updraft.
    from conftest import load_fixture

    oracle = load_fixture("simple_sulfate")
    y0 = jnp.asarray(oracle["y0"])
    Vt = InterpolatedUpdraft(ts=[0.0, 40.0, 120.0], vs=[0.3, 1.2, 1.2])
    args = _args(oracle, Vt)
    ts = jnp.asarray(np.arange(0.0, 120.0, 1.0))
    _, ys, success = integrate_parcel_arrays(y0, args, ts)
    ys = np.asarray(ys)

    assert success
    assert np.all(np.isfinite(ys))
    # z strictly increasing (dz/dt = V(t) > 0 throughout)
    assert np.all(np.diff(ys[:, 0]) > 0)
    # total condensed + vapour water conserved (no ice): wv + wc ~ const
    water = ys[:, 3] + ys[:, 4]
    assert np.max(np.abs(water - water[0])) < 1e-9
    # supersaturation has a single interior peak
    assert 0 < int(np.argmax(ys[:, 6])) < len(ys) - 1
