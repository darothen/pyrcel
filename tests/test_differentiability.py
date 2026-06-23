"""Differentiability + ensemble smoke tests (design §5.3, §7.6; Phase 5).

The locked requirement is gradients *through the integration* w.r.t. the (equilibrated)
initial state ``y0`` -- a capability the numba/CVode stack cannot provide. Here we check
that ``jax.grad`` of ``S_max`` w.r.t. ``y0``, ``accom`` and ``V`` runs, is finite, and
agrees with a finite-difference estimate to a loose tolerance, and that a ``vmap``
ensemble over perturbed parcels matches a Python loop.

All tests use the ``simple_sulfate`` fixture's equilibrated ``y0`` on a fixed horizon
that spans the supersaturation peak (so ``S_max`` is interior and well-defined).
"""

from __future__ import annotations

import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
import pytest

from conftest import load_fixture
from pyrcel.integrator_diffrax import max_supersaturation
from pyrcel.parcel_aux_jax import ParcelVectorField
from pyrcel.updraft import ConstantV

pytestmark = pytest.mark.slow

_TS = jnp.asarray(np.arange(0.0, 80.0, 1.0))  # spans t_smax ~ 52 s for simple_sulfate


def _inputs():
    d = load_fixture("simple_sulfate")
    y0 = jnp.asarray(d["y0"])
    field = ParcelVectorField(
        d["r_drys"], d["Nis"], d["kappas"], float(d["accom"]), ConstantV(float(d["V"]))
    )
    return y0, field


def test_grad_smax_wrt_y0_finite():
    y0, field = _inputs()
    g = jax.grad(lambda yy: max_supersaturation(yy, field.args, _TS))(y0)
    assert g.shape == y0.shape
    assert np.all(np.isfinite(np.asarray(g)))
    assert float(jnp.linalg.norm(g)) > 0.0


def test_grad_smax_wrt_y0_finite_difference():
    # Directional derivative along the initial-supersaturation component (smooth).
    y0, field = _inputs()
    g = jax.grad(lambda yy: max_supersaturation(yy, field.args, _TS))(y0)
    e = jnp.zeros_like(y0).at[6].set(1.0)
    dir_grad = float(jnp.dot(g, e))

    def f(h):
        return float(max_supersaturation(y0 + h * e, field.args, _TS))

    eps = 1e-6
    fd = (f(eps) - f(-eps)) / (2 * eps)
    np.testing.assert_allclose(dir_grad, fd, rtol=1e-2)


def test_grad_smax_wrt_accom_finite_difference():
    y0, field = _inputs()

    def smax_of_accom(a):
        f = eqx.tree_at(lambda fl: fl.accom, field, jnp.asarray(a))
        return max_supersaturation(y0, f.args, _TS)

    a0 = float(field.accom)
    g = float(jax.grad(smax_of_accom)(a0))
    eps = 1e-4
    fd = float((smax_of_accom(a0 + eps) - smax_of_accom(a0 - eps)) / (2 * eps))
    assert np.isfinite(g)
    np.testing.assert_allclose(g, fd, rtol=2e-2)


def test_grad_smax_wrt_V_via_filter_grad():
    y0, field = _inputs()

    def smax_field(fl):
        return max_supersaturation(y0, fl.args, _TS)

    grads = eqx.filter_grad(smax_field)(field)
    gV = float(grads.V.V)
    assert np.isfinite(gV)
    assert np.all(np.isfinite(np.asarray(grads.r_drys)))

    # finite-difference cross-check on V
    v0 = float(field.V.V)

    def smax_of_V(v):
        fl = eqx.tree_at(lambda f: f.V.V, field, jnp.asarray(v))
        return max_supersaturation(y0, fl.args, _TS)

    eps = 1e-4
    fd = float((smax_of_V(v0 + eps) - smax_of_V(v0 - eps)) / (2 * eps))
    np.testing.assert_allclose(gV, fd, rtol=2e-2)


def test_vmap_ensemble_matches_loop():
    y0, field = _inputs()
    rng = np.random.default_rng(0)
    batch = 8
    Y0 = jnp.asarray(
        np.asarray(y0)[None, :] * (1.0 + 1e-4 * rng.standard_normal((batch, y0.shape[0])))
    )

    smax_of = lambda yy: max_supersaturation(yy, field.args, _TS)
    batched = jax.vmap(smax_of)(Y0)
    looped = jnp.asarray([smax_of(Y0[i]) for i in range(batch)])

    assert batched.shape == (batch,)
    assert np.all(np.isfinite(np.asarray(batched)))
    np.testing.assert_allclose(np.asarray(batched), np.asarray(looped), rtol=1e-10)


def test_jit_smax_matches_eager():
    y0, field = _inputs()
    eager = float(max_supersaturation(y0, field.args, _TS))
    jitted = float(jax.jit(lambda yy: max_supersaturation(yy, field.args, _TS))(y0))
    assert jitted == pytest.approx(eager, rel=1e-12)
