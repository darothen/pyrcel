"""Standalone regressions for this fix on unmodified upstream."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from pyrcel.activation import arg2000


@pytest.mark.parametrize("accom", [0.5, 1.0])
def test_arg_accommodation_jit_gradient_and_batch(accom):
    def f(a):
        return arg2000(0.5, 283.15, 85000.0, [0.05], [2.0], [300.0], [0.6], accom=a)[0]

    grad = float(jax.jit(jax.grad(f))(jnp.array(accom)))
    h = 1e-5
    fd = float((f(accom + h) - f(accom - h)) / (2 * h))
    assert grad == pytest.approx(fd, rel=1e-5)
    batch = np.asarray(jax.jit(jax.vmap(f))(jnp.array([0.1, 0.5, 1.0])))
    assert np.isfinite(batch).all() and np.all(np.diff(batch) < 0)
