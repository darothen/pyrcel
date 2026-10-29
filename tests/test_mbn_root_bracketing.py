"""Standalone regressions for this fix on unmodified upstream."""

import jax.numpy as jnp
import pytest
from scipy.optimize import brentq

from pyrcel.activation._mbn2014 import _mbn_residual, mbn2014


def test_mbn_subfloor_root_and_residual():
    s = mbn2014(0.01, 283.15, 85000.0, [2.0], [1.8], [100.0], [0.6])[0]
    reference = brentq(
        lambda value: float(
            _mbn_residual(
                value,
                0.01,
                283.15,
                85000.0,
                jnp.array([2e-6]),
                jnp.array([1.8]),
                jnp.array([100e6]),
                jnp.array([0.6]),
                1.0,
            )
        ),
        1e-9,
        1e-4,
        xtol=1e-16,
    )
    assert float(s) == pytest.approx(reference, rel=1e-6)
    residual = _mbn_residual(
        s,
        0.01,
        283.15,
        85000.0,
        jnp.array([2e-6]),
        jnp.array([1.8]),
        jnp.array([100e6]),
        jnp.array([0.6]),
        1.0,
    )
    assert abs(float(residual)) < 1e-5


def test_mbn_out_of_domain_is_explicit():
    with pytest.raises(RuntimeError, match="bracket"):
        mbn2014(5.0, 283.15, 85000.0, [0.05], [1.8], [0.01], [0.6])
