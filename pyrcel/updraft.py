"""Updraft (``V``) abstraction for the v2 parcel model (design doc §5.4, §6.7).

The master model passes ``V`` either as a constant or as a Python callable closed over
state, which makes the RHS impure and awkward to ``jit``/``grad``. Here ``V`` is a small
[equinox.Module][] -- an immutable, pytree-registered object -- so it threads cleanly
through ``jit``/``vmap``/``grad`` (its array fields become differentiable leaves) while
keeping the vector field a pure function of time.

* [ConstantV][pyrcel.updraft.ConstantV] -- a fixed updraft speed (the common case).
* [InterpolatedUpdraft][pyrcel.updraft.InterpolatedUpdraft] -- a tabulated ``V(t)`` profile,
piecewise-linearly
  interpolated with `jax.numpy.interp` (grad-friendly; no extra dependency, so we
  avoid pulling in ``interpax`` for the linear case).

All updrafts are callables ``V(t) -> speed`` in m/s, which is exactly what
[pyrcel.parcel_aux.parcel_ode_sys][] consumes.
"""

from __future__ import annotations

import equinox as eqx
import jax

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402
from jax import Array  # noqa: E402
from jax.typing import ArrayLike  # noqa: E402


class AbstractUpdraft(eqx.Module):
    """Base class for updraft models. Subclasses implement ``__call__(t) -> V``."""

    def __call__(self, t: ArrayLike) -> Array:  # pragma: no cover - abstract
        raise NotImplementedError


class ConstantV(AbstractUpdraft):
    """A constant updraft speed ``V`` (m/s)."""

    V: jax.Array

    def __init__(self, V: ArrayLike) -> None:
        # Store as a float64 scalar array so it is a differentiable pytree leaf.
        self.V = jnp.asarray(V, dtype=jnp.float64)

    def __call__(self, t: ArrayLike) -> Array:
        return self.V


class InterpolatedUpdraft(AbstractUpdraft):
    """A tabulated, piecewise-linear ``V(t)`` profile.

    Parameters
    ----------
    ts : array
        Strictly increasing knot times (s).
    vs : array
        Updraft speeds (m/s) at ``ts``. Outside ``[ts[0], ts[-1]]`` the endpoint values
        are held constant (the `jax.numpy.interp` default).
    """

    ts: jax.Array
    vs: jax.Array

    def __init__(self, ts: ArrayLike, vs: ArrayLike) -> None:
        self.ts = jnp.asarray(ts, dtype=jnp.float64)
        self.vs = jnp.asarray(vs, dtype=jnp.float64)

    def __call__(self, t: ArrayLike) -> Array:
        return jnp.interp(t, self.ts, self.vs)


def as_updraft(V) -> AbstractUpdraft:
    """Coerce ``V`` to an [AbstractUpdraft][pyrcel.updraft.AbstractUpdraft].

    Accepts an existing updraft (returned as-is) or a scalar (wrapped in
    [ConstantV][pyrcel.updraft.ConstantV]). Plain Python callables are *not* accepted here: model an
    arbitrary profile as an [AbstractUpdraft][pyrcel.updraft.AbstractUpdraft] so the vector field
    stays a pure,
    ``jit``-able pytree.

    Parameters
    ----------
    V : float or AbstractUpdraft
        Updraft speed (m/s) or an existing [AbstractUpdraft][pyrcel.updraft.AbstractUpdraft]
        instance.

    Returns
    -------
    AbstractUpdraft
        The original updraft, or ``V`` wrapped in [ConstantV][pyrcel.updraft.ConstantV].

    Raises
    ------
    TypeError
        If ``V`` is a plain Python callable; use
        [InterpolatedUpdraft][pyrcel.updraft.InterpolatedUpdraft] instead.
    """
    if isinstance(V, AbstractUpdraft):
        return V
    if callable(V):
        raise TypeError(
            "pass V as an AbstractUpdraft (e.g. ConstantV / InterpolatedUpdraft), "
            "not a bare Python callable"
        )
    return ConstantV(V)
