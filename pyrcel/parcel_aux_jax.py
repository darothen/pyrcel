"""JAX implementation of the parcel-model right-hand side (vector field).

This is the v2 counterpart to the numba ``parcel_ode_sys`` in
``pyrcel/_parcel_aux_numba.py``. The per-bin growth loop (a ``numba.prange`` in the
original) is expressed here as vectorized array operations, which is both cleaner and
naturally differentiable/``vmap``-able.

The signature follows the diffrax convention ``f(t, y, args)`` so the function can be
wrapped directly in ``diffrax.ODETerm``. ``args`` is the parameter pytree

    args = (r_drys, Nis, kappas, accom, V)

where ``V`` is the updraft speed: a scalar today, or (in a later phase) a callable /
``equinox.Module`` of time. See ``docs/design/jax-diffrax-migration.md`` §4.2.

Equivalence with the numba derivative is verified against frozen fixtures in
``tests/test_rhs_equivalence.py``.
"""

from __future__ import annotations

import math

import equinox as eqx
import jax

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402

from . import constants as c  # noqa: E402
from .thermo_jax import Seq, dv, es, ka  # noqa: E402
from .updraft import AbstractUpdraft, as_updraft  # noqa: E402

PI = math.pi

#: Number of bulk (non-radius) state variables; mirrors ``pyrcel.constants``.
N_STATE_VARS = c.N_STATE_VARS


def parcel_ode_sys(t, y, args):
    """Instantaneous tendency of the parcel-model state vector.

    Parameters
    ----------
    t : float
        Simulation time, s. Only used if ``V`` is time-dependent.
    y : array, shape ``(N_STATE_VARS + nr,)``
        State vector ``[z, P, T, wv, wc, wi, S, r_0 ... r_{nr-1}]``.
    args : tuple
        ``(r_drys, Nis, kappas, accom, V)``.

    Returns
    -------
    array, shape ``(N_STATE_VARS + nr,)``
        Time derivative ``dy/dt``.
    """
    r_drys, Nis, kappas, accom, V = args
    V_t = V(t) if callable(V) else V

    P = y[1]
    T = y[2]
    wv = y[3]
    S = y[6]
    rs = y[N_STATE_VARS:]

    pv_sat = es(T - 273.15)  # saturation vapor pressure (T in Celsius)
    Tv = (1.0 + 0.61 * wv) * T  # virtual temperature
    e = (1.0 + S) * pv_sat  # ambient vapor pressure
    rho_air = P / c.Rd / Tv
    rho_air_dry = (P - e) / c.Rd / T

    # Per-bin condensational growth (vectorized over all aerosol bins).
    dv_r = dv(T, rs, P, accom)
    ka_r = ka(T, rho_air, rs)
    G_a = (c.rho_w * c.R * T) / (pv_sat * dv_r * c.Mw)
    G_b = (c.L * c.rho_w * ((c.L * c.Mw / (c.R * T)) - 1.0)) / (ka_r * T)
    G = 1.0 / (G_a + G_b)
    delta_S = S - Seq(rs, r_drys, T, kappas)
    drs_dt = (G / rs) * delta_S

    # Liquid water tendency from droplet growth.
    dwc_dt = (4.0 * PI * c.rho_w / rho_air_dry) * jnp.sum(Nis * rs * rs * drs_dt)
    dwi_dt = 0.0  # no ice/freezing yet
    dwv_dt = -1.0 * (dwc_dt + dwi_dt)  # mass balance

    dP_dt = -1.0 * rho_air * c.g * V_t
    dT_dt = -c.g * V_t / c.Cp - c.L * dwv_dt / c.Cp
    dz_dt = V_t

    # Supersaturation tendency, Ghan (2011) form.
    alpha = (c.g * c.Mw * c.L) / (c.Cp * c.R * (T**2)) - (c.g * c.Ma) / (c.R * T)
    gamma = (P * c.Ma) / (c.Mw * pv_sat) + (c.Mw * c.L * c.L) / (c.Cp * c.R * T * T)
    dS_dt = alpha * V_t - gamma * dwc_dt

    bulk = jnp.array([dz_dt, dP_dt, dT_dt, dwv_dt, dwc_dt, dwi_dt, dS_dt])
    return jnp.concatenate([bulk, drs_dt])


#: JIT-compiled alias for performance-sensitive callers.
parcel_ode_sys_jit = jax.jit(parcel_ode_sys)


class ParcelVectorField(eqx.Module):
    """The parcel vector field as a typed, differentiable pytree (design doc §5.4).

    **Mental model: state vs. parameters vs. the rule.** An ODE is
    ``dy/dt = f(t, y, theta)``, with three distinct roles:

    * ``y`` -- the *state* that evolves (``[z, P, T, wv, wc, wi, S, r_0...]``). The
      *solver* owns it and threads it forward step by step.
    * ``theta`` -- the *parameters* that configure the dynamics but do not themselves
      evolve: ``(r_drys, Nis, kappas, accom, V)``.
    * ``f`` -- the *rule* mapping the current ``(t, y)`` to the tendency ``dy/dt``.

    This class **is** ``f`` (the rule), so it holds ``theta`` as fields -- *not* ``y``.
    The evolving state is consumed as a call argument, ``__call__(self, t, y)``: every
    solver step ``diffrax`` calls ``field(t, y_current)`` and gets back ``dy/dt``. Think
    of it as a physics engine configured once for one experiment (this aerosol
    population, this updraft); the state ``y`` is the thing the engine acts on, handed to
    it fresh each tick. Storing ``y`` here would be a category error (and impossible
    anyway -- ``eqx.Module`` is immutable).

    **Why wrap the fixed parameters in a Module at all** (the plain
    :func:`parcel_ode_sys` ``(t, y, args)`` function is still the primary path):

    1. *Named fields* instead of a positional 5-tuple, removing the ``args``/``rhs_args``
       indexing footgun in the master code.
    2. *It is a JAX pytree*, so its array leaves (``r_drys``, ``Nis``, ``kappas``,
       ``accom``, and ``V``'s parameters) are visible to ``jit``/``vmap``/``grad``. That
       is what makes **parameter** sensitivities clean: ``eqx.filter_grad`` over a field
       gives ``d S_max / d accom`` or ``d S_max / d V`` directly. (Gradients w.r.t. the
       *state's initial value* ``y0`` need no Module -- just ``jax.grad`` over the plain
       array.) ``accom`` is stored as an array, not a Python float, precisely so it is a
       differentiable leaf rather than a compile-time constant.

    The instance is callable as ``field(t, y)`` (diffrax ``ODETerm`` convention, with an
    optional ignored ``args``) and exposes :pyattr:`args` to feed the tuple-based
    integrator helpers in :mod:`pyrcel.integrator_diffrax`.
    """

    r_drys: jax.Array
    Nis: jax.Array
    kappas: jax.Array
    accom: jax.Array
    V: AbstractUpdraft

    def __init__(self, r_drys, Nis, kappas, accom, V):
        self.r_drys = jnp.asarray(r_drys, dtype=jnp.float64)
        self.Nis = jnp.asarray(Nis, dtype=jnp.float64)
        self.kappas = jnp.asarray(kappas, dtype=jnp.float64)
        self.accom = jnp.asarray(accom, dtype=jnp.float64)
        self.V = as_updraft(V)

    @classmethod
    def from_args(cls, args) -> "ParcelVectorField":
        """Build from the positional ``(r_drys, Nis, kappas, accom, V)`` tuple."""
        r_drys, Nis, kappas, accom, V = args
        return cls(r_drys, Nis, kappas, accom, V)

    @property
    def args(self) -> tuple:
        """The ``(r_drys, Nis, kappas, accom, V)`` tuple for the integrator helpers."""
        return (self.r_drys, self.Nis, self.kappas, self.accom, self.V)

    def __call__(self, t, y, args=None):
        return parcel_ode_sys(t, y, self.args)
