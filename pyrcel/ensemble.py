"""Vectorized (``vmap``-ed) parcel-model ensembles over updraft velocity.

A direct payoff of the JAX migration (design §5.3): an entire ensemble of parcels can be
integrated in a single compiled, batched call. The canonical use case is propagating a
*distribution of vertical velocities* through the model to estimate the resulting
distributions of peak supersaturation ``S_max`` and activated droplet number.

The core, :func:`smax_nact_ensemble`, ``vmap``s the event-based ``S_max`` solve over a
batch of updraft speeds. Because the equilibrated initial state ``y0`` does not depend on
``V``, it is computed once and shared across the batch. Each member runs to its *own*
supersaturation maximum (the ``dS/dt`` event), so members with very different ``V`` (hence
very different times-to-peak) are handled correctly within one ``vmap``.
"""

from __future__ import annotations

import jax
import numpy as np

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402

from . import constants as c  # noqa: E402
from .equilibrate_jax import equilibrate_initial_state, kohler_crit_approx  # noqa: E402
from .integrator_diffrax import STATE_RTOL, _solve_to_smax, atol_vector  # noqa: E402
from .updraft import ConstantV  # noqa: E402

__all__ = ["sample_gaussian_updrafts", "smax_nact_ensemble", "run_updraft_ensemble"]


def sample_gaussian_updrafts(mean, std, n, *, seed=0, v_min=1e-2):
    """Draw ``n`` updraft speeds from ``Normal(mean, std)`` (m/s).

    Samples are clipped at ``v_min`` so every member is a physical (positive) updraft;
    with a low-mean/high-std request this truncation matters, so the returned array is
    what was actually simulated.
    """
    rng = np.random.default_rng(seed)
    return np.clip(rng.normal(mean, std, size=int(n)), v_min, None)


def smax_nact_ensemble(
    y0,
    r_drys,
    Nis,
    kappas,
    accom,
    V_samples,
    t_end,
    *,
    rtol=STATE_RTOL,
    atol=None,
    max_steps=100_000,
):
    """Peak supersaturation and activated number for a batch of updraft speeds.

    Parameters
    ----------
    y0 : array, shape ``(7 + nr,)``
        Equilibrated initial state (shared across the ensemble).
    r_drys, Nis, kappas : array, shape ``(nr,)``
        Aerosol dry radii (m), number concentrations (m^-3), hygroscopicities.
    accom : float
        Accommodation coefficient.
    V_samples : array, shape ``(n,)``
        Updraft speeds (m/s).
    t_end : float
        Upper bound on integration time; must exceed the slowest member's time-to-peak
        (use ``z_cap / min(V)`` for a target height ``z_cap``).

    Returns
    -------
    dict
        ``S_max`` (n,), ``N_act`` (n, m^-3), ``T_smax`` (n,), and ``activated`` (n, bool)
        -- the last flags members that reached a genuine interior maximum before
        ``t_end``. Activated number uses the equilibrium criterion (modal critical
        supersaturation below ``S_max``) evaluated at the peak temperature.
    """
    y0 = jnp.asarray(y0)
    r_drys = jnp.asarray(r_drys)
    Nis = jnp.asarray(Nis)
    kappas = jnp.asarray(kappas)
    nr = int(y0.shape[0] - c.N_STATE_VARS)
    if atol is None:
        atol = atol_vector(nr)
    V_samples = jnp.asarray(V_samples)

    def _member(V):
        args = (r_drys, Nis, kappas, accom, ConstantV(V))
        sol = _solve_to_smax(y0, args, t_end, rtol, atol, max_steps)
        y_peak = sol.ys[-1]
        s_max = y_peak[6]
        T_peak = y_peak[2]
        t_peak = sol.ts[-1]
        _, s_crit = kohler_crit_approx(T_peak, r_drys, kappas)
        n_act = jnp.sum(jnp.where(s_crit < s_max, Nis, 0.0))
        return s_max, n_act, T_peak, t_peak < t_end

    s_max, n_act, T_smax, activated = jax.jit(jax.vmap(_member))(V_samples)
    return {
        "S_max": np.asarray(s_max),
        "N_act": np.asarray(n_act),
        "T_smax": np.asarray(T_smax),
        "activated": np.asarray(activated),
        "V": np.asarray(V_samples),
    }


def run_updraft_ensemble(
    aerosols,
    T0,
    S0,
    P0,
    *,
    mean,
    std,
    n,
    seed=0,
    v_min=1e-2,
    accom=c.ac,
    z_cap=400.0,
    t_end=None,
    max_steps=100_000,
):
    """Convenience: sample a Gaussian updraft distribution and run the ensemble.

    Equilibrates ``y0`` once for the given aerosols/conditions, samples ``n`` updraft
    speeds from ``Normal(mean, std)`` (clipped at ``v_min``), and returns the per-member
    ``S_max``/``N_act`` (see :func:`smax_nact_ensemble`). ``t_end`` defaults to
    ``z_cap / min(V)`` so even the slowest member reaches its supersaturation maximum.
    """
    species, r_drys, kappas, Nis = [], [], [], []
    for aer in aerosols:
        r_drys.extend(aer.r_drys)
        kappas.extend([aer.kappa] * aer.nr)
        Nis.extend(aer.Nis)
        species.extend([aer.species] * aer.nr)
    r_drys = np.asarray(r_drys)
    kappas = np.asarray(kappas)
    Nis = np.asarray(Nis)

    y0 = np.asarray(equilibrate_initial_state(T0, S0, P0, r_drys, kappas, Nis))
    V_samples = sample_gaussian_updrafts(mean, std, n, seed=seed, v_min=v_min)
    if t_end is None:
        t_end = float(z_cap / np.min(V_samples))

    result = smax_nact_ensemble(
        y0, r_drys, Nis, kappas, float(accom), V_samples, t_end, max_steps=max_steps
    )
    result["t_end"] = t_end
    return result
