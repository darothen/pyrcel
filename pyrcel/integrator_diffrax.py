"""diffrax integrator for the v2 parcel model (design doc §4.4, §4.7).

A single adaptive ``diffeqsolve`` with a stiff ESDIRK solver (``Kvaerno5``) and a
``PIDController``, replacing the chunked Assimulo/CVode loop. Output at an arbitrary
cadence comes from ``SaveAt(ts=...)`` via dense interpolation -- no manual
``solver_dt`` chunking.

Tolerances mirror the CVode setup in ``pyrcel/integrator.py``
(``rtol=1e-7`` and a per-component ``atol`` vector). ``Kvaerno5`` is A-/L-stable and
stiffly accurate, the diffrax analog of CVode's BDF; its default ``VeryChord`` Newton
root finder inherits these tolerances automatically.

This is the differentiable/batchable numerical core: it is a pure function of
``(y0, args)`` and is ``jit``/``vmap``/``grad``-friendly. Interactive niceties
(progress meter, summary tables) are layered on separately (Phase 6).
"""

from __future__ import annotations

import functools

import jax

jax.config.update("jax_enable_x64", True)

import diffrax as dfx  # noqa: E402
import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402
import optimistix as optx  # noqa: E402

from .parcel_aux_jax import N_STATE_VARS, parcel_ode_sys  # noqa: E402

#: CVode-equivalent tolerances (see pyrcel/integrator.py).
STATE_RTOL = 1e-7
STATE_ATOL = [1e-4, 1e-4, 1e-4, 1e-10, 1e-10, 1e-4, 1e-8]
RADIUS_ATOL = 1e-12

_TERM = dfx.ODETerm(parcel_ode_sys)


def atol_vector(nr: int) -> jnp.ndarray:
    """Per-component absolute tolerance: bulk states + one entry per radius."""
    return jnp.asarray(STATE_ATOL + [RADIUS_ATOL] * int(nr))


@functools.partial(jax.jit, static_argnames=("max_steps",))
def _solve(y0, args, ts, rtol, atol, dtmax, max_steps):
    controller = dfx.PIDController(rtol=rtol, atol=atol, dtmax=dtmax)
    return dfx.diffeqsolve(
        _TERM,
        dfx.Kvaerno5(),
        t0=ts[0],
        t1=ts[-1],
        dt0=None,
        y0=y0,
        args=args,
        stepsize_controller=controller,
        saveat=dfx.SaveAt(ts=ts),
        max_steps=max_steps,
        throw=False,
    )


def integrate_parcel(
    y0,
    args,
    ts,
    *,
    rtol: float = STATE_RTOL,
    atol=None,
    dtmax: float | None = None,
    max_steps: int = 100_000,
):
    """Integrate the parcel ODE and return the full ``diffrax`` solution.

    Parameters
    ----------
    y0 : array, shape ``(7 + nr,)``
        Initial (equilibrated) state vector.
    args : tuple
        ``(r_drys, Nis, kappas, accom, V)`` for :func:`parcel_ode_sys`.
    ts : array, shape ``(n_out,)``
        Output times (monotonic). Output is dense-interpolated at these times from a
        single adaptive solve; ``ts[0]`` is ``t0`` and ``ts[-1]`` is ``t1``.
    rtol, atol : float / array
        Solver tolerances. ``atol`` defaults to the per-component CVode vector.
    dtmax : float, optional
        Maximum internal step. ``None`` lets the controller choose freely.
    max_steps : int
        Upper bound on adaptive steps (must be finite under ``jit``).

    Returns
    -------
    diffrax.Solution
        ``sol.ts``, ``sol.ys`` (shape ``(n_out, 7 + nr)``), ``sol.result``.
    """
    y0 = jnp.asarray(y0)
    ts = jnp.asarray(ts)
    nr = int(y0.shape[0] - N_STATE_VARS)
    if atol is None:
        atol = atol_vector(nr)
    return _solve(y0, args, ts, rtol, atol, dtmax, max_steps)


def integrate_parcel_arrays(y0, args, ts, **kwargs):
    """Convenience wrapper returning ``(ts, ys, success)`` as plain arrays."""
    sol = integrate_parcel(y0, args, ts, **kwargs)
    success = bool(sol.result == dfx.RESULTS.successful)
    return sol.ts, sol.ys, success


# --- Event-based S_max termination (design §4.5 option A) -------------------------

def _dS_dt(t, y, args, **kwargs):
    """Continuous event condition: the supersaturation tendency dS/dt."""
    return parcel_ode_sys(t, y, args)[6]


@functools.partial(jax.jit, static_argnames=("max_steps",))
def _solve_to_smax(y0, args, t_end, rtol, atol, max_steps):
    # ``direction=False`` -> trigger only on a downward zero-crossing of dS/dt
    # (the supersaturation maximum), not the (numerically possible) initial rise.
    event = dfx.Event(
        _dS_dt, root_finder=optx.Newton(rtol=1e-8, atol=1e-12), direction=False
    )
    controller = dfx.PIDController(rtol=rtol, atol=atol)
    return dfx.diffeqsolve(
        _TERM,
        dfx.Kvaerno5(),
        t0=0.0,
        t1=t_end,
        dt0=None,
        y0=y0,
        args=args,
        stepsize_controller=controller,
        saveat=dfx.SaveAt(t1=True),
        event=event,
        max_steps=max_steps,
        throw=False,
    )


def find_smax(y0, args, t_end, *, rtol=STATE_RTOL, atol=None, max_steps=100_000):
    """Precisely localize the supersaturation maximum via a ``dS/dt = 0`` event.

    Returns ``(t_smax, smax, y_smax, activated)`` where ``activated`` is False if no
    downward crossing occurred before ``t_end`` (the parcel never reached a maximum).
    Unlike sampling a saved trajectory, ``t_smax`` here is root-found, not quantized to
    the output cadence.
    """
    y0 = jnp.asarray(y0)
    nr = int(y0.shape[0] - N_STATE_VARS)
    if atol is None:
        atol = atol_vector(nr)
    sol = _solve_to_smax(y0, args, t_end, rtol, atol, max_steps)
    t_smax = sol.ts[-1]
    y_smax = sol.ys[-1]
    smax = y_smax[6]
    activated = bool(t_smax < t_end)
    return t_smax, smax, y_smax, activated


def integrate_parcel_terminated(
    y0,
    args,
    t_end,
    output_dt,
    *,
    terminate_depth: float = 100.0,
    rtol: float = STATE_RTOL,
    atol=None,
    max_steps: int = 100_000,
):
    """Integrate, stopping ``terminate_depth`` metres past the supersaturation max.

    Reproduces the ``master`` ``terminate=True`` semantics: locate ``S_max`` (via the
    ``dS/dt`` event), then continue an extra ``terminate_depth`` metres
    (``= terminate_depth / V`` seconds) and stop. Output is produced at ``output_dt``
    cadence over ``[0, t_cutoff]`` from a single adaptive solve.

    This is the interactive / parity path: ``t_cutoff`` is data-dependent so the output
    grid length is dynamic (the outer call runs eagerly; both solves are jitted). For
    ``jit``/``vmap``/``grad`` use the fixed-horizon :func:`integrate_parcel` /
    :func:`find_smax` instead.

    Returns ``(ts, ys, info)`` with ``ts``/``ys`` as numpy arrays and ``info`` a dict
    of ``t_smax``, ``smax``, ``t_cutoff``, ``activated``, ``success``.
    """
    y0 = jnp.asarray(y0)
    nr = int(y0.shape[0] - N_STATE_VARS)
    if atol is None:
        atol = atol_vector(nr)

    V = args[4]
    if callable(V):
        raise NotImplementedError("time-varying V(t) termination is deferred to Phase 5")
    V = float(V)

    t_smax, smax, _, activated = find_smax(
        y0, args, t_end, rtol=rtol, atol=atol, max_steps=max_steps
    )
    t_smax_f = float(t_smax)
    if activated:
        t_cutoff = min(t_smax_f + terminate_depth / V, float(t_end))
    else:
        t_cutoff = float(t_end)

    ts = np.append(np.arange(0.0, t_cutoff, output_dt), t_cutoff)
    ts = jnp.asarray(ts)
    sol = _solve(y0, args, ts, rtol, atol, None, max_steps)
    info = {
        "t_smax": t_smax_f,
        "smax": float(smax),
        "t_cutoff": t_cutoff,
        "activated": activated,
        "success": bool(sol.result == dfx.RESULTS.successful),
    }
    return np.asarray(sol.ts), np.asarray(sol.ys), info
