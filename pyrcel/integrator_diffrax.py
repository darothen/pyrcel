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
from .updraft import ConstantV  # noqa: E402

#: CVode-equivalent tolerances (see pyrcel/integrator.py).
STATE_RTOL = 1e-7
STATE_ATOL = [1e-4, 1e-4, 1e-4, 1e-10, 1e-10, 1e-4, 1e-8]
RADIUS_ATOL = 1e-12

_TERM = dfx.ODETerm(parcel_ode_sys)


def atol_vector(nr: int) -> jnp.ndarray:
    """Build the per-component absolute tolerance vector for the ODE solver.

    Parameters
    ----------
    nr : int
        Number of aerosol radius bins.

    Returns
    -------
    jnp.ndarray, shape ``(7 + nr,)``
        Absolute tolerance vector: CVode-equivalent values for bulk state
        variables followed by ``RADIUS_ATOL`` for each radius bin.
    """
    return jnp.asarray(STATE_ATOL + [RADIUS_ATOL] * int(nr))


@functools.partial(jax.jit, static_argnames=("max_steps", "progress_meter"))
def _solve(y0, args, ts, rtol, atol, dtmax, max_steps, progress_meter=dfx.NoProgressMeter()):
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
        progress_meter=progress_meter,
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
    progress_meter=None,
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
    progress_meter : diffrax.AbstractProgressMeter, optional
        Live progress reporting (e.g. ``diffrax.TextProgressMeter()``) for interactive
        runs. Defaults to ``NoProgressMeter`` -- keep it that way under ``vmap``/large
        batches and in the differentiable core to avoid per-element host syncs.

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
    if progress_meter is None:
        progress_meter = dfx.NoProgressMeter()
    return _solve(y0, args, ts, rtol, atol, dtmax, max_steps, progress_meter)


def integrate_parcel_arrays(y0, args, ts, **kwargs):
    """Integrate the parcel ODE and return raw arrays.

    Convenience wrapper around :func:`integrate_parcel` that unpacks the
    ``diffrax`` solution object into plain arrays.

    Parameters
    ----------
    y0 : array, shape ``(7 + nr,)``
        Initial (equilibrated) state vector.
    args : tuple
        ``(r_drys, Nis, kappas, accom, V)`` for :func:`parcel_ode_sys`.
    ts : array, shape ``(n_out,)``
        Output times (monotonic); see :func:`integrate_parcel`.
    **kwargs
        Forwarded to :func:`integrate_parcel`.

    Returns
    -------
    ts : array
        Output times.
    ys : array, shape ``(n_out, 7 + nr)``
        State trajectory.
    success : bool
        Whether the solve completed successfully.
    """
    sol = integrate_parcel(y0, args, ts, **kwargs)
    success = bool(sol.result == dfx.RESULTS.successful)
    return sol.ts, sol.ys, success


# --- Differentiable diagnostics (design §5.3, §7.6) ------------------------------

def max_supersaturation(y0, args, ts, *, rtol=STATE_RTOL, atol=None, max_steps=100_000):
    """Peak supersaturation over the output grid, as a differentiable scalar.

    A fixed-horizon solve (no data-dependent event) so the whole computation is a
    pure function compatible with ``jax.grad`` / ``jax.jacfwd`` / ``jax.vmap``.
    ``diffrax``'s default ``RecursiveCheckpointAdjoint`` provides reverse-mode
    gradients through the solve; ``ts`` must span the supersaturation peak.

    Parameters
    ----------
    y0 : array, shape ``(7 + nr,)``
        Initial (equilibrated) state vector.
    args : tuple
        ``(r_drys, Nis, kappas, accom, V)`` for :func:`parcel_ode_sys`.
    ts : array, shape ``(n_out,)``
        Output times (monotonic). Must span the supersaturation peak.
    rtol, atol : float or array, optional
        Solver tolerances.
    max_steps : int, optional
        Upper bound on adaptive steps.

    Returns
    -------
    float
        Maximum supersaturation ``S_max`` over ``ts``.
    """
    y0 = jnp.asarray(y0)
    ts = jnp.asarray(ts)
    nr = int(y0.shape[0] - N_STATE_VARS)
    if atol is None:
        atol = atol_vector(nr)
    sol = _solve(y0, args, ts, rtol, atol, None, max_steps)
    return jnp.max(sol.ys[:, 6])


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


@functools.partial(jax.jit, static_argnames=("max_steps",))
def _solve_to_depth(y0, args, t_end, z_target, rtol, atol, max_steps):
    # Stop when altitude z (= y[0]) first rises through z_target. Works for any updraft
    # (constant or V(t)); z is strictly increasing since dz/dt = V > 0.
    def cond_fn(t, y, args, **kwargs):
        return y[0] - z_target

    event = dfx.Event(
        cond_fn, root_finder=optx.Newton(rtol=1e-8, atol=1e-10), direction=True
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

    Unlike sampling a saved trajectory, ``t_smax`` is root-found rather than
    quantized to the output cadence.

    Parameters
    ----------
    y0 : array, shape ``(7 + nr,)``
        Initial (equilibrated) state vector.
    args : tuple
        ``(r_drys, Nis, kappas, accom, V)`` for :func:`parcel_ode_sys`.
    t_end : float
        Upper bound on integration time, s.
    rtol, atol : float or array, optional
        Solver tolerances.
    max_steps : int, optional
        Upper bound on adaptive steps.

    Returns
    -------
    t_smax : float
        Time of the supersaturation maximum, s.
    smax : float
        Peak supersaturation value.
    y_smax : array, shape ``(7 + nr,)``
        State vector at ``t_smax``.
    activated : bool
        ``False`` if no downward zero-crossing of ``dS/dt`` occurred before
        ``t_end`` (the parcel never reached a maximum in the integration window).
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


def terminate_cutoff_time(
    y0,
    args,
    t_end,
    *,
    terminate_depth: float,
    t_smax,
    y_smax,
    activated: bool,
    rtol: float = STATE_RTOL,
    atol=None,
    max_steps: int = 100_000,
) -> float:
    """Compute the integration stop time for altitude-based cutoff past ``S_max``.

    Parameters
    ----------
    y0 : array, shape ``(7 + nr,)``
        Initial state vector.
    args : tuple
        ``(r_drys, Nis, kappas, accom, V)`` for :func:`parcel_ode_sys`.
    t_end : float
        Upper bound on integration time, s.
    terminate_depth : float
        Extra vertical distance (m) past ``z_smax`` before stopping.
    t_smax : float
        Time of the supersaturation maximum, s.
    y_smax : array
        State vector at ``t_smax``.
    activated : bool
        If ``False``, returns ``t_end`` without any further integration.
    rtol, atol : float or array, optional
        Solver tolerances.
    max_steps : int, optional
        Upper bound on adaptive steps.

    Returns
    -------
    float
        Cutoff time, s. Always ``<= t_end``.
    """
    if not activated:
        return float(t_end)
    t_smax_f = float(t_smax)
    V = args[4]
    if isinstance(V, ConstantV):
        return min(t_smax_f + terminate_depth / float(V.V), float(t_end))
    if not callable(V):
        return min(t_smax_f + terminate_depth / float(V), float(t_end))
    z_target = float(y_smax[0]) + terminate_depth
    nr = int(jnp.asarray(y0).shape[0] - N_STATE_VARS)
    if atol is None:
        atol = atol_vector(nr)
    sol_d = _solve_to_depth(y0, args, t_end, z_target, rtol, atol, max_steps)
    return min(float(sol_d.ts[-1]), float(t_end))


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
    progress_meter=None,
    phase_timer=None,
):
    """Integrate, stopping ``terminate_depth`` metres past the supersaturation max.

    Reproduces the ``master`` ``terminate=True`` semantics: locate ``S_max`` (via the
    ``dS/dt`` event), then continue an extra ``terminate_depth`` metres
    (``= terminate_depth / V`` seconds) and stop. Output is produced at ``output_dt``
    cadence over ``[0, t_cutoff]`` from a single adaptive solve.

    The cutoff is altitude-based (stop ``terminate_depth`` metres above ``z_smax``), so
    it is correct for both a constant updraft and a time-varying ``V(t)``: for constant
    ``V`` the cutoff time is computed analytically, otherwise it is root-found with a
    second ``z``-crossing event.

    This is the interactive / parity path: ``t_cutoff`` is data-dependent so the output
    grid length is dynamic (the outer call runs eagerly; the solves are jitted). For
    ``jit``/``vmap``/``grad`` use the fixed-horizon :func:`integrate_parcel` /
    :func:`find_smax` / :func:`max_supersaturation` instead.

    Returns ``(ts, ys, info)`` with ``ts``/``ys`` as numpy arrays and ``info`` a dict
    of ``t_smax``, ``smax``, ``t_cutoff``, ``activated``, ``success``.
    """
    y0 = jnp.asarray(y0)
    nr = int(y0.shape[0] - N_STATE_VARS)
    if atol is None:
        atol = atol_vector(nr)

    def _find():
        return find_smax(y0, args, t_end, rtol=rtol, atol=atol, max_steps=max_steps)

    if phase_timer is not None:
        (t_smax, smax, y_smax, activated), _ = phase_timer.run(
            ("smax", nr), "S_max event solve", _find
        )
    else:
        t_smax, smax, y_smax, activated = _find()
    t_smax_f = float(t_smax)
    if not activated:
        t_cutoff = float(t_end)
    else:
        t_cutoff = terminate_cutoff_time(
            y0,
            args,
            t_end,
            terminate_depth=terminate_depth,
            t_smax=t_smax,
            y_smax=y_smax,
            activated=True,
            rtol=rtol,
            atol=atol,
            max_steps=max_steps,
        )

    ts = np.append(np.arange(0.0, t_cutoff, output_dt), t_cutoff)
    ts = jnp.asarray(ts)
    if progress_meter is None:
        progress_meter = dfx.NoProgressMeter()

    def _traj():
        return _solve(y0, args, ts, rtol, atol, None, max_steps, progress_meter)

    if phase_timer is not None:
        sol, _ = phase_timer.run(("traj", nr, len(ts)), "trajectory solve", _traj)
    else:
        sol = _traj()
    info = {
        "t_smax": t_smax_f,
        "smax": float(smax),
        "t_cutoff": t_cutoff,
        "activated": activated,
        "success": bool(sol.result == dfx.RESULTS.successful),
        "z_smax": float(y_smax[0]),
        "z_end": float(np.asarray(sol.ys)[-1, 0]),
    }
    return np.asarray(sol.ts), np.asarray(sol.ys), info


def integrate_parcel_chunked(
    y0,
    args,
    t_final,
    output_dt,
    *,
    chunk_dt=10.0,
    rtol: float = STATE_RTOL,
    atol=None,
    max_steps: int = 100_000,
    on_step=None,
):
    """Interactive-only integration in ``chunk_dt`` slices with optional per-chunk callback.

    Unlike :func:`integrate_parcel`, this runs a Python loop over successive
    ``diffeqsolve`` calls so host code can print live diagnostics between chunks. It is
    **not** ``jit``/``vmap``/``grad``-safe and must not be used on the differentiable
    path. Adaptive solver state is re-initialised at each chunk boundary (acceptable for
    CLI feedback; see design doc §4.7).

    Parameters
    ----------
    on_step : callable, optional
        ``on_step(step, t, z, T, S_pct, wall_total, wall_delta)`` invoked immediately
        before each chunk integrate (master ``CVODEIntegrator`` semantics).

    Returns
    -------
    ts, ys, success
        Concatenated output grid and whether every chunk reported success.
    """
    import time

    y0 = jnp.asarray(y0)
    t_final = float(t_final)
    output_dt = float(output_dt)
    chunk_dt = float(chunk_dt)
    nr = int(y0.shape[0] - N_STATE_VARS)
    if atol is None:
        atol = atol_vector(nr)

    y_curr = y0
    t_curr = 0.0
    ts_parts: list[np.ndarray] = []
    ys_parts: list[np.ndarray] = []
    step = 0
    wall_total = 0.0
    wall_prev = time.perf_counter()
    success = True

    while t_curr < t_final - 1e-12:
        step += 1
        t_next = min(t_curr + chunk_dt, t_final)
        n_out = max(1, int(round((t_next - t_curr) / output_dt)))
        chunk_ts = np.linspace(t_curr, t_next, n_out + 1)

        if on_step is not None:
            state = np.asarray(y_curr)
            now = time.perf_counter()
            wall_delta = now - wall_prev
            wall_total += wall_delta
            on_step(
                step,
                t_curr,
                float(state[0]),
                float(state[2]),
                float(state[6]) * 100.0,
                wall_total,
                wall_delta,
            )
            wall_prev = now

        sol = integrate_parcel(
            y_curr,
            args,
            jnp.asarray(chunk_ts),
            rtol=rtol,
            atol=atol,
            max_steps=max_steps,
        )
        jax.block_until_ready(sol.ys)
        if not bool(sol.result == dfx.RESULTS.successful):
            success = False
            break

        ts = np.asarray(sol.ts)
        ys = np.asarray(sol.ys)
        if ts_parts:
            ts_parts.append(ts[1:])
            ys_parts.append(ys[1:])
        else:
            ts_parts.append(ts)
            ys_parts.append(ys)

        y_curr = sol.ys[-1]
        t_curr = float(ts[-1])

    if not ts_parts:
        return np.array([0.0]), np.asarray(y0)[None, :], False
    return np.concatenate(ts_parts), np.concatenate(ys_parts, axis=0), success
