"""Formatted console output for the interactive [ParcelModel][pyrcel.model.ParcelModel] layer.

Default UX (``console=True``): single adaptive solve, compile/integrate phase banners,
rich setup and integration-plan tables, a post-hoc trajectory sample, and the activation
summary. ``live=True`` adds a master-style per-chunk z/T/S table during integration
(interactive-only; kept out of the ``jit``/``grad`` path).

Uses the ``logging`` module for one-line phase messages (``[pyrcel] ...``) and ``print``
for aligned tables so interactive runs stay readable without extra dependencies.
"""

from __future__ import annotations

import logging
import time
from collections.abc import Callable

import numpy as np

from .distributions import Lognorm, MultiModeLognorm
from .updraft import AbstractUpdraft, ConstantV

log = logging.getLogger("pyrcel")

#: Rows shown in the post-hoc trajectory table (subsampled if longer).
MAX_TRAJECTORY_ROWS = 24
#: Phase durations above this likely include one-time XLA compilation.
COMPILE_HINT_SECONDS = 1.0


class _PyrcelFormatter(logging.Formatter):
    """Format log lines without ``%``-style interpolation on the message body."""

    def format(self, record: logging.LogRecord) -> str:
        return f"[pyrcel] {record.getMessage()}"


def configure_logging() -> None:
    """Attach a stdout handler for ``pyrcel`` if none is configured yet."""
    if log.handlers:
        return
    handler = logging.StreamHandler()
    handler.setFormatter(_PyrcelFormatter())
    log.addHandler(handler)
    log.setLevel(logging.INFO)
    log.propagate = False


def backend_banner() -> str:
    """One-line backend / device summary."""
    try:
        import jax

        x64 = bool(jax.config.read("jax_enable_x64"))
        device = jax.devices()[0]
        return f"JAX {jax.__version__} | {device} | float64={'on' if x64 else 'off'}"
    except Exception:
        return "JAX backend (version unavailable)"


class PhaseTimer:
    """Track phase wall times and whether XLA compile was likely involved."""

    def __init__(self) -> None:
        self._seen: set[tuple] = set()

    def run(self, key: tuple, label: str, fn: Callable, *args, **kwargs):
        """Run ``fn``, log duration, and note likely first-time compilation."""
        configure_logging()
        first = key not in self._seen
        if first:
            log.info("Compiling %s...", label)
        else:
            log.info("%s...", label)
        t0 = time.perf_counter()
        result = fn(*args, **kwargs)
        # Block JAX async dispatch so wall time reflects the real work.
        _block_until_ready(result)
        elapsed = time.perf_counter() - t0
        self._seen.add(key)
        if first and elapsed >= COMPILE_HINT_SECONDS:
            log.info("  finished in %.2f s (includes one-time XLA compile)", elapsed)
        else:
            log.info("  finished in %.2f s", elapsed)
        return result, elapsed


def _block_until_ready(result) -> None:
    try:
        if hasattr(result, "block_until_ready"):
            result.block_until_ready()
        elif isinstance(result, tuple):
            for item in result:
                _block_until_ready(item)
        elif isinstance(result, dict):
            for item in result.values():
                _block_until_ready(item)
    except Exception:
        pass


def _hline(width: int = 72, char: str = "-") -> str:
    return char * width


def _v_desc(V) -> str:
    if isinstance(V, ConstantV):
        return f"{float(V.V):.4g} m/s"
    if isinstance(V, AbstractUpdraft):
        return "V(t) profile"
    return f"{float(V):.4g} m/s"


def _aerosol_rows(aerosols) -> list[tuple[str, float, float, int, str, str]]:
    """Build table rows: species, κ, N [cm⁻³], bins, size spec, notes."""
    rows = []
    for aer in aerosols:
        dist = aer.distribution
        if isinstance(dist, Lognorm):
            rows.append(
                (
                    aer.species,
                    aer.kappa,
                    float(aer.total_N),
                    aer.nr,
                    f"μ={dist.mu:g} μm",
                    f"σ={dist.sigma:g}",
                )
            )
        elif isinstance(dist, MultiModeLognorm):
            rows.append(
                (
                    aer.species,
                    aer.kappa,
                    float(aer.total_N),
                    aer.nr,
                    f"{len(dist.mus)} modes",
                    "",
                )
            )
        elif isinstance(dist, dict):
            rows.append((aer.species, aer.kappa, float(aer.total_N), aer.nr, "monodisperse", ""))
        else:
            rows.append((aer.species, aer.kappa, float(aer.total_N), aer.nr, "custom", ""))
    return rows


def print_setup(
    *,
    V,
    T0,
    S0,
    P0,
    accom,
    nr,
    aerosols,
    y0,
    equil_elapsed: float,
    equil_residual: float,
) -> None:
    """Print the parcel configuration, aerosol modes, and equilibrated initial state.

    Parameters
    ----------
    V : float or AbstractUpdraft
        Updraft speed specification.
    T0, S0, P0 : float
        Initial temperature (K), supersaturation, and pressure (Pa).
    accom : float
        Condensation coefficient.
    nr : int
        Total number of aerosol size bins.
    aerosols : list of AerosolSpecies
        Aerosol population.
    y0 : array
        Equilibrated initial state vector.
    equil_elapsed : float
        Wall time for equilibration, s.
    equil_residual : float
        Maximum ``|Seq - S0|`` across all bins after equilibration.
    """
    configure_logging()
    print()
    print("Parcel model (JAX / diffrax)")
    print(_hline())
    print(f"  Backend   {backend_banner()}")
    print()
    print("  Configuration")
    print(f"    {'V':<12} {_v_desc(V):<18} {'accom':<12} {accom:.4g}")
    print(f"    {'T0':<12} {T0:.2f} K{'':<10} {'P0':<12} {P0 / 100.0:.1f} hPa")
    print(f"    {'S0':<12} {S0:+.5f} ({100 * S0:+.2f} %){'':<3} {'bins':<12} {nr}")
    print()
    rows = _aerosol_rows(aerosols)
    print(f"  {'species':<12} {'κ':>6} {'N [cm⁻³]':>12} {'bins':>6}  size")
    print(f"  {_hline(68, '-')}")
    for sp, kap, n, nb, spec, note in rows:
        extra = f"  {note}" if note else ""
        print(f"  {sp:<12} {kap:6.3f} {n:12.1f} {nb:6d}  {spec}{extra}")
    print()
    z, P, T, wv, wc, wi, S = y0[:7]
    print("  Equilibrated initial state")
    print(f"    equilibration   {equil_elapsed * 1e3:.1f} ms  |Seq-S0|_max = {equil_residual:.2e}")
    print(
        f"    {'P':>4} {P / 100.0:8.1f} hPa   {'T':>3} {T:7.2f} K   "
        f"{'wv':>3} {wv * 1e3:8.3g} g/kg   {'wc':>3} {wc * 1e3:8.2e} g/kg   "
        f"{'S':>3} {S:+.5f}"
    )
    print(_hline())


def print_integration_plan(
    *,
    t_end,
    output_dt,
    terminate,
    terminate_depth,
    max_steps,
    rtol,
    progress,
    live=False,
    live_chunk_dt=10.0,
) -> None:
    """Print a summary of the planned integration settings.

    Parameters
    ----------
    t_end : float
        Maximum integration time cap, s.
    output_dt : float
        Output cadence, s.
    terminate : bool
        Whether to stop past ``S_max``.
    terminate_depth : float
        Extra depth (m) to integrate past ``S_max``.
    max_steps : int
        Adaptive-step cap for the solver.
    rtol : float
        Relative tolerance for the ODE solver.
    progress : bool
        Whether a text progress meter is active.
    live : bool, optional
        Whether the live chunk-loop mode is active.
    live_chunk_dt : float, optional
        Chunk size for live mode, s.
    """
    print()
    print("  Integration plan")
    term = f"yes (+{terminate_depth:g} m past S_max)" if terminate else "no (integrate to t_end)"
    if live:
        prog = f"live chunk loop ({live_chunk_dt:g} s chunks)"
    elif progress:
        prog = "TextProgressMeter"
    else:
        prog = "none"
    print(f"    t_end cap     {t_end:g} s")
    print(f"    output_dt     {output_dt:g} s")
    print(f"    terminate     {term}")
    print(f"    solver        Kvaerno5 + PIDController  rtol={rtol:.0e}  max_steps={max_steps}")
    print(f"    progress      {prog}")
    print(_hline())


class LiveStepPrinter:
    """Master-style integration loop table (``CVODEIntegrator.integrate``)."""

    _HEADER = (
        "\n  Integration loop\n\n"
        "    step     time  walltime  Δwalltime |     z       T       S\n"
        "   ------------------------------------|----------------------"
    )
    _ROW = "   {:5d} {:7.2f}s  {:7.2f}s  {:8.2f}s | {:5.1f} {:7.2f} {:6.2f}%"
    _FOOTER = "   ---- end of integration loop ----"

    def __init__(self) -> None:
        self._started = False

    def __call__(
        self,
        step: int,
        t: float,
        z: float,
        T: float,
        S_pct: float,
        wall_total: float,
        wall_delta: float,
    ) -> None:
        if not self._started:
            print(self._HEADER)
            self._started = True
        print(self._ROW.format(step, t, wall_total, wall_delta, z, T, S_pct))

    def finish(self) -> None:
        if self._started:
            print(self._FOOTER)
            print()


def print_termination_narrative(info: dict) -> None:
    """Log a one-line termination summary after the solve completes.

    Parameters
    ----------
    info : dict
        Run-info dict returned by
        [integrate_parcel_terminated][pyrcel.integrator.integrate_parcel_terminated]
        or the ``live`` integration path; expected keys ``activated``, ``smax``,
        ``t_smax``, ``t_cutoff``, and optionally ``z_smax`` / ``z_end``.
    """
    if not info.get("activated", True):
        log.warning("No interior supersaturation maximum before t_end — ran to horizon.")
        return
    smax_pct = 100.0 * info["smax"]
    z_smax = info.get("z_smax")
    z_end = info.get("z_end")
    parts = [
        f"S_max = {smax_pct:.4f} % at t = {info['t_smax']:.2f} s",
    ]
    if z_smax is not None:
        parts[0] += f" (z = {z_smax:.1f} m)"
    parts.append(f"stopped at t = {info['t_cutoff']:.2f} s")
    if z_end is not None:
        parts[-1] += f" (z = {z_end:.1f} m)"
    msg = f"Termination: {parts[0]}; {parts[1]}."
    log.info(msg)


def print_trajectory_table(time, x, *, max_rows: int = MAX_TRAJECTORY_ROWS) -> None:
    """Post-hoc sample of z, T, S, wc along the saved trajectory."""
    time = np.asarray(time)
    x = np.asarray(x)
    n = len(time)
    if n == 0:
        return
    if n > max_rows:
        idx = np.unique(np.linspace(0, n - 1, max_rows, dtype=int))
    else:
        idx = np.arange(n)
    print()
    print("  Trajectory (post-hoc sample)")
    print(f"  {'t [s]':>8} {'z [m]':>8} {'T [K]':>8} {'S [%]':>9} {'wc [g/kg]':>11}")
    print(f"  {_hline(48, '-')}")
    for i in idx:
        row = x[i]
        print(
            f"  {time[i]:8.2f} {row[0]:8.1f} {row[2]:8.2f} {100 * row[6]:9.4f} {row[4] * 1e3:11.3e}"
        )
    print(_hline())


def print_summary(summary: dict) -> None:
    """Print the post-solve activation summary table.

    Parameters
    ----------
    summary : dict
        Summary dict as returned by
        [_compute_summary][pyrcel.model.ParcelModel._compute_summary]; expected keys
        ``S_max``, ``t_smax``, ``T_smax``, ``z_smax``, ``per_species``, and
        ``total_act_frac``.
    """
    s = summary
    print()
    print("  Simulation summary")
    print(_hline())
    print(
        f"  S_max = {s['S_max'] * 100:.4f} %  at t = {s['t_smax']:.2f} s "
        f"(T = {s['T_smax']:.2f} K, z = {s.get('z_smax', float('nan')):.1f} m)"
    )
    print(f"  {'species':>12} {'eq_act':>8} {'kn_act':>8} {'N_act':>10} {'N':>10}")
    for p in s["per_species"]:
        print(
            f"  {p['species']:>12} {p['eq_act_frac']:8.3f} "
            f"{p['kn_act_frac']:8.3f} {p['N_act']:10.1f} {p['N']:10.1f}"
        )
    print(f"  total activated fraction = {s['total_act_frac']:.3f}")
    print(_hline())
    print()


def equilibration_residual(T0, S0, r_drys, kappas, r0s) -> float:
    """Maximum absolute equilibration residual ``|Seq - S0|`` across all bins.

    Parameters
    ----------
    T0 : float
        Initial temperature, K.
    S0 : float
        Initial supersaturation.
    r_drys : array
        Dry radii, m.
    kappas : array
        Hygroscopicities.
    r0s : array
        Equilibrated wet radii, m.

    Returns
    -------
    float
        Maximum ``|Seq(r0, r_dry, T0, kappa) - S0|`` across all bins.
    """
    from .thermo import Seq

    r_drys = np.asarray(r_drys)
    kappas = np.asarray(kappas)
    r0s = np.asarray(r0s)
    if r0s.size == 0:
        return 0.0
    seq = np.array([float(Seq(r, rd, T0, k)) for r, rd, k in zip(r0s, r_drys, kappas)])
    return float(np.max(np.abs(seq - S0)))
