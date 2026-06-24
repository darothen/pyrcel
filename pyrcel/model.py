"""High-level v2 parcel-model interface (design doc §4.7 mode 2, Phase 6).

A thin, user-facing wrapper over the differentiable JAX/diffrax core that mirrors the
spirit of the legacy parcel model without numba or Assimulo. It ties together the three v2 pieces:

* equilibration of the initial wet radii (:mod:`pyrcel.equilibrate`),
* a single adaptive ``diffrax`` solve, with optional event-based ``S_max`` termination
  (:mod:`pyrcel.integrator`), and
* post-solve activation diagnostics (:func:`pyrcel.activation.binned_activation`).

There are deliberately *two* presentation modes over the same numerical core (design
§4.7): the functions in :mod:`pyrcel.integrator` are the ``jit``/``vmap``/``grad``
core, while this class is the **interactive** layer -- an optional live progress meter and
a post-solve summary table -- and is intentionally kept out of the differentiable path.

:class:`ParcelModel` is a plain mutable Python class (per the locked decision); only the
inner vector field / updraft are Equinox modules.
"""

from __future__ import annotations

import contextlib
from pathlib import Path
from typing import Any

import jax
import numpy as np

from . import constants as c
from .activation import binned_activation
from .console_report import (
    LiveStepPrinter,
    PhaseTimer,
    equilibration_residual,
    print_integration_plan,
    print_setup,
    print_summary,
    print_termination_narrative,
    print_trajectory_table,
)
from .equilibrate import equilibrate_initial_state
from .integrator import (
    STATE_RTOL,
    find_smax,
    integrate_parcel,
    integrate_parcel_chunked,
    integrate_parcel_terminated,
    terminate_cutoff_time,
)
from .model_output import ModelOutput
from .updraft import AbstractUpdraft

__all__ = ["ParcelModel"]


def _resolve_device(device: jax.Device | str | None) -> jax.Device | None:
    """Normalise *device* to a ``jax.Device`` or ``None``."""
    if device is None:
        return None
    if isinstance(device, str):
        try:
            return jax.devices(device)[0]
        except (RuntimeError, IndexError) as exc:
            raise ValueError(f"No JAX device of type {device!r} found: {exc}") from exc
    return device


class ParcelModel:
    """Set up and run a parcel-model simulation on the JAX/diffrax backend.

    Parameters
    ----------
    aerosols : sequence of :class:`pyrcel.aerosol.AerosolSpecies`
        The aerosol population in the parcel.
    V : float or :class:`pyrcel.updraft.AbstractUpdraft`
        Updraft speed (m/s). A scalar is a constant updraft; pass a
        :class:`~pyrcel.updraft.InterpolatedUpdraft` for a time-varying ``V(t)``.
    T0, S0, P0 : float
        Initial temperature (K), supersaturation (0.0 == 100% RH), pressure (Pa).
    accom : float, optional
        Condensation/accommodation coefficient (default :data:`pyrcel.constants.ac`).
    console : bool, optional
        Print an initial-conditions and post-solve summary table.
    device : :class:`jax.Device`, str, or None, optional
        JAX device on which to run the integration.  ``None`` (default) uses
        JAX's current default device (typically the first available GPU when a
        CUDA-capable GPU is present, otherwise CPU).  Pass ``"gpu"`` or
        ``"cpu"`` as a shorthand, or an explicit :class:`jax.Device` obtained
        from :func:`jax.devices`.

        Example::

            model = ParcelModel(..., device="gpu")
            model = ParcelModel(..., device=jax.devices("gpu")[1])  # second GPU

    Notes
    -----
    Equilibration runs at construction (like the legacy ``_setup_run``), so ``y0`` is
    available immediately as :attr:`y0`.  All JAX computation — both equilibration and
    integration — is dispatched to ``device``.  Output arrays (:attr:`x`, :attr:`time`)
    are always returned as NumPy arrays on CPU regardless of ``device``.
    """

    def __init__(
        self,
        aerosols: list,
        V: float | AbstractUpdraft,
        T0: float,
        S0: float,
        P0: float,
        accom: float = c.ac,
        console: bool = False,
        device: jax.Device | str | None = None,
    ) -> None:
        self.aerosols = list(aerosols)
        self.V = V
        self.T0 = float(T0)
        self.S0 = float(S0)
        self.P0 = float(P0)
        self.accom = float(accom)
        self.console = console
        self.device: jax.Device | None = _resolve_device(device)

        species, r_drys, kappas, Nis = [], [], [], []
        for aer in self.aerosols:
            r_drys.extend(aer.r_drys)
            kappas.extend([aer.kappa] * aer.nr)
            Nis.extend(aer.Nis)
            species.extend([aer.species] * aer.nr)
        self.species = np.asarray(species)
        self._r_drys = np.asarray(r_drys)
        self._kappas = np.asarray(kappas)
        self._Nis = np.asarray(Nis)
        self._nr = len(r_drys)

        import time

        t0 = time.perf_counter()
        _equil_ctx = (
            jax.default_device(self.device) if self.device is not None else contextlib.nullcontext()
        )
        with _equil_ctx:
            y0 = equilibrate_initial_state(
                self.T0, self.S0, self.P0, self._r_drys, self._kappas, self._Nis
            )
        self._equil_elapsed = time.perf_counter() - t0
        self.y0 = np.asarray(y0)
        self._equil_residual = equilibration_residual(
            self.T0, self.S0, self._r_drys, self._kappas, self.y0[c.N_STATE_VARS :]
        )

        # Outputs, populated by run().
        self.x = None
        self.time = None
        self.heights = None
        self._summary = None
        self._phase_timer = PhaseTimer()
        self._run_info: dict | None = None

        if self.console:
            print_setup(
                V=self.V,
                T0=self.T0,
                S0=self.S0,
                P0=self.P0,
                accom=self.accom,
                nr=self._nr,
                aerosols=self.aerosols,
                y0=self.y0,
                equil_elapsed=self._equil_elapsed,
                equil_residual=self._equil_residual,
            )

    @property
    def args(self) -> tuple:
        """The ``(r_drys, Nis, kappas, accom, V)`` parameter tuple for the integrator."""
        V = self.V if isinstance(self.V, AbstractUpdraft) else float(self.V)
        return (self._r_drys, self._Nis, self._kappas, self.accom, V)

    def run(
        self,
        t_end: float,
        output_dt: float = 1.0,
        *,
        terminate: bool = True,
        terminate_depth: float = 10.0,
        max_steps: int = 100_000,
        progress: bool = False,
        live: bool = False,
        live_chunk_dt: float = 10.0,
        trajectory_table: bool | None = None,
        mode: str = "full",
    ) -> ModelOutput | float:
        """Run the simulation.

        Parameters
        ----------
        t_end : float
            Maximum integration time (s). With ``terminate=True`` the run usually stops
            earlier, ``terminate_depth`` metres above the supersaturation maximum.
        output_dt : float
            Output cadence (s); the trajectory is dense-interpolated to this grid.
        terminate : bool
            Stop shortly after ``S_max`` (event-based; reproduces legacy behaviour).
        terminate_depth : float
            Extra depth (m) to integrate past ``S_max`` before stopping.
        max_steps : int
            Adaptive-step cap for the solver.
        progress : bool
            Show a live text progress meter during the solve (``False`` by default).
            Mutually exclusive with ``live``.
        live : bool
            Print a legacy-style z/T/S table after each integration chunk (``False`` by
            default). Uses an interactive-only Python chunk loop; mutually exclusive
            with ``progress``.
        live_chunk_dt : float
            Simulation-time length of each live-integration chunk (s). Default ``10``.
        trajectory_table : bool or None
            Print a post-hoc trajectory sample after integration. ``None`` (default)
            follows ``console`` (on when ``console=True``).
        mode : {'full', 'smax', 'nd'}
            * ``'full'`` -- :class:`~pyrcel.model_output.ModelOutput` containing the
              full trajectory; call ``.to_pandas()``, ``.to_polars()``,
              ``.to_xarray()``, ``.to_netcdf()``, ``.to_csv()``, or
              ``.to_parquet()`` on the result.
            * ``'smax'`` -- the scalar peak supersaturation (primary differentiable
              path; no trajectory stored).
            * ``'nd'`` -- total activated droplet number concentration (cm⁻³),
              evaluated at the last trajectory step via a hard radius threshold.
              Uses the same integration as ``'full'``; only the return value differs.
              For a differentiable analog see issue #67.

        Returns
        -------
        :class:`~pyrcel.model_output.ModelOutput` or float
            ``ModelOutput`` for ``mode='full'``; ``float`` for ``mode='smax'``.
        """
        if mode not in ("full", "smax", "nd"):
            raise ValueError(f"invalid mode {mode!r}")
        if live and progress:
            raise ValueError("live and progress are mutually exclusive")

        _run_ctx = (
            jax.default_device(self.device) if self.device is not None else contextlib.nullcontext()
        )
        with _run_ctx:
            import diffrax as dfx

            if trajectory_table is None:
                trajectory_table = self.console and not live

            meter = (
                dfx.NoProgressMeter()
                if live
                else (dfx.TextProgressMeter() if progress else dfx.NoProgressMeter())
            )

            if self.console:
                print_integration_plan(
                    t_end=t_end,
                    output_dt=output_dt,
                    terminate=terminate,
                    terminate_depth=terminate_depth,
                    max_steps=max_steps,
                    rtol=STATE_RTOL,
                    progress=progress,
                    live=live,
                    live_chunk_dt=live_chunk_dt,
                )

            peak = None
            run_info = None
            timer = self._phase_timer if self.console and not live else None
            live_printer = LiveStepPrinter() if live else None

            if live:
                t_final = float(t_end)
                activated = True
                t_smax_f = float("nan")
                smax = float("nan")
                y_smax = None

                if terminate:

                    def _find():
                        return find_smax(
                            self.y0,
                            self.args,
                            t_end,
                            max_steps=max_steps,
                        )

                    if timer is not None:
                        (t_smax, smax, y_smax, activated), _ = timer.run(
                            ("smax", self._nr), "S_max event solve", _find
                        )
                    else:
                        t_smax, smax, y_smax, activated = _find()
                    t_smax_f = float(t_smax)
                    if activated:
                        peak = (float(smax), t_smax_f)
                        t_final = terminate_cutoff_time(
                            self.y0,
                            self.args,
                            t_end,
                            terminate_depth=terminate_depth,
                            t_smax=t_smax,
                            y_smax=y_smax,
                            activated=True,
                            max_steps=max_steps,
                        )

                ts, ys, success = integrate_parcel_chunked(
                    self.y0,
                    self.args,
                    t_final,
                    output_dt,
                    chunk_dt=live_chunk_dt,
                    max_steps=max_steps,
                    on_step=live_printer,
                )
                if live_printer is not None:
                    live_printer.finish()
                ts, ys = np.asarray(ts), np.asarray(ys)
                run_info = {
                    "success": success,
                    "activated": activated,
                    "t_smax": t_smax_f,
                    "smax": smax,
                    "t_cutoff": float(ts[-1]) if len(ts) else t_final,
                    "z_smax": float(y_smax[0]) if y_smax is not None else float("nan"),
                    "z_end": float(ys[-1, c.STATE_VAR_MAP["z"]]) if len(ys) else float("nan"),
                }
            elif terminate:
                ts, ys, run_info = integrate_parcel_terminated(
                    self.y0,
                    self.args,
                    t_end,
                    output_dt,
                    terminate_depth=terminate_depth,
                    max_steps=max_steps,
                    progress_meter=meter,
                    phase_timer=timer,
                )
                ts, ys = np.asarray(ts), np.asarray(ys)
                success = run_info["success"]
                if run_info["activated"]:
                    peak = (run_info["smax"], run_info["t_smax"])
            else:
                ts_arr = np.append(np.arange(0.0, float(t_end), output_dt), float(t_end))

                def _integrate():
                    return integrate_parcel(
                        self.y0, self.args, ts_arr, max_steps=max_steps, progress_meter=meter
                    )

                key = ("fixed", self._nr, len(ts_arr))
                if timer is not None:
                    sol, _ = timer.run(key, "integration (fixed horizon)", _integrate)
                else:
                    sol = _integrate()
                ys = np.asarray(sol.ys)
                ts = np.asarray(sol.ts)
                success = bool(sol.result == dfx.RESULTS.successful)
                run_info = {"success": success, "activated": True, "t_cutoff": float(ts[-1])}

            if not success:
                from .util import ParcelModelError

                raise ParcelModelError("diffrax integration failed to complete.")

            self.x = np.asarray(ys)
            self.time = np.asarray(ts)
            self.heights = self.x[:, c.STATE_VAR_MAP["z"]]
            self._run_info = run_info
            self._summary = self._compute_summary(peak)

            if self.console:
                if run_info is not None:
                    print_termination_narrative(run_info)
                if trajectory_table:
                    print_trajectory_table(self.time, self.x)
                print_summary(self._summary)

            if mode == "smax":
                return float(self._summary["S_max"])
            if mode == "nd":
                return float(self._summary["total_Nd"])
            return ModelOutput(
                time=self.time,
                state=self.x,
                aerosols=self.aerosols,
                summary=self._summary,
                V=self.V,
                T0=self.T0,
                S0=self.S0,
                P0=self.P0,
                accom=self.accom,
            )

    # --- diagnostics --------------------------------------------------------------

    def _compute_summary(self, peak: object = None) -> dict:
        assert self.x is not None
        assert self.time is not None
        S = self.x[:, c.STATE_VAR_MAP["S"]]
        if peak is not None:
            # Use the event-localized (precise) S_max/t_smax; take droplet radii from
            # the nearest output sample for the activation diagnostics.
            S_max, t_smax = float(peak[0]), float(peak[1])
            si = int(np.argmin(np.abs(self.time - t_smax)))
        else:
            si = int(np.argmax(S))
            S_max = float(S[si])
            t_smax = float(self.time[si])
        T_smax = float(self.x[si, c.STATE_VAR_MAP["T"]])
        rs_smax = self.x[si, c.N_STATE_VARS :]

        # Nd snapshot: wet radii at the last trajectory step (terminate_depth m past
        # S_max when terminate=True), compared against per-bin critical radii.
        i_nd = len(self.time) - 1
        T_nd = float(self.x[i_nd, c.STATE_VAR_MAP["T"]])
        rs_nd = self.x[i_nd, c.N_STATE_VARS :]
        nd_t_eval = float(self.time[i_nd])

        per_species = []
        total_N = 0.0
        total_act = 0.0
        total_Nd = 0.0
        offset = 0
        for aer in self.aerosols:
            nr = aer.nr
            eq, kn, _, _ = binned_activation(S_max, T_smax, rs_smax[offset : offset + nr], aer)
            # Suppress numpy divide-by-zero warning from alpha=N_kn/N_eq when
            # S_max is below all s_crits (N_eq=0). We only use kn_frac; alpha
            # is discarded. The same edge case exists in the S_max snapshot call above.
            with np.errstate(divide="ignore", invalid="ignore"):
                _, nd_frac, _, _ = binned_activation(S_max, T_nd, rs_nd[offset : offset + nr], aer)
            offset += nr
            N = float(np.sum(aer.Nis))
            Nd_mode = float(nd_frac) * N
            per_species.append(
                {
                    "species": aer.species,
                    "eq_act_frac": float(eq),
                    "kn_act_frac": float(kn),
                    "N": N,
                    "N_act": float(eq) * N,
                    "nd_frac": float(nd_frac),
                    "Nd": Nd_mode,
                }
            )
            total_N += N
            total_act += float(eq) * N
            total_Nd += Nd_mode

        return {
            "S_max": S_max,
            "t_smax": t_smax,
            "T_smax": T_smax,
            "z_smax": float(self.x[si, c.STATE_VAR_MAP["z"]]),
            "per_species": per_species,
            "total_act_frac": (total_act / total_N) if total_N > 0 else float("nan"),
            "total_Nd": total_Nd,
            "total_nd_frac": (total_Nd / total_N) if total_N > 0 else float("nan"),
            "nd_t_eval": nd_t_eval,
        }

    def summary(self) -> dict:
        """Return the post-solve summary (``S_max``, ``t_smax``, per-species activation)."""
        if self._summary is None:
            raise RuntimeError("call run() before summary()")
        return self._summary

    # --- output formatting (delegates to ModelOutput) ----------------------------

    def _make_output(self) -> ModelOutput:
        if self.x is None or self.time is None or self._summary is None:
            raise RuntimeError("call run() before accessing output")
        return ModelOutput(
            time=self.time,
            state=self.x,
            aerosols=self.aerosols,
            summary=self._summary,
            V=self.V,
            T0=self.T0,
            S0=self.S0,
            P0=self.P0,
            accom=self.accom,
        )

    def to_dataset(self) -> Any:
        """Return a CF-flavoured :class:`xarray.Dataset`.

        Delegates to :meth:`~pyrcel.model_output.ModelOutput.to_xarray`.
        """
        return self._make_output().to_xarray()

    def save_netcdf(self, filename: str | Path) -> str | Path:
        """Write the run to a NetCDF file (see :meth:`to_dataset`)."""
        self._make_output().to_netcdf(filename)
        if self.console:
            from .console_report import configure_logging

            configure_logging()
            log = __import__("logging").getLogger("pyrcel")
            log.info("Saved output to %s", filename)
        return filename


# Backward-compatible alias — prefer ParcelModel.
ParcelModelJAX = ParcelModel
