"""High-level v2 parcel-model interface (design doc §4.7 mode 2, Phase 6).

A thin, user-facing wrapper over the differentiable JAX/diffrax core that mirrors the
spirit of the master :class:`pyrcel.parcel.ParcelModel` without numba or Assimulo. It
ties together the three v2 pieces:

* equilibration of the initial wet radii (:mod:`pyrcel.equilibrate_jax`),
* a single adaptive ``diffrax`` solve, with optional event-based ``S_max`` termination
  (:mod:`pyrcel.integrator_diffrax`), and
* post-solve activation diagnostics (:func:`pyrcel.activation.binned_activation`).

There are deliberately *two* presentation modes over the same numerical core (design
§4.7): the functions in :mod:`pyrcel.integrator_diffrax` are the ``jit``/``vmap``/``grad``
core, while this class is the **interactive** layer -- an optional live progress meter and
a post-solve summary table -- and is intentionally kept out of the differentiable path.

``ParcelModel`` itself is a plain mutable Python class (per the locked decision); only the
inner vector field / updraft are Equinox modules.
"""

from __future__ import annotations

import numpy as np

from . import constants as c
from .activation import binned_activation
from .equilibrate_jax import equilibrate_initial_state
from .integrator_diffrax import (
    integrate_parcel,
    integrate_parcel_terminated,
)
from .updraft import AbstractUpdraft, ConstantV

__all__ = ["ParcelModelJAX"]


class ParcelModelJAX:
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

    Notes
    -----
    Equilibration runs at construction (like ``master``'s ``_setup_run``), so ``y0`` is
    available immediately as :pyattr:`y0`.
    """

    def __init__(self, aerosols, V, T0, S0, P0, accom=c.ac, console=False):
        self.aerosols = list(aerosols)
        self.V = V
        self.T0 = float(T0)
        self.S0 = float(S0)
        self.P0 = float(P0)
        self.accom = float(accom)
        self.console = console

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

        y0 = equilibrate_initial_state(
            self.T0, self.S0, self.P0, self._r_drys, self._kappas, self._Nis
        )
        self.y0 = np.asarray(y0)

        # Outputs, populated by run().
        self.x = None
        self.time = None
        self.heights = None
        self._summary = None

        if self.console:
            self._print_setup()

    @property
    def args(self) -> tuple:
        """The ``(r_drys, Nis, kappas, accom, V)`` parameter tuple for the integrator."""
        V = self.V if isinstance(self.V, AbstractUpdraft) else float(self.V)
        return (self._r_drys, self._Nis, self._kappas, self.accom, V)

    def run(
        self,
        t_end,
        output_dt=1.0,
        *,
        terminate=True,
        terminate_depth=10.0,
        max_steps=100_000,
        progress=False,
        output_fmt="dataframes",
    ):
        """Run the simulation.

        Parameters
        ----------
        t_end : float
            Maximum integration time (s). With ``terminate=True`` the run usually stops
            earlier, ``terminate_depth`` metres above the supersaturation maximum.
        output_dt : float
            Output cadence (s); the trajectory is dense-interpolated to this grid.
        terminate : bool
            Stop shortly after ``S_max`` (event-based; reproduces ``master``'s behaviour).
        terminate_depth : float
            Extra depth (m) to integrate past ``S_max`` before stopping.
        max_steps : int
            Adaptive-step cap for the solver.
        progress : bool
            Show a live text progress meter during the solve.
        output_fmt : {'dataframes', 'arrays', 'smax'}
            * ``'dataframes'`` -- ``(parcel_df, {species: aerosol_df})`` (pandas),
            * ``'arrays'`` -- ``(state_array, heights)``,
            * ``'smax'`` -- the scalar peak supersaturation.

        Returns
        -------
        Depends on ``output_fmt``; see above.
        """
        if output_fmt not in ("dataframes", "arrays", "smax"):
            raise ValueError(f"invalid output_fmt {output_fmt!r}")

        import diffrax as dfx

        meter = dfx.TextProgressMeter() if progress else None

        peak = None
        if terminate:
            ts, ys, info = integrate_parcel_terminated(
                self.y0, self.args, t_end, output_dt,
                terminate_depth=terminate_depth, max_steps=max_steps,
                progress_meter=meter,
            )
            success = info["success"]
            if info["activated"]:
                # Event-localized peak (precise), not the output-grid argmax.
                peak = (info["smax"], info["t_smax"])
        else:
            ts = np.append(np.arange(0.0, float(t_end), output_dt), float(t_end))
            sol = integrate_parcel(
                self.y0, self.args, ts, max_steps=max_steps, progress_meter=meter
            )
            ys = np.asarray(sol.ys)
            ts = np.asarray(sol.ts)
            success = bool(sol.result == dfx.RESULTS.successful)

        if not success:
            from .util import ParcelModelError

            raise ParcelModelError("diffrax integration failed to complete.")

        self.x = np.asarray(ys)
        self.time = np.asarray(ts)
        self.heights = self.x[:, c.STATE_VAR_MAP["z"]]
        self._summary = self._compute_summary(peak)

        if self.console:
            self._print_summary()

        if output_fmt == "smax":
            return float(self._summary["S_max"])
        if output_fmt == "arrays":
            return self.x, self.heights
        return self._to_dataframes()

    # --- diagnostics --------------------------------------------------------------

    def _compute_summary(self, peak=None) -> dict:
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
        rs = self.x[si, c.N_STATE_VARS:]

        per_species = []
        total_N = 0.0
        total_act = 0.0
        offset = 0
        for aer in self.aerosols:
            nr = aer.nr
            eq, kn, _, _ = binned_activation(S_max, T_smax, rs[offset : offset + nr], aer)
            offset += nr
            N = float(np.sum(aer.Nis))
            per_species.append(
                {"species": aer.species, "eq_act_frac": float(eq),
                 "kn_act_frac": float(kn), "N": N, "N_act": float(eq) * N}
            )
            total_N += N
            total_act += float(eq) * N

        return {
            "S_max": S_max,
            "t_smax": t_smax,
            "T_smax": T_smax,
            "per_species": per_species,
            "total_act_frac": (total_act / total_N) if total_N > 0 else float("nan"),
        }

    def summary(self) -> dict:
        """Return the post-solve summary (``S_max``, ``t_smax``, per-species activation)."""
        if self._summary is None:
            raise RuntimeError("call run() before summary()")
        return self._summary

    # --- output formatting --------------------------------------------------------

    def _to_dataframes(self):
        import pandas as pd

        parcel = pd.DataFrame(
            {var: self.x[:, i] for i, var in enumerate(c.STATE_VARS)},
            index=pd.Index(self.time, name="time"),
        )
        aerosol = {}
        offset = c.N_STATE_VARS
        for aer in self.aerosols:
            nr = aer.nr
            cols = {f"r{j:03d}": self.x[:, offset + j] for j in range(nr)}
            aerosol[aer.species] = pd.DataFrame(
                cols, index=pd.Index(self.time, name="time")
            )
            offset += nr
        return parcel, aerosol

    # --- console output -----------------------------------------------------------

    def _print_setup(self):
        print("Parcel model (JAX/diffrax) -- initial conditions")
        print("-" * 52)
        V_desc = "V(t)" if isinstance(self.V, AbstractUpdraft) and not isinstance(
            self.V, ConstantV
        ) else f"{float(self.V.V) if isinstance(self.V, ConstantV) else self.V:.3g} m/s"
        print(f"  V = {V_desc}   T0 = {self.T0:.2f} K   "
              f"P0 = {self.P0 / 100.0:.1f} hPa   S0 = {self.S0:+.4f}")
        print(f"  aerosol bins: {self._nr}   accom = {self.accom:.3g}")
        z, P, T, wv, wc, wi, S = self.y0[:7]
        print(f"  y0: P={P / 100:.1f} hPa  T={T:.2f} K  "
              f"wv={wv * 1e3:.3g} g/kg  wc={wc * 1e3:.2e} g/kg  S={S:+.4f}")
        print("-" * 52)

    def _print_summary(self):
        s = self._summary
        print("\nSimulation summary")
        print("-" * 52)
        print(f"  S_max = {s['S_max'] * 100:.4f} %  at t = {s['t_smax']:.2f} s "
              f"(T = {s['T_smax']:.2f} K)")
        print(f"  {'species':>12} {'eq_act':>8} {'kn_act':>8} {'N_act':>10} {'N':>10}")
        for p in s["per_species"]:
            print(f"  {p['species']:>12} {p['eq_act_frac']:8.3f} "
                  f"{p['kn_act_frac']:8.3f} {p['N_act']:10.1f} {p['N']:10.1f}")
        print(f"  total activated fraction = {s['total_act_frac']:.3f}")
        print("-" * 52)
