"""Structured output container for a completed :class:`~pyrcel.model.ParcelModel` run.

:class:`ModelOutput` is a plain Python dataclass (not a JAX pytree) that wraps the
raw numpy arrays produced by the integrator and exposes them through a set of
format-conversion methods:

* :meth:`to_pandas` — ``(parcel_df, {species: aerosol_df})`` pandas DataFrames
* :meth:`to_polars` — same structure in polars
* :meth:`to_xarray` — ``xr.Dataset`` with CF-flavoured coordinates and metadata
* :meth:`to_netcdf` — write the xarray Dataset to a NetCDF4 file
* :meth:`to_csv` — write the flat parcel trajectory as CSV
* :meth:`to_parquet` — write the flat parcel trajectory as Parquet
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from . import constants as c
from .thermo import rho_air
from .updraft import AbstractUpdraft, ConstantV

if TYPE_CHECKING:
    import pandas as pd
    import polars as pl
    import xarray as xr


@dataclass
class ModelOutput:
    """Output from a single :class:`~pyrcel.model.ParcelModel` run.

    Attributes
    ----------
    time : np.ndarray, shape ``(n_time,)``
        Simulation time (s).
    state : np.ndarray, shape ``(n_time, 7 + nr)``
        Full state trajectory.  The first seven columns are the bulk parcel
        variables (see :data:`pyrcel.constants.STATE_VARS`); the remaining
        columns are per-bin wet radii (m) ordered as in ``aerosols``.
    aerosols : list of :class:`~pyrcel.aerosol.AerosolSpecies`
        Aerosol modes, in the same order as the radius columns in ``state``.
    summary : dict
        Post-solve diagnostics: ``S_max``, ``t_smax``, ``T_smax``, ``z_smax``,
        ``per_species`` (list of per-mode dicts), ``total_act_frac``.
    V : float or :class:`~pyrcel.updraft.AbstractUpdraft`
        Updraft used for the run (stored for dataset metadata).
    T0, S0, P0, accom : float
        Initial conditions (stored for dataset metadata).
    """

    time: np.ndarray
    state: np.ndarray
    aerosols: list
    summary: dict
    V: float | AbstractUpdraft
    T0: float
    S0: float
    P0: float
    accom: float

    # ------------------------------------------------------------------
    # Convenience accessors for the bulk state variables
    # ------------------------------------------------------------------

    @property
    def heights(self) -> np.ndarray:
        return self.state[:, c.STATE_VAR_MAP["z"]]

    @property
    def T(self) -> np.ndarray:
        return self.state[:, c.STATE_VAR_MAP["T"]]

    @property
    def P(self) -> np.ndarray:
        return self.state[:, c.STATE_VAR_MAP["P"]]

    @property
    def S(self) -> np.ndarray:
        return self.state[:, c.STATE_VAR_MAP["S"]]

    @property
    def wv(self) -> np.ndarray:
        return self.state[:, c.STATE_VAR_MAP["wv"]]

    @property
    def wc(self) -> np.ndarray:
        return self.state[:, c.STATE_VAR_MAP["wc"]]

    @property
    def wi(self) -> np.ndarray:
        return self.state[:, c.STATE_VAR_MAP["wi"]]

    @property
    def Nd(self) -> float:
        """Total activated droplet number concentration (cm⁻³) at the last trajectory step.

        Evaluated by comparing wet radii against per-bin critical radii at
        ``summary["nd_t_eval"]`` (the final output time, which is
        ``terminate_depth`` m above S_max when ``terminate=True``).
        This is a hard-threshold diagnostic; for the differentiable analog see issue #67.
        """
        return float(self.summary["total_Nd"])

    @property
    def nd_frac(self) -> float:
        """Total activated fraction at the last trajectory step (see :attr:`Nd`)."""
        return float(self.summary["total_nd_frac"])

    # ------------------------------------------------------------------
    # Format conversions
    # ------------------------------------------------------------------

    def to_pandas(self) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
        """Return ``(parcel_df, {species: aerosol_df})`` as pandas DataFrames.

        ``parcel_df`` has columns ``["z", "P", "T", "wv", "wc", "wi", "S"]``
        indexed by simulation time.  Each ``aerosol_df`` has columns
        ``["r000", "r001", …]`` (wet radii in metres) with the same time index.
        """
        import pandas as pd

        idx = pd.Index(self.time, name="time")
        parcel = pd.DataFrame(
            {var: self.state[:, i] for i, var in enumerate(c.STATE_VARS)},
            index=idx,
        )
        aerosol: dict[str, pd.DataFrame] = {}
        offset = c.N_STATE_VARS
        for aer in self.aerosols:
            cols = {f"r{j:03d}": self.state[:, offset + j] for j in range(aer.nr)}
            aerosol[aer.species] = pd.DataFrame(cols, index=idx)
            offset += aer.nr
        return parcel, aerosol

    def to_polars(self) -> tuple[pl.DataFrame, dict[str, pl.DataFrame]]:
        """Return ``(parcel_df, {species: aerosol_df})`` as polars DataFrames.

        Columns and layout mirror :meth:`to_pandas`; the time index becomes an
        explicit ``"time"`` column (polars does not have a named index).
        """
        import polars as pl

        parcel = pl.DataFrame(
            {"time": self.time} | {var: self.state[:, i] for i, var in enumerate(c.STATE_VARS)}
        )
        aerosol: dict[str, pl.DataFrame] = {}
        offset = c.N_STATE_VARS
        for aer in self.aerosols:
            cols = {"time": self.time} | {
                f"r{j:03d}": self.state[:, offset + j] for j in range(aer.nr)
            }
            aerosol[aer.species] = pl.DataFrame(cols)
            offset += aer.nr
        return parcel, aerosol

    def to_xarray(self) -> xr.Dataset:
        """Return a CF-flavoured :class:`xarray.Dataset`.

        Mirrors the variable layout of the legacy NetCDF writer: a ``time``
        coordinate, per-species ``<species>_bins`` coordinates with dry radii /
        kappa / number concentration, wet-radius histories ``<species>_size``,
        parcel thermodynamic profiles, and the post-solve summary as scalar
        variables / global attributes.
        """
        import xarray as xr

        from . import __version__ as ver

        ds = xr.Dataset(attrs={"Conventions": "CF-1.0", "source": f"pyrcel v{ver} (JAX/diffrax)"})
        ds.coords["time"] = (
            "time",
            self.time,
            {"units": "seconds", "long_name": "simulation time"},
        )

        offset = c.N_STATE_VARS
        for aer in self.aerosols:
            nr = aer.nr
            sp = aer.species
            bins = f"{sp}_bins"
            ds.coords[bins] = (
                bins,
                np.arange(1, nr + 1, dtype=np.int32),
                {"long_name": f"{sp} size bin number"},
            )
            ds[f"{sp}_rdry"] = (
                (bins,),
                np.asarray(aer.r_drys) * 1e6,
                {"units": "micron", "long_name": f"{sp} bin dry radii"},
            )
            ds[f"{sp}_kappas"] = (
                (bins,),
                np.full(nr, aer.kappa),
                {"long_name": f"{sp} bin kappa-kohler hygroscopicity"},
            )
            ds[f"{sp}_Nis"] = (
                (bins,),
                np.asarray(aer.Nis) * 1e-6,
                {"units": "cm-3", "long_name": f"{sp} bin number concentration"},
            )
            ds[f"{sp}_size"] = (
                ("time", bins),
                self.state[:, offset : offset + nr] * 1e6,
                {"units": "micron", "long_name": f"{sp} bin wet radii"},
            )
            offset += nr

        rho = rho_air(self.T, self.P, self.S + 1.0)
        profiles = {
            "S": (self.S * 100.0, {"units": "%", "long_name": "Supersaturation"}),
            "T": (self.T, {"units": "K", "long_name": "Temperature"}),
            "P": (self.P, {"units": "Pa", "long_name": "Pressure"}),
            "wv": (self.wv, {"units": "kg/kg", "long_name": "Water vapor mixing ratio"}),
            "wc": (self.wc, {"units": "kg/kg", "long_name": "Liquid water mixing ratio"}),
            "wi": (self.wi, {"units": "kg/kg", "long_name": "Ice water mixing ratio"}),
            "height": (self.heights, {"units": "meters", "long_name": "Parcel height above start"}),
            "rho": (np.asarray(rho), {"units": "kg/m3", "long_name": "Air density"}),
            "wtot": (
                self.wv + self.wc,
                {"units": "kg/kg", "long_name": "Total water mixing ratio"},
            ),
        }
        for name, (data, attrs) in profiles.items():
            ds[name] = (("time",), np.asarray(data), attrs)

        s = self.summary
        ds["S_max"] = (
            (),
            s["S_max"] * 100.0,
            {"units": "%", "long_name": "Maximum supersaturation"},
        )
        ds["t_smax"] = (
            (),
            s["t_smax"],
            {"units": "seconds", "long_name": "Time of maximum supersaturation"},
        )
        for p in s["per_species"]:
            ds[f"{p['species']}_eq_act_frac"] = (
                (),
                p["eq_act_frac"],
                {"long_name": f"{p['species']} equilibrium activated fraction"},
            )
            ds[f"{p['species']}_Nd"] = (
                (),
                p["Nd"],
                {"units": "cm-3", "long_name": f"{p['species']} activated droplet number"},
            )
        ds["Nd"] = (
            (),
            s["total_Nd"],
            {"units": "cm-3", "long_name": "Total activated droplet number concentration"},
        )
        ds["nd_t_eval"] = (
            (),
            s["nd_t_eval"],
            {"units": "s", "long_name": "Time at which Nd snapshot was evaluated"},
        )

        if isinstance(self.V, ConstantV):
            v_attr = float(self.V.V)
        elif isinstance(self.V, AbstractUpdraft):
            v_attr = "V(t)"
        else:
            v_attr = float(self.V)
        ds.attrs.update(
            {
                "V": v_attr,
                "T0": self.T0,
                "S0": self.S0,
                "P0": self.P0,
                "accom": self.accom,
                "total_act_frac": s["total_act_frac"],
            }
        )
        return ds

    def to_netcdf(self, path: str | Path) -> str:
        """Write to a NetCDF4 file and return the path."""
        path = str(path)
        self.to_xarray().to_netcdf(path)
        return path

    def to_csv(self, path: str | Path) -> str:
        """Write the flat parcel trajectory (all state columns + radius bins) to CSV.

        Columns: ``time``, then all state variables, then per-bin wet radii
        prefixed by species (e.g. ``sulfate_r000``).
        """
        parcel, aerosol = self.to_pandas()
        # Flatten aerosol DFs with species-prefixed column names
        import pandas as pd

        parts = [parcel.reset_index()]
        for species, df in aerosol.items():
            parts.append(
                df.reset_index(drop=True).rename(
                    columns={col: f"{species}_{col}" for col in df.columns}
                )
            )
        flat = pd.concat(parts, axis=1)
        flat.to_csv(path, index=False)
        return str(path)

    def to_parquet(self, path: str | Path) -> str:
        """Write the flat parcel trajectory (all state columns + radius bins) to Parquet."""
        parcel_pl, aerosol_pl = self.to_polars()
        import polars as pl

        parts = [parcel_pl]
        for species, df in aerosol_pl.items():
            parts.append(
                df.drop("time").rename(
                    {col: f"{species}_{col}" for col in df.columns if col != "time"}
                )
            )
        flat = pl.concat(parts, how="horizontal")
        flat.write_parquet(path)
        return str(path)
