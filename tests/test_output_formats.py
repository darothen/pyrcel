"""Tests for ModelOutput format-conversion methods.

Covers .to_pandas(), .to_polars(), .to_xarray(), .to_netcdf(), .to_csv(), .to_parquet()
and the convenience property accessors (heights, T, P, S, wv, wc, wi).
"""

from __future__ import annotations

import numpy as np
import pytest
import scenarios as scn

from pyrcel.model import ParcelModel
from pyrcel.model_output import ModelOutput

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def output() -> ModelOutput:
    sc = scn.get_scenario("simple_sulfate")
    ic, run = sc["initial"], sc["run"]
    m = ParcelModel(
        scn.build_aerosols(sc),
        V=ic["V"],
        T0=ic["T0"],
        S0=ic["S0"],
        P0=ic["P0"],
        accom=ic["accom"],
    )
    return m.run(run["t_end"], run["output_dt"])


# --- type -----------------------------------------------------------------------


def test_run_full_returns_model_output(output):
    assert isinstance(output, ModelOutput)


# --- property accessors ---------------------------------------------------------


def test_property_accessors(output):
    nt = len(output.time)
    for attr in ("heights", "T", "P", "S", "wv", "wc", "wi"):
        arr = getattr(output, attr)
        assert arr.shape == (nt,), attr
    # all temperatures must be physically sane
    assert np.all(output.T > 200) and np.all(output.T < 310)
    assert np.all(output.P > 5e4) and np.all(output.P < 1.1e5)


# --- to_pandas ------------------------------------------------------------------


def test_to_pandas_structure(output):
    parcel, aerosol = output.to_pandas()
    assert list(parcel.columns) == ["z", "P", "T", "wv", "wc", "wi", "S"]
    assert parcel.index.name == "time"
    assert "sulfate" in aerosol
    assert aerosol["sulfate"].index.name == "time"
    assert aerosol["sulfate"].shape[1] == output.aerosols[0].nr
    # column naming pattern
    assert all(c.startswith("r") for c in aerosol["sulfate"].columns)


def test_to_pandas_values_match_state(output):
    parcel, _ = output.to_pandas()
    np.testing.assert_array_equal(parcel["S"].values, output.S)
    np.testing.assert_array_equal(parcel["T"].values, output.T)
    np.testing.assert_array_equal(parcel.index.values, output.time)


# --- to_polars ------------------------------------------------------------------


def test_to_polars_structure(output):
    import polars as pl

    parcel, aerosol = output.to_polars()
    assert isinstance(parcel, pl.DataFrame)
    assert "time" in parcel.columns
    assert "S" in parcel.columns
    assert "sulfate" in aerosol
    assert "time" in aerosol["sulfate"].columns
    assert aerosol["sulfate"].width - 1 == output.aerosols[0].nr  # -1 for time col


def test_to_polars_values_match_state(output):
    parcel, _ = output.to_polars()
    np.testing.assert_array_almost_equal(parcel["S"].to_numpy(), output.S)
    np.testing.assert_array_almost_equal(parcel["time"].to_numpy(), output.time)


def test_pandas_polars_agree(output):
    pd_parcel, _ = output.to_pandas()
    pl_parcel, _ = output.to_polars()
    np.testing.assert_allclose(
        pd_parcel["T"].values,
        pl_parcel["T"].to_numpy(),
        rtol=1e-12,
    )


# --- to_xarray ------------------------------------------------------------------


def test_to_xarray_structure(output):
    ds = output.to_xarray()
    for var in ("S", "T", "P", "wv", "wc", "wi", "height", "rho", "wtot", "S_max", "t_smax"):
        assert var in ds, var
    assert "time" in ds.coords
    assert "sulfate_bins" in ds.coords
    assert ds["sulfate_size"].dims == ("time", "sulfate_bins")


def test_to_xarray_units(output):
    ds = output.to_xarray()
    # S is stored in percent
    assert float(ds["S"].max()) <= float(ds["S_max"]) + 1e-9
    np.testing.assert_allclose(
        float(ds["S_max"]),
        output.summary["S_max"] * 100.0,
        rtol=1e-12,
    )


def test_to_xarray_global_attrs(output):
    ds = output.to_xarray()
    for attr in ("V", "T0", "S0", "P0", "accom", "total_act_frac", "Conventions", "source"):
        assert attr in ds.attrs, attr


# --- to_netcdf ------------------------------------------------------------------


def test_to_netcdf_roundtrip(output, tmp_path):
    import xarray as xr

    path = tmp_path / "run.nc"
    returned = output.to_netcdf(path)
    assert returned == str(path)
    assert path.exists()

    reloaded = xr.open_dataset(str(path))
    np.testing.assert_allclose(
        reloaded["S"].values,
        output.to_xarray()["S"].values,
    )
    reloaded.close()


# --- to_csv ---------------------------------------------------------------------


def test_to_csv_roundtrip(output, tmp_path):
    import pandas as pd

    path = tmp_path / "run.csv"
    returned = output.to_csv(path)
    assert returned == str(path)
    assert path.exists()

    flat = pd.read_csv(path)
    assert "time" in flat.columns
    for var in ("z", "P", "T", "wv", "wc", "wi", "S"):
        assert var in flat.columns, var
    # per-species radius columns are prefixed
    assert any(col.startswith("sulfate_r") for col in flat.columns)
    np.testing.assert_allclose(flat["S"].values, output.S, rtol=1e-12)


# --- to_parquet -----------------------------------------------------------------


def test_to_parquet_roundtrip(output, tmp_path):
    import polars as pl

    path = tmp_path / "run.parquet"
    returned = output.to_parquet(path)
    assert returned == str(path)
    assert path.exists()

    flat = pl.read_parquet(path)
    assert "time" in flat.columns
    for var in ("z", "P", "T", "wv", "wc", "wi", "S"):
        assert var in flat.columns, var
    assert any(col.startswith("sulfate_r") for col in flat.columns)
    np.testing.assert_allclose(flat["S"].to_numpy(), output.S, rtol=1e-12)
