"""NetCDF / xarray output writer for ParcelModelJAX (follow-up)."""

from __future__ import annotations

import numpy as np
import pytest

import scenarios as scn
from pyrcel.model_jax import ParcelModelJAX


def _run_simple():
    sc = scn.get_scenario("simple_sulfate")
    ic, run = sc["initial"], sc["run"]
    m = ParcelModelJAX(
        scn.build_aerosols(sc),
        V=ic["V"], T0=ic["T0"], S0=ic["S0"], P0=ic["P0"], accom=ic["accom"],
    )
    m.run(run["t_end"], run["output_dt"], terminate=True, terminate_depth=run["terminate_depth"])
    return m


@pytest.mark.slow
def test_to_dataset_structure():
    m = _run_simple()
    ds = m.to_dataset()
    n_time = m.x.shape[0]
    for var in ("S", "T", "P", "wv", "wc", "wi", "height", "rho", "wtot"):
        assert var in ds, var
        assert ds[var].shape == (n_time,)
    assert "sulfate_bins" in ds.coords
    assert ds["sulfate_size"].shape == (n_time, m._nr)
    assert ds["sulfate_rdry"].shape == (m._nr,)
    # S stored in percent; matches the (event-localized) summary S_max
    assert float(ds["S"].max()) <= float(ds["S_max"]) + 1e-6
    np.testing.assert_allclose(float(ds["S_max"]), m.summary()["S_max"] * 100.0, rtol=1e-12)
    assert ds.attrs["total_act_frac"] == m.summary()["total_act_frac"]


def test_to_dataset_requires_run():
    sc = scn.get_scenario("simple_sulfate")
    ic = sc["initial"]
    m = ParcelModelJAX(
        scn.build_aerosols(sc),
        V=ic["V"], T0=ic["T0"], S0=ic["S0"], P0=ic["P0"], accom=ic["accom"],
    )
    try:
        m.to_dataset()
    except RuntimeError:
        return
    raise AssertionError("expected RuntimeError before run()")


@pytest.mark.slow
def test_save_and_reload_netcdf(tmp_path):
    import xarray as xr

    m = _run_simple()
    path = tmp_path / "run.nc"
    m.save_netcdf(str(path))
    assert path.exists()
    reloaded = xr.open_dataset(str(path))
    np.testing.assert_allclose(reloaded["S"].values, m.to_dataset()["S"].values)
    np.testing.assert_allclose(reloaded["sulfate_size"].values, m.x[:, 7:] * 1e6)
    reloaded.close()
