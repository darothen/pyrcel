"""Standalone regressions for this fix on unmodified upstream."""

import numpy as np
import pytest

from pyrcel import AerosolSpecies
from pyrcel.aerosol import aerosol_mode_ids
from pyrcel.model_output import ModelOutput


@pytest.mark.parametrize("fmt", ["pandas", "polars", "xarray", "csv", "parquet", "netcdf"])
def test_duplicate_material_modes_survive_exports(fmt, tmp_path):
    aerosols = [
        AerosolSpecies("sulfate", {"r_drys": [0.05, 0.1], "Nis": [100.0, 200.0]}, 0.6),
        AerosolSpecies("sulfate", {"r_drys": [0.5], "Nis": [50.0]}, 0.3),
    ]
    state = np.array([[0.0, 85000.0, 283.15, 0.009, 0.0, 0.0, -0.01, 1e-7, 2e-7, 3e-7]])
    summary = dict(
        S_max=0.01,
        t_smax=0.0,
        total_Nd=3e8,
        nd_t_eval=0.0,
        total_act_frac=1.0,
        per_species=[
            dict(species="sulfate", eq_act_frac=1.0, Nd=1e8),
            dict(species="sulfate", eq_act_frac=1.0, Nd=2e8),
        ],
    )
    out = ModelOutput(np.array([0.0]), state, aerosols, summary, 0.5, 283.15, -0.01, 85000.0, 1.0)
    mode1, mode2 = out.mode_ids
    if fmt in ("pandas", "polars"):
        _, frames = getattr(out, f"to_{fmt}")()
        assert set(frames) == {mode1, mode2}
        assert frames[mode1].width == 3 if fmt == "polars" else frames[mode1].shape == (1, 2)
    elif fmt in ("xarray", "netcdf"):
        ds = out.to_xarray()
        if fmt == "netcdf":
            import xarray as xr

            out.to_netcdf(tmp_path / "modes.nc")
            ds = xr.load_dataset(tmp_path / "modes.nc", decode_timedelta=False)
        assert ds[f"{mode1}_size"].shape == (1, 2)
        assert ds[f"{mode2}_size"].shape == (1, 1)
    else:
        import pandas as pd

        path = tmp_path / f"modes.{fmt}"
        getattr(out, f"to_{fmt}")(path)
        df = pd.read_csv(path) if fmt == "csv" else __import__("polars").read_parquet(path)
        assert f"{mode1}_r001" in df.columns and f"{mode2}_r000" in df.columns


def test_generated_mode_ids_cannot_collide_with_supplied_names():
    def a(name):
        return AerosolSpecies(name, {"r_drys": [0.05], "Nis": [1.0]}, 0.6)

    labels = aerosol_mode_ids([a("sulfate"), a("sulfate"), a("sulfate__mode1")])
    assert len(set(labels)) == 3 and labels[-1] == "sulfate__mode1"
