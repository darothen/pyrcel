# Basic Run

A single parcel simulation from initial conditions through droplet activation.
Demonstrates the core API: constructing an aerosol population, running the model,
and accessing the trajectory and activation diagnostics.

**Script:** [`examples/basic_run.py`](https://github.com/darothen/pyrcel/blob/master/examples/basic_run.py)

```bash
python examples/basic_run.py
python examples/basic_run.py --V 0.5 --N 2000 --mu 0.05 --kappa 0.54 \
    --out output/run.nc --plot output/run.png
```

## Setup

```python
import pyrcel as pm

aerosol = pm.AerosolSpecies(
    "sulfate",
    pm.Lognorm(mu=0.05e-6, sigma=2.0, N=1000.0),
    kappa=0.54,
    bins=100,
)

model = pm.ParcelModel(
    [aerosol], V=1.0, T0=283.15, S0=-0.02, P0=85000.0, console=True
)
out = model.run(t_end=300.0, output_dt=1.0)
```

## Activation summary

```
--8<-- "docs/assets/output/basic_run.txt"
```

## Output figure

![Parcel trajectory and aerosol size evolution](../assets/figures/basic_run.png)

**Left panel** — supersaturation $S$ (blue, bottom axis) and temperature $T$ (red, top
axis) as functions of height. The supersaturation rises as the parcel lifts, peaks at
$S_\text{max}$ (dashed gold line), then falls as condensation depletes water vapor faster
than adiabatic cooling can generate it.

**Right panel** — size evolution of a sub-sample of aerosol bins from both modes.
Bold solid traces are activated droplets ($s_\text{crit} < S_\text{max}$); dashed traces
remain as haze. The kink in each activated trace marks the moment the droplet crosses its
Köhler critical radius and begins growing freely.

## Saving output

```python
# NetCDF (CF-flavoured xarray dataset)
out.to_netcdf("run.nc")

# Pandas DataFrames (mirrors v1 tuple output)
parcel_df, aer_dfs = out.to_pandas()

# Parquet for fast downstream analysis
out.to_parquet("run.parquet")
```

See [`ModelOutput`](../api/model_output.md) for the full list of format conversions.
