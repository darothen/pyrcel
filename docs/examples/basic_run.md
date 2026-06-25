# Basic Run

A single parcel simulation from initial conditions through droplet activation,
using a two-mode aerosol population (accumulation + Aitken). Demonstrates the
core API: constructing an aerosol population, running the model, and accessing
the trajectory and activation diagnostics.

**Script:** [`examples/basic_run.py`](https://github.com/darothen/pyrcel/blob/master/examples/basic_run.py)

```bash
python examples/basic_run.py
python examples/basic_run.py --V 0.5 --T0 283.0 --t-end 200.0 \
    --out output/run.nc --plot output/run.png
```

## Setup

```python
import pyrcel as pm

accumulation = pm.AerosolSpecies(
    "accumulation",
    pm.Lognorm(mu=0.05, sigma=2.0, N=1000.0),  # mu in µm
    kappa=0.54,
    bins=50,
)
aitken = pm.AerosolSpecies(
    "aitken",
    pm.Lognorm(mu=0.01, sigma=1.6, N=2000.0),  # mu in µm
    kappa=0.3,
    bins=50,
)

model = pm.ParcelModel(
    [accumulation, aitken], V=1.0, T0=283.0, S0=-0.02, P0=85000.0, console=True
)
out = model.run(100.0, output_dt=1.0, terminate=False, progress=True)
```

The accumulation mode ($\mu = 50$ nm, $\sigma = 2.0$, $\kappa = 0.54$) represents
sulfate-like particles that activate readily. The Aitken mode ($\mu = 10$ nm,
$\sigma = 1.6$, $\kappa = 0.3$) is smaller and less hygroscopic; only a small
fraction of Aitken bins exceed their critical supersaturation at this updraft speed.

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

**Right panel** — size evolution of a sub-sample of aerosol bins from both modes
(teal = accumulation, orange = Aitken). Bold solid traces are activated droplets
($S_\text{crit} < S_\text{max}$); dashed traces remain as haze. The kink in each
activated trace marks the moment the droplet crosses its Köhler critical radius and
begins growing freely.

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
