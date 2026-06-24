# Quick Start

This page walks through the core workflow: constructing an aerosol population,
running a parcel simulation, accessing the output, and computing gradients.

## 1. Construct an aerosol population

```python
import pyrcel as pm

# Single lognormal sulfate mode: μ = 0.05 μm, σ = 2.0, N = 1000 cm⁻³
aerosol = pm.AerosolSpecies(
    "sulfate",
    pm.Lognorm(mu=0.05e-6, sigma=2.0, N=1000.0),
    kappa=0.6,
    bins=100,
)
```

`AerosolSpecies` discretizes the continuous lognormal into 100 size bins and
converts the number concentration to m⁻³ internally. The constructor also accepts
a pre-discretized `dict` with `"r_drys"` (µm) and `"Nis"` (cm⁻³) keys for
custom distributions.

## 2. Run the parcel model

```python
model = pm.ParcelModel(
    [aerosol],
    V=1.0,          # updraft speed (m/s)
    T0=283.15,      # initial temperature (K)
    S0=-0.02,       # initial supersaturation (−2 % = 98 % RH)
    P0=85000.0,     # initial pressure (Pa)
    console=True,   # print setup and activation summary
)
out = model.run(t_end=300.0, output_dt=1.0)
```

`run()` returns a [`ModelOutput`](../api/model_output.md) object. The integration
uses an event-based termination: it stops `terminate_depth` metres (default 10 m)
above the supersaturation maximum, which is precisely located by solving $dS/dt = 0$.

## 3. Access the output

```python
# Scalar diagnostics
print(f"S_max  = {out.summary['S_max']*100:.3f} %")
print(f"t_smax = {out.summary['t_smax']:.1f} s")
print(f"N_act  = {out.Nd:.3e} m⁻³")
print(f"f_act  = {out.nd_frac:.3f}")

# Full trajectory as pandas DataFrames (mirrors v1 output)
parcel_df, aer_dfs = out.to_pandas()

# CF-flavoured xarray Dataset (for NetCDF / analysis)
ds = out.to_xarray()
ds.to_netcdf("run.nc")

# Polars for fast tabular work
parcel_pl, aer_pls = out.to_polars()
```

## 4. Compare with an activation parameterization

```python
from pyrcel.activation import ARG2000, MBN2014

arg = ARG2000()
smax_arg, nact_arg, frac_arg = arg(
    V=1.0, T=283.15, P=85000.0,
    mus=[0.05e-6], sigmas=[2.0], Ns=[1000e6], kappas=[0.6],
)
print(f"ARG2000  S_max = {smax_arg*100:.3f} %,  N_act = {nact_arg:.3e} m⁻³")

mbn = MBN2014()
smax_mbn, nact_mbn, frac_mbn = mbn(
    V=1.0, T=283.15, P=85000.0,
    mus=[0.05e-6], sigmas=[2.0], Ns=[1000e6], kappas=[0.6],
)
print(f"MBN2014  S_max = {smax_mbn*100:.3f} %,  N_act = {nact_mbn:.3e} m⁻³")
```

Both parameterizations are JAX-traceable and support `jax.grad` / `jax.vmap`.

## 5. Compute gradients

Gradients flow through the full ODE integration via the adjoint method. The
low-level `max_supersaturation` function is the differentiable entry point:

```python
import jax
import jax.numpy as jnp
from pyrcel.integrator import max_supersaturation, equilibrate_initial_state

# Build the initial state vector (same as ParcelModel does internally)
r_drys  = aerosol.r_drys
Nis     = aerosol.Nis
kappas  = jnp.full(aerosol.nr, 0.6)
y0      = equilibrate_initial_state(283.15, -0.02, 85000.0, r_drys, kappas, Nis)

ts      = jnp.linspace(0.0, 300.0, 301)

# ∂S_max/∂V — sensitivity to updraft speed
def smax_of_V(V):
    args = (r_drys, Nis, kappas, 1.0, pm.ConstantV(V))
    return max_supersaturation(y0, args, ts)

dSmax_dV = float(jax.grad(smax_of_V)(jnp.float64(1.0)))
print(f"∂S_max/∂V = {dSmax_dV:.4f}  (units: 1 / (m/s))")
```

See the [Sensitivity Analysis example](../examples/differentiable_smax.md) for a
complete gradient table with finite-difference verification.

## Next steps

- [Scientific Description](../user_guide/sci_descr.md) — the ODE system and physics
- [Numerical Methods](../user_guide/numerical_methods.md) — solver, tolerances, adjoints
- [Examples](../examples/basic_run.md) — worked scripts with rendered output
- [API Reference](../api/model.md) — full class and function documentation
