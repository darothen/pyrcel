# Migrating from pyrcel v1 to v2

pyrcel v2 replaces the numba + Assimulo/SUNDIALS backend with a pure-Python
[JAX](https://github.com/google/jax) + [diffrax](https://github.com/patrick-kidger/diffrax)
stack. The headline benefits of this replacement (`pip`-installable, differentiable,
GPU-capable, and `jax.vmap`-batchable) required retiring some internal machinery, but
the public interface was mostly left unchanged.

---

## What changed at a glance

| Area | v1 | v2 |
|---|---|---|
| ODE solver | SUNDIALS CVode (BDF) via Assimulo | diffrax `Kvaerno5` (ESDIRK-5/4) |
| RHS JIT | numba `njit` + `pycc` AOT | `jax.jit` (XLA compilation) |
| Installation | `pip install pyrcel` + `conda install assimulo` | `pip install pyrcel` (no native deps) |
| GPU support | — | `device="gpu"` on `ParcelModel` |
| `jax.grad` | — | ✓ through the full integration |
| `jax.vmap` | — | ✓ for ensemble runs |
| `ParcelModel.run()` return | tuple `(parcel_df, aer_dfs)` | `ModelOutput` object |
| Activation parameterizations | module-level functions | `ARG2000()`, `MBN2014()` classes |
| Numerical precision | numpy float64 (effective) | explicit `jax_enable_x64=True` |

Legacy NumPy/SciPy reference implementations for thermodynamics and activation
remain in `pyrcel.legacy` for cross-checking.

---

## Installation

v2 has no compiled native dependencies:

```bash
pip install pyrcel          # CPU (JAX with XLA)
pip install "pyrcel[gpu]"   # CUDA 12 (jax[cuda12])
```

The old Assimulo and numba dependencies are gone entirely. If you previously
installed via a conda environment, a clean `pip install` into a fresh environment
is the cleanest path.

---

## Import changes

Top-level imports are unchanged:

```python
import pyrcel as pm
from pyrcel import AerosolSpecies, Lognorm, ParcelModel
```

`ParcelModelJAX` is a backwards compatible alias for `ParcelModel` and will
continue to work, but new code should use `ParcelModel` directly.

---

## `ParcelModel` constructor

The constructor signature is the same. Two new optional keyword arguments have
been added; existing call-sites that do not pass them are unaffected:

```python
# v1 / still valid in v2
model = pm.ParcelModel(
    [aerosol], V=1.0, T0=283.15, S0=0.0, P0=85000.0
)

# v2 additions
model = pm.ParcelModel(
    [aerosol], V=1.0, T0=283.15, S0=0.0, P0=85000.0,
    console=True,    # print setup / summary tables (was implicit in v1)
    device="gpu",    # dispatch integration to GPU (new in v2)
)
```

Equilibration of the initial wet-radius spectrum still happens at construction
time (mirroring `_setup_run` in v1), but now uses
[optimistix](https://github.com/patrick-kidger/optimistix) bisection internally
rather than `scipy.optimize.bisect`. Numerically this algorithm is identical, but
the new machinery allows us to compute gradients with respect to the analytic form
of the input aerosol size distribution, regardless of binning.

---

## `ParcelModel.run()` — return type and new parameters

The most visible API change: `run()` now returns a `ModelOutput` object rather
than the v1 `(parcel_df, aer_dfs)` tuple.

### v1 pattern
```python
parcel_df, aer_dfs = model.run(
    t_end=300.0, output_dt=1.0, solver="cvode"
)
s_max = parcel_df["S"].max()
```

### v2 pattern
```python
out = model.run(t_end=300.0, output_dt=1.0)

# Access the state trajectory
parcel_df, aer_dfs = out.to_pandas()
s_max = out.summary["S_max"]           # precisely located via dS/dt event
```

The `solver=` parameter is removed (there is only one solver). The `terminate`
and `terminate_depth` parameters are unchanged:

```python
out = model.run(
    t_end=300.0,
    output_dt=1.0,
    terminate=True,          # stop terminate_depth m above S_max (default True)
    terminate_depth=10.0,    # metres of extra ascent after S_max (default 10)
    progress=True,           # diffrax text progress meter
)
```

### Short-circuit modes

Two new `mode` values return a scalar directly, bypassing the full `ModelOutput`
construction — useful on the differentiable path:

```python
smax  = model.run(t_end=300.0, output_dt=1.0, mode="smax")   # float, S_max
nd    = model.run(t_end=300.0, output_dt=1.0, mode="nd")      # float, N_d (cm⁻³)
```

---

## Accessing output

`ModelOutput` replaces the v1 `(parcel_df, aer_dfs)` tuple and adds several
output-format conversions:

```python
out = model.run(t_end=300.0, output_dt=1.0)

# Pandas (same structure as v1)
parcel_df, aer_dfs = out.to_pandas()

# Polars
parcel_pl, aer_pls = out.to_polars()

# xarray Dataset (CF-flavoured, with metadata)
ds = out.to_xarray()

# Write to disk
out.to_netcdf("run.nc")
out.to_csv("run.csv")
out.to_parquet("run.parquet")

# Convenience properties
print(out.S)          # supersaturation trajectory, shape (n_time,)
print(out.T)          # temperature trajectory
print(out.Nd)         # activated droplet number at trajectory end (m⁻³)
print(out.nd_frac)    # activated fraction at trajectory end
```

The post-solve summary dict is accessible through the model directly:

```python
s = model.summary()
print(s["S_max"])          # peak supersaturation (decimal)
print(s["t_smax"])         # time of S_max (s)
print(s["total_act_frac"]) # total activated fraction
```

---

## Activation parameterizations

v1 exposed activation as module-level functions (`pyrcel.activation.arg2000`,
`pyrcel.activation.mbn2014`). v2 wraps each scheme in a callable class:

```python
from pyrcel.activation import ARG2000, MBN2014

arg = ARG2000()
smax, nact, act_frac = arg(V=1.0, T=283.15, P=85000.0,
                            mus=[0.05e-6], sigmas=[2.0],
                            Ns=[1000e6], kappas=[0.54])

mbn = MBN2014()
smax, nact, act_frac = mbn(V=1.0, T=283.15, P=85000.0,
                            mus=[0.05e-6], sigmas=[2.0],
                            Ns=[1000e6], kappas=[0.54])
```

Both classes are thin wrappers over fully JAX-traceable implementations and are
therefore compatible with `jax.grad` and `jax.vmap`. The underlying functions
`pyrcel.activation.arg2000` and `pyrcel.activation.mbn2014` remain importable
for backward compatibility.

---

## Updraft specification

Constant updraft speeds still work as plain floats. Time-varying updrafts now
use the `AbstractUpdraft` hierarchy rather than bare callables:

```python
import numpy as np
from pyrcel import ConstantV, InterpolatedUpdraft, as_updraft

# v1: callable
model = pm.ParcelModel([aer], V=lambda t: 1.0 + 0.5 * np.sin(t / 30.0), ...)

# v2: InterpolatedUpdraft (JAX-traceable, vmap-safe)
ts = np.linspace(0, 300, 1000)
Vs = 1.0 + 0.5 * np.sin(ts / 30.0)
model = pm.ParcelModel([aer], V=InterpolatedUpdraft(ts, Vs), ...)

# Convenience helper: wraps a scalar or (ts, Vs) pair
model = pm.ParcelModel([aer], V=as_updraft(1.0), ...)
```

The `as_updraft` helper accepts a scalar (returns `ConstantV`), a tuple
`(ts, Vs)` (returns `InterpolatedUpdraft`), or an existing `AbstractUpdraft`
(returned as-is).

---

## CLI runner

The v1 `run_parcel <yaml>` CLI is replaced by a v2 equivalent with compatible
YAML format and two new keys:

```yaml
model_control:
  dt_output: 1.0
  t_end: 300.0
  terminate: true
  terminate_depth: 10.0
  device: gpu        # NEW: dispatch to GPU
```

```bash
run_parcel config.yml
run_parcel config.yml --device gpu   # CLI flag overrides YAML
```

---

## `pyrcel.legacy` — what is preserved and why

`pyrcel.legacy` contains the original NumPy/SciPy implementations of
thermodynamics and activation that existed before the v2 rewrite:

| Module | Contents | Use case |
|---|---|---|
| `pyrcel.legacy.thermo` | NumPy thermodynamic functions | Cross-checking, reference |
| `pyrcel.legacy.activation` | NumPy activation parameterizations (`arg2000`, `mbn2014`, `binned_activation`) | Cross-checking v2 JAX output |

These modules are **not** part of the v2 numerical path — they are preserved
solely for validation and cross-checking. They are intentionally not
differentiable and should not be used in new application code.

The v1 CVode/Assimulo integrator is **not** preserved; the SUNDIALS dependency
has been removed entirely. You can still download older versions of the model
directly from GitHub or through [PyPI](https://pypi.org/project/pyrcel/) if you
need to run with the original numerical solver for any reason.

---

## Known numerical differences

The switch from CVode (BDF) to `diffrax.Kvaerno5` (ESDIRK-5/4) produces results
that are not bit-identical to v1 but agree within the solver tolerances used by
both methods ($\text{rtol} = 10^{-7}$). Observed differences in practice:

- **$S_\text{max}$**: typically $|\Delta S_\text{max}| < 10^{-5}$ (absolute) across
  standard test scenarios.
- **Activated fraction**: typically agrees to three significant figures.
- **Trajectory timing**: the precise time of $S_\text{max}$ may differ by $O(0.1)$ s
  due to the different BDF vs. ESDIRK step sequences.

These differences are within the physical uncertainty of the model (accommodation
coefficient, surface tension parameterization). The test suite cross-validates v2
output against frozen golden data generated from v1.

The equilibration step also changes slightly: v1 used `scipy.optimize.bisect`
with a tolerance of $10^{-6}$; v2 uses `optimistix.Bisection` with a tighter
tolerance. Initial wet radii therefore agree to $\lesssim 10^{-10}$ m rather
than $\sim 10^{-6}$ m, which slightly affects the early trajectory but not the
activated fraction.

---

## New capabilities in v2

### Gradient-based sensitivity analysis

Every output of the parcel model is now a differentiable function of its inputs
via `jax.grad`:

```python
import jax
import jax.numpy as jnp
from pyrcel.integrator import max_supersaturation

# ∂S_max/∂V — how sensitive is the peak supersaturation to updraft speed?
grad_fn = jax.grad(lambda V: max_supersaturation(
    y0, (r_drys, Nis, kappas, accom, pm.ConstantV(V)), ts
))
dsmax_dV = float(grad_fn(jnp.float64(1.0)))
```

See `examples/differentiable_smax.py` for a complete worked example
computing a full sensitivity table.

### Batched ensemble runs

`jax.vmap` maps a single model solve over an array of inputs without a Python
loop, using a single JIT-compiled kernel:

```python
smax_arr, nact_arr = pm.smax_nact_ensemble(y0_batch, args_batch, t_end)
```

The `run_updraft_ensemble` convenience function handles the full pipeline for
Gaussian updraft-speed ensembles:

```python
result = pm.run_updraft_ensemble([aerosol], T0=283.15, S0=0.0, P0=85000.0,
                                  mean=0.5, std=0.2, n=1024)
print(result["S_max"].mean())
```

### GPU acceleration

Pass `device="gpu"` to `ParcelModel` or set `JAX_PLATFORM_NAME=gpu` in the
environment. See [issue #41](https://github.com/darothen/pyrcel/issues/41) for
background.

```python
model = pm.ParcelModel([aerosol], V=1.0, T0=283.15, S0=0.0, P0=85000.0,
                        device="gpu")
```
