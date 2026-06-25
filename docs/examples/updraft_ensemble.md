# Updraft Ensemble

Maps a single parcel solve over an ensemble of updraft speeds sampled from a
Gaussian distribution using `jax.vmap`. All ensemble members run in a single
JIT-compiled batched kernel — no Python loop.

**Script:** [`examples/updraft_ensemble.py`](https://github.com/darothen/pyrcel/blob/master/examples/updraft_ensemble.py)

```bash
python examples/updraft_ensemble.py
python examples/updraft_ensemble.py --mean 1.0 --std 0.4 --n 1024 \
    --out output/ensemble.nc --plot output/ensemble.png
```

## Setup

```python
import jax
import jax.numpy as jnp
import pyrcel as pm
from pyrcel.integrator import max_supersaturation, equilibrate_initial_state

# Build shared initial state
y0   = equilibrate_initial_state(T0, S0, P0, r_drys, kappas, Nis)
ts   = jnp.linspace(0.0, t_end, n_out)

# Sample updraft speeds
V_samples = pm.sample_gaussian_updrafts(mean=0.5, std=0.25, n=256, seed=42)

# Batch over updraft speeds with vmap
def single_run(V):
    args = (r_drys, Nis, kappas, accom, pm.ConstantV(V))
    return max_supersaturation(y0, args, ts)

smax_batch = jax.vmap(single_run)(V_samples)
```

## Console output

```
--8<-- "docs/assets/output/updraft_ensemble.txt"
```

## Output distributions

![S_max and N_act distribution across ensemble](../assets/figures/updraft_ensemble.png)

The figure shows the distributions of $S_\text{max}$ (left) and $N_\text{act}$ (right)
across the 256-member ensemble. The spread reflects the sensitivity of activation to
updraft speed: stronger updrafts cool faster, reach higher supersaturation, and activate
more droplets.

## Convenience wrapper

`run_updraft_ensemble` handles the full pipeline including aerosol setup, equilibration,
and ensemble execution:

```python
result = pm.run_updraft_ensemble(
    [aerosol], T0=283.15, S0=-0.02, P0=85000.0,
    mean=0.5, std=0.25, n=256,
)
# result keys: "S_max", "N_act", "T_smax", "activated", "V"
print(result["S_max"].mean())
```

## Performance notes

- The first `vmap` call triggers JIT compilation (~10–30 s for 100 bins).
  Subsequent calls with the same shapes are fast (< 1 s for n=256 on CPU).
- On GPU, `vmap` runs all members in parallel — typically 50–100× faster than
  a sequential CPU loop for large ensembles.
- `t_end` must be long enough for the *slowest* member to reach its supersaturation
  maximum. The script auto-computes `t_end = z_cap / min(V)`.
