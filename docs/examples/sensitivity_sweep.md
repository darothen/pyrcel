# Sensitivity Sweep

Maps $S_\text{max}$ and its gradients $\partial S_\text{max}/\partial V$ and
$\partial S_\text{max}/\partial \mu$ across a dense 2D grid of updraft speed $V$
and lognormal median radius $\mu$ for a single-mode sulfate aerosol population.
Three methods are compared side by side:

| Method | Implementation | Gradients |
|--------|---------------|-----------|
| Parcel model | Full ODE (JAX/diffrax) | $\partial S/\partial V$ via `jax.grad`; $\partial S/\partial \mu$ via `numpy.gradient` on grid |
| ARG2000 | Closed-form | `jax.grad` (exact) |
| MBN2014 | Iterative | `jax.grad` (exact) |

**Script:** [`examples/sensitivity_sweep.py`](https://github.com/darothen/pyrcel/blob/master/examples/sensitivity_sweep.py)

```bash
python examples/sensitivity_sweep.py              # compute, cache, plot
python examples/sensitivity_sweep.py --recompute  # force fresh sweep
python examples/sensitivity_sweep.py --no-plot    # cache only
python examples/sensitivity_sweep.py --cache PATH --plot PATH
```

## Caching

The parcel model sweep (25 × 25 grid × 2 quantities per point) takes several minutes on
the first run. Results are saved to a compressed `.npz` cache so the figure can be
regenerated instantly without re-running the integrator:

```python
import numpy as np

data = np.load("output/sensitivity_sweep_cache.npz")
smax_parcel    = data["smax_parcel"]     # shape (25, 25)
dsmax_dV_arg   = data["dsmax_dV_arg"]
dsmax_dmu_mbn  = data["dsmax_dmu_mbn"]
V_grid         = data["V_grid"]
mu_grid        = data["mu_grid"]
```

The cache is invalidated automatically if the grid dimensions or aerosol parameters
change between runs.

## Computing gradients

For the parameterizations, gradients at any $(V, \mu)$ point are exact and cheap:

```python
import jax
import jax.numpy as jnp
from pyrcel.activation import arg2000

def smax_arg(V, mu):
    s, _, _ = arg2000(
        V, 283.0, 85000.0,
        jnp.array([mu]), jnp.array([2.0]),
        jnp.array([1000.0]), jnp.array([0.54]),
    )
    return s

# Gradients at a single point
dS_dV, dS_dmu = jax.grad(smax_arg, argnums=(0, 1))(1.0, 0.05)

# Gradients over a full V × μ grid via double-vmap
grid_grads = jax.vmap(
    jax.vmap(jax.grad(smax_arg, argnums=(0, 1)), in_axes=(None, 0)),
    in_axes=(0, None),
)(V_grid, mu_grid)
```

For the parcel model, $\partial S_\text{max}/\partial V$ is computed the same way,
replacing `smax_arg` with a wrapper around
[`max_supersaturation`](../api/integrator.md):

```python
from pyrcel.integrator import max_supersaturation
import pyrcel as pm

def smax_parcel(V_val):
    args = (r_drys, Nis, kappas, accom, pm.ConstantV(V_val))
    return max_supersaturation(y0, args, ts)

dS_dV_parcel = jax.grad(smax_parcel)(V_val)
```

## Output figure

![Smax and sensitivity maps across (V, μ) space](../assets/figures/sensitivity_sweep.png)

Each row shows a different quantity; each column a different method. All panels share
a common colour scale within their row, so differences are directly comparable.

**Row 1 — $S_\text{max}$.**  All three methods predict higher peak supersaturation at
larger $V$ and smaller $\mu$ (smaller particles are harder to activate, so more
supersaturation builds before the condensation sink catches up). ARG2000 consistently
underestimates $S_\text{max}$, most visibly at high $V$ and intermediate $\mu$.

**Row 2 — $\partial S_\text{max}/\partial V$.**  The sensitivity to updraft speed is
positive everywhere. It peaks at small $\mu$ and large $V$, where the condensation
sink responds most slowly. ARG2000 and MBN2014 agree well with the parcel model in
the bulk of the domain; both tend to underestimate the sensitivity at the small-$\mu$,
high-$V$ corner where non-equilibrium growth is most important.

**Row 3 — $\partial S_\text{max}/\partial \mu$.**  The sensitivity to median radius is
negative: larger particles have a lower critical supersaturation and a larger
condensation sink, both of which suppress $S_\text{max}$. The sign and spatial
structure are consistent across all three methods, though the magnitude differs,
particularly for MBN2014 at low $V$.
