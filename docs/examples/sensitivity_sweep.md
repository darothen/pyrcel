# Sensitivity Sweep: Exact vs Numerical ∂Smax/∂V

Compares the sensitivity of peak supersaturation $\partial S_\text{max}/\partial V$
computed two ways — via the **exact adjoint (jax.grad)** and via **numerical
central-difference** approximation — for each of three methods:

| Method | Implementation | Exact gradient source |
|--------|---------------|-----------|
| Parcel model | Full ODE (JAX/diffrax) | `jax.grad` through the integrator (adjoint) |
| ARG2000 | Closed-form | `jax.grad` (analytical) |
| MBN2014 | Iterative | `jax.grad` (analytical) |

The numerical gradient is computed identically for all three methods:
`numpy.gradient` applied to the pre-computed $S_\text{max}$ grid.  Comparing
the two columns reveals where finite-difference approximation on a coarse grid
diverges from the true gradient — and highlights why differentiable
parameterizations are valuable even when a grid-based approximation seems
sufficient.

**Script:** [`examples/sensitivity_sweep.py`](https://github.com/darothen/pyrcel/blob/master/examples/sensitivity_sweep.py)

```bash
python examples/sensitivity_sweep.py              # compute, cache, plot
python examples/sensitivity_sweep.py --recompute  # force fresh sweep
python examples/sensitivity_sweep.py --no-plot    # cache only
python examples/sensitivity_sweep.py --cache PATH --plot PATH
```

## Caching

The parcel model requires one JIT-compiled forward pass and one adjoint call
per grid point (100 points on the default 10 × 10 grid). Results are saved to
a compressed `.npz` cache so the figure can be regenerated without re-running
the integrator:

```python
import numpy as np

data = np.load("output/sensitivity_sweep_cache.npz")
smax_parcel            = data["smax_parcel"]            # shape (n_V, n_mu)
dsmax_dV_parcel_exact  = data["dsmax_dV_parcel_exact"]  # exact adjoint
dsmax_dV_parcel_num    = data["dsmax_dV_parcel_num"]    # numpy.gradient
V_grid                 = data["V_grid"]
mu_grid                = data["mu_grid"]
```

The cache is invalidated automatically if the grid dimensions or aerosol
parameters differ between runs.

## How adjoint and numerical gradients are computed

For the parameterizations, the exact gradient is a single `jax.grad` call:

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

dS_dV_exact = jax.grad(smax_arg, argnums=0)(1.0, 0.05)
```

For the parcel model, the adjoint is computed the same way but the function
wraps the ODE integrator.  To reuse a single JAX compilation across the full
grid, all array-valued inputs are passed as explicit arguments (not closed
over):

```python
from pyrcel.integrator import max_supersaturation
from pyrcel.updraft import ConstantV

def smax_parcel(V_val, y0, r_drys, Nis, kappas_arr):
    return max_supersaturation(
        y0, (r_drys, Nis, kappas_arr, accom, ConstantV(V_val)), ts
    )

grad_parcel = jax.jit(jax.grad(smax_parcel, argnums=0))
dS_dV_exact = grad_parcel(V_val, y0, r_drys, Nis, kappas_arr)
```

The JIT compilation happens once on the first call; all subsequent grid points
reuse the compiled kernel.

The numerical approximation is then identical for every method:

```python
import numpy as np

dsmax_dV_num = np.gradient(smax_grid, V_grid, axis=0)
```

## Output figure

![Exact vs numerical ∂Smax/∂V for parcel model, ARG2000, and MBN2014](../assets/figures/sensitivity_sweep.png)

**Layout:** 3 rows × 2 columns. Each row is one method; the left column shows
the exact adjoint gradient and the right column shows the numerical
finite-difference estimate. Colorbars are shared within each row so that
discrepancies between exact and numerical are immediately visible.

**Parcel model (row 1).** The adjoint gradient is obtained by backpropagating
through the full ODE solver — an expensive but exact operation. The numerical
estimate (right) recovers the same broad pattern from the forward-pass grid
alone.  Discrepancies are largest near steep gradients at the boundaries of the
domain, where the 10-point grid is coarsest.

**ARG2000 (row 2).** The closed-form parameterization is fully differentiable
at negligible cost. For this method the exact and numerical gradients agree
closely across the entire domain, providing a useful sanity check for the
grid-based approach.

**MBN2014 (row 3).** Similar to ARG2000. Any residual mismatch between columns
reflects the finite grid spacing rather than errors in either gradient method.
