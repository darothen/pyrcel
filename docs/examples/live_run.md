# Live Integration

Demonstrates `live=True` mode, which prints a z/T/S diagnostic table after each
integration chunk — replicating the interactive output of the v1 CVode integrator.
Useful for watching long runs or debugging initial conditions.

**Script:** [`examples/live_run.py`](https://github.com/darothen/pyrcel/blob/master/examples/live_run.py)

```bash
python examples/live_run.py
python examples/live_run.py --chunk-dt 5 --no-console
```

## Usage

```python
import pyrcel as pm

model = pm.ParcelModel([aerosol], V=1.0, T0=283.15, S0=-0.02, P0=85000.0)

out = model.run(
    t_end=300.0,
    output_dt=1.0,
    live=True,           # print z/T/S after each chunk
    live_chunk_dt=10.0,  # chunk length in simulation seconds
)
```

Each chunk completes one JAX solve covering `live_chunk_dt` seconds of simulation
time, then prints a summary row before starting the next chunk.

## Console output

```
--8<-- "docs/assets/output/live_run.txt"
```

## Output figure

![Live integration trajectory](../assets/figures/live_run.png)

The figure shows the same supersaturation and temperature profile as the basic run,
confirming that chunked live mode produces an identical trajectory to a single compiled
solve.

## When to use `live` vs. `progress`

| Mode | How it works | Best for |
|---|---|---|
| `live=True` | Python loop over chunks; prints after each | Watching long runs; debugging |
| `progress=True` | diffrax text meter; single compiled call | Interactive sessions wanting a progress bar |
| Neither | Silent single call | Production / `jax.vmap` / `jax.grad` |

!!! warning "live mode is not differentiable"
    Because `live=True` splits the solve across Python-level chunk boundaries,
    it is incompatible with `jax.grad` and `jax.vmap`. Use the low-level
    [`max_supersaturation`](../api/integrator.md) or
    [`nd_from_parcel`](../api/integrator.md) functions for the differentiable path.
