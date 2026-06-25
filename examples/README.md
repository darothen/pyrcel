# pyrcel examples

Demonstration scripts for the v2 (JAX/diffrax) parcel-model backend.

```bash
pip install "pyrcel[examples]"   # adds matplotlib for figures
```

---

## `basic_run.py` — single simulation

Equilibrates and integrates one parcel, prints the activation summary, writes a
NetCDF file, and (with matplotlib) saves a supersaturation/temperature-vs-height figure.

```bash
python examples/basic_run.py
python examples/basic_run.py --V 1.0 --N 1000 --bins 100 \
    --out output/basic_run.nc --plot output/basic_run.png
```

## `live_run.py` — live integration output

Prints altitude / temperature / supersaturation after each integration chunk, mirroring
the master-branch CVode interactive output. Useful for debugging or demonstrating
the model interactively.

```bash
python examples/live_run.py
python examples/live_run.py --chunk-dt 5 --no-console
```

## `updraft_ensemble.py` — Gaussian updraft ensemble via `jax.vmap`

Samples `n` vertical velocities from `Normal(mean, std)` and runs a single
`jax.vmap`-ed ensemble of parcel simulations (one compiled, batched call) to estimate
the output distributions of peak supersaturation $S_\text{max}$ and activated droplet
number $N_\text{act}$. Each member integrates to its own supersaturation maximum via the
`dS/dt` event.

```bash
python examples/updraft_ensemble.py --mean 0.5 --std 0.25 --n 512
python examples/updraft_ensemble.py --mean 1.0 --std 0.4 --n 1024 \
    --out output/ensemble.nc --plot output/ensemble.png
```

Notes:
- `t_end` defaults to `z_cap / min(V)` so even the slowest sampled member reaches its
  supersaturation maximum before the integration stops.
- Increase `--bins` to smooth the $N_\text{act}$ histogram, which is banded at coarse
  resolution because activation is counted over discrete size bins.

## `differentiable_smax.py` — sensitivity analysis via `jax.grad`

Computes $\partial S_\text{max}/\partial\{V, N, \kappa\}$ by propagating gradients
through the full ODE integration, then verifies the autodiff gradient against a
central-difference finite-difference check. This is the showpiece for pyrcel v2's
differentiability.

```bash
python examples/differentiable_smax.py
python examples/differentiable_smax.py --V 0.5 --N 500 --bins 50
```
