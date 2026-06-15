# JAX / diffrax demos

Demonstration scripts for the v2 (JAX/diffrax) parcel-model backend. Install the optional
dependencies first:

```bash
pip install "pyrcel[jax]"          # core JAX backend
pip install "pyrcel[jax,examples]" # + matplotlib for the figures
```

## `basic_run.py` — single simulation

Equilibrates and integrates one parcel, prints the summary table, writes a NetCDF file,
and (with matplotlib) saves a supersaturation/temperature-vs-height figure.

```bash
python examples/jax/basic_run.py --V 1.0 --N 1000 --bins 100 \
    --out output/jax_basic_run.nc --plot output/jax_basic_run.png
```

## `updraft_ensemble.py` — Gaussian updraft ensemble (`vmap`)

Samples `n` vertical velocities from `Normal(mean, std)` and runs a single `jax.vmap`-ed
ensemble of parcel simulations (one compiled, batched call) to estimate the output
distributions of peak supersaturation `S_max` and activated droplet number `N_act`. Each
member integrates to its own supersaturation maximum via the `dS/dt` event.

```bash
python examples/jax/updraft_ensemble.py --mean 0.5 --std 0.25 --n 512 \
    --out output/ensemble.nc --plot output/ensemble.png
```

Notes:

- `N_act` is the activated droplet number concentration (cm⁻³), counted with the
  equilibrium criterion (modal critical supersaturation below `S_max`) at the peak
  temperature. Its histogram is *banded* at coarse `--bins` because activation is counted
  over discrete size bins; increase `--bins` to smooth it.
- `t_end` defaults to `z_cap / min(V)` so even the slowest sampled member reaches its
  maximum; event termination bounds the actual work per member regardless.
- The programmatic entry point is `pyrcel.run_updraft_ensemble(...)` (see
  `pyrcel/ensemble.py`), with `pyrcel.smax_nact_ensemble(...)` as the lower-level
  `vmap` core.
```
