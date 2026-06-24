# pyrcel

**pyrcel** is a zero-dimensional adiabatic cloud parcel model for studying aerosol–cloud
interactions. It simulates the condensational growth of an arbitrary aerosol population as
a parcel rises adiabatically, predicting the peak supersaturation and the number of
droplets activated at cloud base.

[![DOI](https://zenodo.org/badge/12927551.svg)](https://zenodo.org/badge/latestdoi/12927551)
[![PyPI](https://badge.fury.io/py/pyrcel.svg)](https://badge.fury.io/py/pyrcel)
[![CI](https://github.com/darothen/pyrcel/actions/workflows/ci.yml/badge.svg)](https://github.com/darothen/pyrcel/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/pyrcel/badge/?version=stable)](https://pyrcel.readthedocs.io/en/stable/)

## Highlights

**Differentiable.** Every output is a smooth function of its inputs. Compute
$\partial S_\text{max}/\partial V$, $\partial N_d/\partial \kappa$, or full Jacobians
through the ODE using `jax.grad` and `jax.jacfwd`.

**GPU-capable.** Pass `device="gpu"` to `ParcelModel` and the integration dispatches
to CUDA with no code changes. Particularly useful for ensemble or batch simulations!

**Batchable.** Use `jax.vmap` to map a single compiled kernel over an ensemble of
initial conditions — thousands of parcel trajectories in a single call.

**Pure Python.** No compiled extensions or conda dependencies. Install with `pip`.

## Quick install

[`uv`](https://docs.astral.sh/uv/) is the recommended way to manage Python projects
that use pyrcel.

**Add to an existing `uv` project:**

```bash
uv add pyrcel          # CPU
uv add "pyrcel[gpu]"   # CUDA 12
```

**New project from scratch:**

```bash
uv init myproject && cd myproject
uv add pyrcel
uv run python main.py
```

**Directly from GitHub (latest unreleased):**

```bash
uv add "pyrcel @ git+https://github.com/darothen/pyrcel.git"
```

**Editable install from a local clone:**

```bash
git clone https://github.com/darothen/pyrcel.git && cd pyrcel
uv sync          # installs all core deps into an isolated .venv
uv run python    # runs Python inside that environment
```

**pip (without uv):**

```bash
pip install pyrcel
pip install "pyrcel[gpu]"   # CUDA 12
```

---

## Minimal example

```python
import pyrcel as pm

aerosol = pm.AerosolSpecies(
    "sulfate",
    pm.Lognorm(mu=0.05e-6, sigma=2.0, N=1000.0),
    kappa=0.6,
    bins=100,
)
model = pm.ParcelModel([aerosol], V=1.0, T0=283.15, S0=-0.02, P0=85000.0)
out   = model.run(t_end=300.0, output_dt=1.0)

print(f"S_max  = {out.summary['S_max']*100:.3f} %")
print(f"N_act  = {out.Nd:.3e} m⁻³")
```

---

## Cite

If you use pyrcel for research, please cite the original manuscrip where
the model was detailed:

> Rothenberg, D., & Wang, C. (2016). Metamodeling of droplet activation for global
> climate models. *J. Atmos. Sci.*, **73**(4), 1255–1272.
> doi:[10.1175/JAS-D-15-0223.1](https://doi.org/10.1175/JAS-D-15-0223.1)

Additionally, please consider citing the bespoke DOI for the
[specific release version of pyrcel](https://zenodo.org/records/20693507)
that you used during your research (or the base version you modified). This
allows us to track adoption and use of specific model versions over time.

---

## Contents

<div class="grid cards" markdown>

- :material-rocket-launch: **[Getting Started](getting_started/installation.md)**

    Install pyrcel and run your first simulation in minutes.

- :material-book-open-variant: **[User Guide](user_guide/sci_descr.md)**

    Scientific background, migration from v1, numerical methods, and GPU setup.

- :material-code-braces: **[Examples](examples/basic_run.md)**

    Worked examples with rendered output: basic runs, live integration,
    ensemble sweeps, and gradient-based sensitivity analysis.

- :material-api: **[API Reference](api/model.md)**

    Complete reference for all public classes and functions.

</div>
