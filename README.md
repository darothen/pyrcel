pyrcel: cloud parcel model
==========================

![sample parcel model run](docs/assets/figures/basic_run.png)

[![DOI](https://zenodo.org/badge/12927551.svg)](https://zenodo.org/badge/latestdoi/12927551)[![PyPI Version](https://badge.fury.io/py/pyrcel.svg)](https://badge.fury.io/py/pyrcel)[![CI](https://github.com/darothen/pyrcel/actions/workflows/ci.yml/badge.svg)](https://github.com/darothen/pyrcel/actions/workflows/ci.yml)[![Documentation Status](https://readthedocs.org/projects/pyrcel/badge/?version=stable)](https://pyrcel.readthedocs.io/en/latest/)


`pyrcel` is a simple, adiabatic cloud parcel model for use in
aerosol-cloud interaction studies. [Rothenberg and Wang (2016)](http://journals.ametsoc.org/doi/full/10.1175/JAS-D-15-0223.1) discuss the model in detail and its improvements over [Nenes et al (2001)][nenes2001]:

* Implementation of κ-Köhler theory for condensation physics ([Petters and Kreidenweis, 2007][pk2007])
* Extension of model to handle arbitrary sectional representations of aerosol populations, based on user-controlled empirical or parameterized size distributions
* JAX/[diffrax][diffrax]-based numerical core — differentiable, batchable, GPU-ready, with no Fortran/SUNDIALS dependency

[Detailed documentation is available](https://pyrcel.readthedocs.io/en/latest/), including a [scientific description](https://pyrcel.readthedocs.io/en/latest/user_guide/sci_descr/), [installation details](https://pyrcel.readthedocs.io/en/latest/getting_started/installation/), and a [basic example](https://pyrcel.readthedocs.io/en/latest/examples/basic_run/).

Quick Start
-----------

The easiest way to run `pyrcel` from source is with [`uv`](https://docs.astral.sh/uv/):

```shell
$ git clone https://github.com/darothen/pyrcel.git && cd pyrcel
$ uv run python examples/basic_run.py
```

`uv` will automatically create an isolated environment and install all dependencies.
The first call compiles JAX kernels; subsequent calls are fast.

Usage
-----

```python
import pyrcel as pm

sulfate = pm.AerosolSpecies(
    "sulfate", pm.Lognorm(mu=0.05, sigma=2.0, N=1000.0), kappa=0.54, bins=50
)
model = pm.ParcelModel([sulfate], V=1.0, T0=283.15, S0=-0.02, P0=85000.0, console=True)
output = model.run(t_end=200.0, output_dt=1.0, terminate=True)
print(output.summary["S_max"])
```

Key capabilities:

* **Autodiff** — exact gradients of `S_max` w.r.t. updraft speed, initial conditions,
  accommodation coefficient, and aerosol properties via `jax.grad`.
* **Batching / GPU** — `jax.vmap` runs ensembles of parcels in one compiled call;
  pass `device="gpu"` to `ParcelModel` for CUDA acceleration.
* **Time-varying updraft** — pass a `pyrcel.InterpolatedUpdraft(ts=..., vs=...)` as `V`.
* **Flexible output** — `output.to_pandas()`, `.to_xarray()`, `.to_netcdf()`, `.to_parquet()`.

The differentiable core (`pyrcel.integrator`, `pyrcel.equilibrate`) is usable directly
for `jit`/`grad`/`vmap`; `ParcelModel` is the interactive convenience layer (console
output, progress meter, post-solve summary table).

Installation
------------

**From PyPI (recommended):**

```shell
$ pip install pyrcel
```

JAX (CPU) is included by default. No extras needed for standard use.

**From source:**

```shell
$ git clone https://github.com/darothen/pyrcel.git && cd pyrcel
$ uv sync
$ uv run python examples/basic_run.py
```

**GPU support (CUDA 12):**

```shell
$ pip install "pyrcel[gpu]"
```

Requirements
------------

* Python >= 3.11
* [JAX](https://docs.jax.dev/) >= 0.4.38
* [diffrax](https://docs.kidger.site/diffrax/) >= 0.6.2
* [equinox](https://docs.kidger.site/equinox/) >= 0.11.10
* [optimistix](https://docs.kidger.site/optimistix/) >= 0.0.7
* NumPy, SciPy, pandas, polars, xarray

Development
-----------

Clone the repo and install with dev dependencies:

```shell
$ git clone https://github.com/darothen/pyrcel.git && cd pyrcel
$ uv sync --extra test
$ prek install   # installs the git pre-commit hook (requires prek: https://prek.j178.dev)
```

Run the fast test suite:

```shell
$ uv run pytest tests/ -m "not slow"
```

Lint and format are handled automatically by `prek` on commit, or run manually:

```shell
$ prek run --all-files
```

Please fork this repository if you intend to develop the model further so that the
code's provenance can be maintained.

License / Usage
---------------

[All scientific code should be licensed](http://www.astrobetter.com/the-whys-and-hows-of-licensing-scientific-code/). This code is released under the New BSD (3-clause) [license](LICENSE.md).

If you use this for any scientific work resulting in a publication, please cite our
original publication detailing the model:

```
@article {
      author = "Daniel Rothenberg and Chien Wang",
      title = "Metamodeling of Droplet Activation for Global Climate Models",
      journal = "Journal of the Atmospheric Sciences",
      year = "2016",
      publisher = "American Meteorological Society",
      address = "Boston MA, USA",
      volume = "73",
      number = "3",
      doi = "10.1175/JAS-D-15-0223.1",
      pages= "1255 - 1272",
      url = "https://journals.ametsoc.org/view/journals/atsc/73/3/jas-d-15-0223.1.xml"
}
```


[author_email]: mailto:daniel@danielrothenberg.com
[nenes2001]: https://onlinelibrary.wiley.com/doi/abs/10.1034/j.1600-0889.2001.d01-12.x
[pk2007]: http://www.atmos-chem-phys.net/7/1961/2007/acp-7-1961-2007.html
[diffrax]: https://docs.kidger.site/diffrax/
[jax]: https://docs.jax.dev/
