# Installation

## Requirements

pyrcel v2 requires **Python ≥ 3.11**. There are no compiled native dependencies —
everything installs with `pip`.

## CPU install

```bash
pip install pyrcel
```

This installs the JAX CPU backend (`jax[cpu]`) automatically along with
diffrax, equinox, and optimistix.

## GPU install (CUDA 12)

```bash
pip install "pyrcel[gpu]"
```

Swaps in `jax[cuda12]` and the matching CUDA-enabled diffrax/equinox builds.
See the [GPU setup guide](../user_guide/gpu.md) for float64 and device-placement
details.

## From source

```bash
git clone https://github.com/darothen/pyrcel.git
cd pyrcel
pip install -e .
```

For development (includes test and linting tools):

```bash
pip install -e ".[test]"
uv run prek run --all-files   # lint + type-check
uv run pytest                 # test suite
```

## Optional extras

| Extra | Contents |
|---|---|
| `gpu` | CUDA 12 JAX backend |
| `test` | pytest, hypothesis, pyrefly |
| `examples` | matplotlib for example scripts |
| `docs` | mkdocs + Material + mkdocstrings for building these docs |

## Verifying the install

```python
import pyrcel as pm
print(pm.__version__)

import jax
print(jax.devices())   # should list cpu:0 (or gpu:0 if CUDA install)
```

## uv quick start

The fastest way to run the examples without a permanent install is with
[uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/darothen/pyrcel.git && cd pyrcel
uv run python examples/basic_run.py
```

`uv` creates an isolated environment and installs all dependencies automatically.
The first call compiles JAX kernels; subsequent calls are fast.
