# Installation

## Requirements

pyrcel v2 requires **Python ≥ 3.11**. There are no compiled native dependencies —
everything installs from PyPI.

## Recommended: uv

[`uv`](https://docs.astral.sh/uv/) is the recommended tool for managing Python
environments and projects. Install it once:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### Add to an existing project

```bash
uv add pyrcel          # CPU
uv add "pyrcel[gpu]"   # CUDA 12
```

### New project from scratch

```bash
uv init myproject && cd myproject
uv add pyrcel
uv run python main.py
```

## pip

```bash
pip install pyrcel          # CPU
pip install "pyrcel[gpu]"   # CUDA 12
```

## GPU install (CUDA 12)

Swap in the CUDA-enabled JAX backend with the `gpu` extra:

```bash
uv add "pyrcel[gpu]"
# or
pip install "pyrcel[gpu]"
```

See the [GPU setup guide](../user_guide/gpu.md) for float64 configuration and
device-placement details.

## From GitHub (latest unreleased)

Install directly from the `master` branch to get unreleased changes:

```bash
uv add "pyrcel @ git+https://github.com/darothen/pyrcel.git"
# or
pip install "git+https://github.com/darothen/pyrcel.git"
```

A specific branch or tag:

```bash
uv add "pyrcel @ git+https://github.com/darothen/pyrcel.git@feat/v2-implementation"
```

## Editable install (local development)

Clone the repo and install in editable mode so source changes take effect immediately:

```bash
git clone https://github.com/darothen/pyrcel.git && cd pyrcel
uv sync                        # installs core deps into an isolated .venv
uv run python examples/basic_run.py
```

With development tools (tests, linting, type checking):

```bash
uv sync --extra test
uv run prek run --all-files    # lint + type-check
uv run pytest -m "not slow"    # fast test suite
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
print(jax.devices())   # cpu:0, or gpu:0 with a CUDA install
```
