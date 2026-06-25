# GPU Setup

pyrcel v2 dispatches the parcel integration to GPU via JAX's XLA backend.
The same Python code runs on CPU and GPU — only the `device=` argument
(or an environment variable) changes.

## Installation

```bash
pip install "pyrcel[gpu]"   # installs jax[cuda12]
```

This replaces the CPU JAX build with the CUDA 12 variant. Ensure your CUDA
driver is ≥ 525 and that `nvidia-smi` shows a compatible GPU.

For CUDA 11 or custom CUDA versions, install JAX manually before pyrcel:

```bash
pip install "jax[cuda11_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
pip install pyrcel
```

## Using the GPU

Pass `device="gpu"` to `ParcelModel`:

```python
import pyrcel as pm

model = pm.ParcelModel(
    [aerosol], V=1.0, T0=283.15, S0=-0.02, P0=85000.0,
    device="gpu",
)
out = model.run(t_end=300.0, output_dt=1.0)
```

Or set the environment variable to make all calls use the GPU without changing
any code:

```bash
JAX_PLATFORM_NAME=gpu python examples/basic_run.py
```

## float64 on CUDA

JAX defaults to float32 on GPU for performance. pyrcel forces float64 at import
time (`jax.config.update("jax_enable_x64", True)`) because the parcel ODE
requires double precision — the $r^3 - r_d^3$ cancellation in the condensation
equation loses significant digits at single precision for small droplets.

On NVIDIA GPUs, float64 throughput is typically 1/32 of float32 (consumer GPUs)
to 1/2 (data-centre GPUs like A100/H100). If throughput is a concern, profile
on your specific hardware; for most parcel-model workloads the bottleneck is
memory bandwidth, not FLOPS.

To confirm float64 is active on your device:

```python
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)
x = jnp.float64(1.0)
print(x.dtype)   # float64
```

## Ensemble runs on GPU

`jax.vmap` + GPU is the highest-throughput path for large ensembles:

```python
import jax
import jax.numpy as jnp
from pyrcel.integrator import max_supersaturation

V_samples = jnp.linspace(0.1, 3.0, 1024)

def single_run(V):
    args = (r_drys, Nis, kappas, accom, pm.ConstantV(V))
    return max_supersaturation(y0, args, ts)

smax_batch = jax.vmap(single_run)(V_samples)   # 1024 parcel runs, one kernel
```

On a modern A100, this runs ~100× faster than a sequential CPU loop.

## Device selection via `default_device`

For more fine-grained control, use JAX's `default_device` context manager:

```python
import jax

gpu = jax.devices("gpu")[0]
with jax.default_device(gpu):
    out = model.run(t_end=300.0, output_dt=1.0)
```

## Troubleshooting

**`RuntimeError: No GPU/TPU found`** — JAX cannot detect a CUDA device. Check
that `nvidia-smi` works and that you installed `jax[cuda12]` (not plain `jax`).

**float32 outputs despite `jax_enable_x64=True`** — this usually means the
config was set too late (after JAX already initialized). Import pyrcel before
any other JAX imports, or set `JAX_ENABLE_X64=1` in the environment before
launching Python.

**Out of memory on large ensembles** — the adjoint backward pass
(`RecursiveCheckpointAdjoint`) uses $O(\sqrt{N})$ memory in the number of
solver steps. For very large ensembles, reduce `max_steps` or run in smaller
batches.
