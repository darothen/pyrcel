# Best Practices

Short, opinionated tips for the most common pyrcel v2 workflows. Each section
targets a specific question; follow the cross-links for full derivations.

---

## Choosing output times for `max_supersaturation`

The output-time array `ts` passed to `max_supersaturation` is used both as
the ODE output grid and to locate the coarse peak bracket; the actual maximum
is then found analytically via a Hermite cubic quadratic solve (see
[Numerical Methods — Differentiable $S_\text{max}$](numerical_methods.md#differentiable-s_max-hermite-cubic-interpolation)).
Three rules apply:

**1. The peak must lie strictly inside `[ts[0], ts[-1]]`.**
If `ts[-1]` is before the peak, `max_supersaturation` returns the value at
`ts[-1]`, not the true $S_\text{max}$, and the gradient is wrong. A reliable
upper bound is $z_\text{target} \approx 1500$ m above cloud base:

$$
t_\text{end} = \frac{1500\,\text{m}}{V}, \quad \text{clamped to }
[200\,\text{s},\, 15\,000\,\text{s}]
$$

**2. Fix the array *shape* across calls to share the JIT cache.**
The XLA kernel is respecialised whenever `ts.shape[0]` changes. Keep
`n_output` constant and vary only `t_end`:

```python
N_TS = 600   # constant: one JIT compilation for all V

def ts_for_V(V: float) -> jnp.ndarray:
    t_end = float(max(200.0, min(1500.0 / V, 15000.0)))
    return jnp.linspace(0.0, t_end, N_TS)
```

**3. `t_end ∝ 1/V` does not bias gradients.**
By the envelope theorem, $\partial S_\text{max}/\partial t_\text{end} = 0$
whenever the peak is in the interior, so scaling `t_end` with $V$ introduces
no gradient error.

---

## When to use each diagnostic function

| Task | Function | Notes |
|---|---|---|
| Forward simulation + trajectory | `integrate_parcel` / `integrate_parcel_arrays` | Returns full state history |
| Interactive run with termination | `ParcelModel.run(terminate=True)` | Uses event detection; not differentiable |
| Differentiable $S_\text{max}$ | `max_supersaturation(y0, args, ts)` | Dense interpolant + Newton; `jax.grad`-able |
| Differentiable $N_d$ | `nd_from_parcel(y0, args, t_end)` | Sigmoid soft threshold; `jax.grad`-able |
| Precise $t_{S_\text{max}}$ location | `find_smax(y0, args, t_end)` | Event-based; **not** differentiable |
| MBN2014 parameterisation gradient | `mbn2014_smax(...)` | IFT custom VJP; exact gradient |

For gradient-based workflows (`jax.grad`, `jax.jacfwd`, optimisation loops)
always use `max_supersaturation` or `nd_from_parcel`. `find_smax` and the
`terminate=True` path contain data-dependent branches that JAX cannot
differentiate.

---

## Setting up a gradient computation

A minimal checklist for a correct `jax.grad` call:

```python
import jax
import jax.numpy as jnp
import pyrcel as pm
from pyrcel.integrator import max_supersaturation, atol_vector
from pyrcel.equilibrate import equilibrate_initial_state

# 1. Build initial state (must be JAX-traceable)
y0 = equilibrate_initial_state(T0, S0, P0, r_drys, kappas, Nis)

# 2. Build args tuple
V  = pm.ConstantV(v_val)
args = (r_drys, Nis, kappas, accom, V)

# 3. Choose ts: fixed shape, t_end past S_max
ts = jnp.linspace(0.0, 1500.0 / v_val, 600)

# 4. Warm up JIT (first call compiles; subsequent calls are fast)
_ = jax.block_until_ready(max_supersaturation(y0, args, ts))

# 5. Differentiate
grad_fn = jax.grad(max_supersaturation, argnums=1)
grad_args = grad_fn(y0, args, ts)
dSmax_dV = grad_args[4].V   # ConstantV wraps V as an Equinox module
```

Common mistakes:

- **Passing Python floats instead of JAX arrays** for `r_drys`, `Nis`, or
  `kappas` — wrap with `jnp.asarray(...)` before building `args`.
- **Using `nd_from_parcel` without sufficient `t_end`** — if `t_end` is before
  the activated droplets grow past their critical radii the soft $N_d$ will
  undercount. A rule of thumb is `t_end ≥ 2 × t_{S_\text{max}}`.
- **Calling `jax.grad` inside a `vmap` on the first execution** — always warm
  up the JIT kernel with a scalar call first.

---

## Epsilon choice for `nd_from_parcel`

The sigmoid half-width $\varepsilon$ controls the accuracy / smoothness
trade-off for the differentiable droplet number:

$$
N_d^{\text{soft}} = \sum_i N_i\, \sigma\!\left(\frac{r_i - r_{\text{crit},i}}{\varepsilon}\right)
$$

| $\varepsilon$ (m) | Typical error vs. hard threshold | Gradient magnitude | When to use |
|---|---|---|---|
| $10^{-7}$ (100 nm) | $\lesssim 1\%$ | small near threshold | Smooth landscapes, far from threshold |
| $10^{-8}$ (10 nm) | $\lesssim 10\%$ | moderate | **Default** — balances accuracy and smoothness |
| $10^{-9}$ (1 nm) | $\lesssim 20\%$ | large near threshold | Optimisation requiring sharp $N_d$ |

The residual discrepancy vs. the hard threshold comes from the approximate
Köhler formula used for $r_\text{crit}$, not from $\varepsilon$; reducing
$\varepsilon$ below $10^{-9}$ increases gradient magnitude without improving
accuracy relative to exact activation.

---

## Tolerance tuning

The defaults (`rtol=1e-7`, `atol` per-component vector) match the legacy CVode
configuration and are correct for the vast majority of simulations. Tighten
tolerances when:

- **Running `jax.grad` at high updraft speed or small particle size** — if you
  see `EquinoxRuntimeError: A linear solver returned non-finite output` during
  the adjoint backward pass, try `rtol=1e-9` and radius `atol=1e-14`. See
  [Adjoint conditioning failures](numerical_methods.md#adjoint-conditioning-failures).

- **Very high aerosol concentration** ($N \gtrsim 10^4$ cm$^{-3}$ with $\geq
  200$ bins) — the stiffness ratio grows near activation; tighter tolerances
  suppress spurious oscillations in the radius trajectories.

To override tolerances:

```python
from pyrcel.integrator import atol_vector
import jax.numpy as jnp

nr = r_drys.shape[0]
tight_atol = jnp.concatenate([
    jnp.array([1e-4, 1e-4, 1e-4, 1e-10, 1e-10, 1e-4, 1e-8]),
    jnp.full(nr, 1e-14),
])
smax = max_supersaturation(y0, args, ts, rtol=1e-9, atol=tight_atol)
```

---

## JIT warm-up and timing

The first call to any JAX-compiled function triggers trace + XLA compilation,
which takes 10–30 s for a 200-bin simulation. Always warm up before timing:

```python
import time, jax
from pyrcel.integrator import max_supersaturation

# Warm up — discards result but fills the JIT cache
_ = jax.block_until_ready(max_supersaturation(y0, args, ts))

# Now time the actual computation
t0 = time.perf_counter()
smax = float(jax.block_until_ready(max_supersaturation(y0, args, ts)))
print(f"S_max = {smax:.5f}  [{time.perf_counter() - t0:.3f} s]")
```

The compiled kernel is shape-specialised: any change in `ts.shape[0]` or
`y0.shape[0]` (i.e. bin count) triggers a new compilation. Plan your bin
counts before production sweeps.

---

## Ensemble runs with `jax.vmap`

For a batch of simulations that share the same ODE shape (same number of bins)
but differ in parameters, `jax.vmap` compiles a single batched kernel:

```python
import jax

V_samples = jnp.linspace(0.1, 3.0, 256)

def smax_for_V(V_val):
    ts_v = ts_for_V(float(V_val))  # V-scaled horizon, fixed shape
    new_args = (r_drys, Nis, kappas, accom, pm.ConstantV(V_val))
    return max_supersaturation(y0, new_args, ts_v)

# First warm up a scalar call
_ = jax.block_until_ready(smax_for_V(V_samples[0]))

# Then batch
smax_batch = jax.jit(jax.vmap(smax_for_V))(V_samples)  # shape (256,)
```

`jax.vmap` over `ConstantV` works because `ConstantV` is an
[Equinox](https://github.com/patrick-kidger/equinox) module and its `V`
attribute is a JAX array leaf. For `InterpolatedUpdraft`, batch over the
amplitude scalar rather than the full module.

For convenience, `ParcelModel.run_updraft_ensemble` wraps this pattern for
updraft-speed ensembles; see [Updraft Ensemble](../examples/updraft_ensemble.md).

---

## Gradient vs. finite-difference verification

Before trusting a gradient, verify it against a central finite-difference
estimate:

```python
eps = 1e-4  # perturbation size for the parameter being tested

smax_plus  = float(max_supersaturation(y0, args_plus(eps),  ts))
smax_minus = float(max_supersaturation(y0, args_minus(eps), ts))
fd_grad = (smax_plus - smax_minus) / (2 * eps)

adjoint_grad = float(jax.grad(max_supersaturation, argnums=1)(y0, args, ts)[4].V)

ratio = adjoint_grad / fd_grad
print(f"adjoint / FD = {ratio:.4f}")   # should be 0.999–1.001
```

A ratio outside $[0.99, 1.01]$ indicates one of: insufficient `ts` coverage
(peak not bracketed), adjoint conditioning failure, or a numerical issue with
the perturbation size.  See
[Numerical Methods](numerical_methods.md#adjoint-conditioning-failures)
for a diagnosis checklist.
