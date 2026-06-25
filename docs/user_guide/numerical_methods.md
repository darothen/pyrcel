# Numerical Methods

This page documents the ODE formulation, solver configuration, step-size control,
event detection, and automatic differentiation strategy used by the v2
JAX/diffrax backend. For the physical model equations see
[Scientific Description](sci_descr.md); for practical usage see the
[API reference](../api/model.md).

---

## The ODE system

The parcel model integrates a system of $7 + n_r$ ordinary differential equations,
where $n_r$ is the total number of aerosol size bins across all species. The state
vector is

$$
\mathbf{y} = \bigl[z,\; P,\; T,\; w_v,\; w_c,\; w_i,\; S,\; r_1,\dots,r_{n_r}\bigr]^\top
$$

with the bulk parcel equations (altitude $z$, pressure $P$, temperature $T$, water
mixing ratios $w_v$/$w_c$/$w_i$, supersaturation $S$) driving per-bin droplet
growth equations for the wet radii $r_i$.  The right-hand side
$\mathbf{f}(t, \mathbf{y}, \boldsymbol{\theta})$ — parameterised by aerosol
properties $\boldsymbol{\theta} = (\mathbf{r}_\text{dry}, \mathbf{N}_i,
\boldsymbol{\kappa}, \alpha_c, V)$ — is implemented as a single JAX-native
array operation in `pyrcel.parcel_aux` and is `jax.jit`-compiled on first call.

The bulk state equations follow [[Nenes2001]](#Nenes2001) and [[Ghan2011]](#Ghan2011); see
[Scientific Description](sci_descr.md) for the full derivation. The per-bin
growth equation is

$$
\frac{dr_i}{dt} = \frac{G_i}{r_i}\bigl(S - S_{\mathrm{eq},i}\bigr),
\qquad
G_i = \left(\frac{\rho_w R T}{e_s D_v' M_w} + \frac{L \rho_w \bigl(\tfrac{L M_w}{RT} - 1\bigr)}{k_a' T}\right)^{-1}
$$

where $D_v'$ and $k_a'$ are the non-continuum-corrected vapour diffusivity and
thermal conductivity [[Seinfeld2006]](#Seinfeld2006), and $S_{\mathrm{eq},i}$ is the $\kappa$-Köhler
equilibrium supersaturation [[Petters2007]](#Petters2007).

---

## Solver selection

### Why `Kvaerno5` (ESDIRK-5/4)?

The parcel ODE is mildly stiff near activation, where the rapid growth of
critical-size droplets drives a fast $r_i$ relaxation that is coupled back to the
supersaturation tendency $dS/dt$. The v1 backend used SUNDIALS CVode with a BDF
scheme of order up to 5 and an internally adaptive Newton iteration. The v2
backend uses `diffrax.Kvaerno5` [[Kvaerno2004]](#Kvaerno2004), an Explicit Singly Diagonally
Implicit Runge–Kutta (ESDIRK) method of order 5 with an embedded order-4
error estimator.

`Kvaerno5` is A-stable and L-stable ("stiffly accurate"), making it robust
for the stiff components of the system without requiring a variable-order BDF
controller. The algebraic stage equations are solved with `diffrax`'s default
`VeryChord` Newton root-finder, which inherits the same `rtol`/`atol` vector
used by the step-size controller — there is no separate Newton tolerance to tune.

For a brief comparison with the legacy solver:

| Property | v1 (CVode BDF) | v2 (`Kvaerno5`) |
|---|---|---|
| Method family | Linear multistep (BDF) | Implicit Runge–Kutta (ESDIRK) |
| Maximum order | 5 (adaptive) | 5 (fixed) |
| A-stability | ✓ (BDF ≥ order 3 weakly) | ✓ (all stages) |
| L-stability | — | ✓ |
| Step-size control | Internal (SUNDIALS) | PID controller |
| `jit`/`vmap`/`grad` | — | ✓ |

The accuracy difference between the two solvers, at matching tolerances, is
within the physical uncertainty of the model; see
[Numerical accuracy](#numerical-accuracy) below.

### When to consider a different solver

For the vast majority of parcel simulations `Kvaerno5` with the default
tolerances is the right choice. Two situations where you might deviate:

- **Very high aerosol loading** ($N \gtrsim 10^4\,\text{cm}^{-3}$ with many
  bins): the stiffness ratio grows, and reducing `rtol` to $10^{-8}$ and
  `atol` for radii to $10^{-14}$ m may be needed to suppress spurious oscillation
  near activation.
- **Coarse output cadence / speed over accuracy**: `diffrax.Tsit5` (explicit
  Runge–Kutta 5/4, non-stiff) runs faster at reduced accuracy when the simulation
  is below the cloud base where the ODE is effectively non-stiff. In practice
  the adaptive step-size controller already handles this automatically with
  `Kvaerno5`, so switching solvers is rarely necessary.
- **Adjoint failure at extreme conditions**: when `jax.grad` is called at
  very high updraft speeds ($V \gtrsim 2$ m/s) combined with very small median
  radii ($\mu \lesssim 0.03\,\mu$m), the implicit solver's Newton step
  Jacobian can become nearly singular, causing the adjoint backward pass to
  raise `EquinoxRuntimeError: A linear solver returned non-finite output`. See
  [Adjoint conditioning failures](#adjoint-conditioning-failures) below for
  diagnosis and workarounds.

---

## Solver tolerances

The default tolerances match the CVode configuration from the legacy model:

$$
\text{rtol} = 10^{-7}, \qquad
\mathbf{atol} = \bigl[\underbrace{10^{-4},\,10^{-4},\,10^{-4},\,10^{-10},\,10^{-10},\,10^{-4},\,10^{-8}}_{7\text{ bulk vars}},\;\underbrace{10^{-12},\dots,10^{-12}}_{n_r\text{ radii}}\bigr]
$$

The mixed absolute tolerance vector reflects the disparate scales of the state
variables: mixing ratios near $10^{-5}$ kg/kg demand tighter absolute control
than temperature or pressure, and wet radii near $10^{-8}$ m demand even tighter.
The `atol` for radii ($10^{-12}$ m) is deliberately more stringent than the
initial equilibrium tolerance so that the radius trajectory does not accumulate
truncation error during the ascent.

You can override tolerances per-call through the integrator functions in
`pyrcel.integrator`:

```python
from pyrcel.integrator import integrate_parcel, atol_vector
import jax.numpy as jnp

# Tighten radius tolerance for high-concentration runs
custom_atol = jnp.concatenate([
    jnp.array([1e-4, 1e-4, 1e-4, 1e-10, 1e-10, 1e-4, 1e-8]),
    jnp.full(nr, 1e-14),  # tighter radius atol
])
sol = integrate_parcel(y0, args, ts, rtol=1e-8, atol=custom_atol)
```

---

## Step-size control

The v2 backend uses `diffrax.PIDController` with proportional, integral, and
derivative terms acting on the local error estimate [[Soderlind2003]](#Soderlind2003):

$$
h_{n+1} = h_n \cdot \min\!\left(h_\text{max},\, \max\!\left(h_\text{min},\,
\text{safety} \cdot e_n^{-k_I} \cdot e_{n-1}^{k_P} \cdot e_{n-2}^{-k_D}
\right)\right)
$$

where $e_n = \|\mathbf{y}_n^\text{err}\|_\text{rms}$ is the normalised local
error.  The default `diffrax` PID coefficients are $k_I = 0.3/q$, $k_P = 0.4/q$,
$k_D = 0$ (pure PI), where $q = \min(\text{order}_\text{sol},
\text{order}_\text{err}) + 1 = 5$. This is the standard Gustafsson PI controller
widely used in stiff ODE software.

The `dtmax` parameter (forwarded to `PIDController`) can be used to cap the
maximum internal step size, which is occasionally useful for ensuring adequate
density of the trajectory near activation:

```python
sol = integrate_parcel(y0, args, ts, dtmax=0.5)  # cap at 0.5 s
```

---

## Dense output and trajectory storage

Rather than integrating in fixed-size chunks and manually accumulating output
(the v1 `solver_dt`/`output_dt` loop), the v2 backend uses `diffrax`'s
`SaveAt(ts=...)` facility. A single adaptive solve covers the entire time
interval $[0, t_\text{end}]$; the solution is dense-interpolated at the
requested output times using the solver's built-in continuous extension. This
removes the chunking bookkeeping entirely and reduces the number of JAX
compilations to one per unique `(nr, n_output)` configuration.

For the `nd_from_parcel` diagnostic, the backend uses `SaveAt(t1=True)` to
retain only the final state, minimising the memory footprint of the adjoint
checkpoints. `max_supersaturation` uses `SaveAt(ts=ts)` — the same kernel as
`integrate_parcel` — so no second JIT compilation is triggered for gradient
workflows.

The `live=True` mode on `ParcelModel.run()` is the one exception to single-call
integration: it deliberately re-introduces a Python chunk loop to print per-step
diagnostics, at the cost of JIT overhead per chunk and incompatibility with
`jax.grad`. See the [migration guide](migration.md) for the full display-option
comparison (`progress`, `live`, `trajectory_table`).

For interactive runs without early termination (`terminate=False`), $S_\text{max}$
is located post-hoc from the saved trajectory rather than by a second ODE solve.
The interactive layer constructs cubic Hermite polynomials over the two
sub-intervals bracketing the discrete argmax, using `parcel_ode_sys` to supply
the exact endpoint derivatives $\dot{S}$ — the same polynomial family diffrax
stores per solver step when `SaveAt(dense=True)` is used — and finds the
analytic maximum by solving the resulting quadratic $dp/du = 0$. The cost is
three RHS evaluations; the accuracy is O($\Delta t^4$) where $\Delta t$ is the
output spacing.

---

## Supersaturation-maximum event detection

The parcel model's termination criterion is defined by the supersaturation
maximum: the point at which $dS/dt$ changes sign from positive to negative,
marking the end of primary droplet activation. After this point the model is 
integrated a further `terminate_depth` metres before stopping.

Event detection is implemented as a continuous-root-finding solve in
`diffrax.Event`:

```python
def dS_dt_event(t, y, args, **kwargs):
    return parcel_ode_sys(t, y, args)[6]   # the S component of dy/dt

event = dfx.Event(dS_dt_event,
                  root_finder=optx.Newton(rtol=1e-8, atol=1e-12),
                  direction=False)  # trigger only on a downward zero-crossing
```

The `direction=False` flag ensures that only the *downward* zero-crossing
(supersaturation maximum) triggers the event, not the initial upward sweep.
The root is located precisely by Newton iteration to tolerances of $10^{-8}$
(relative) and $10^{-12}$ (absolute), giving $t_\text{smax}$ accurate to
$\lesssim 10^{-6}$ s.

### Why the event-based solve is *not* differentiable

The `diffrax.Event` mechanism makes the stop time $t_\text{smax}$ a
data-dependent function of the initial conditions: JAX cannot differentiate
through a data-dependent control-flow decision. As a result,
`integrate_parcel_terminated` (the interactive path) is **not** suitable for
`jax.grad`. For gradient-based workflows, use the fixed-horizon functions
instead:

```python
# Differentiable: Hermite cubic peak finder (exact via envelope theorem)
from pyrcel.integrator import max_supersaturation
smax = max_supersaturation(y0, args, ts)
grad_smax = jax.grad(max_supersaturation, argnums=(0, 1))(y0, args, ts)

# Differentiable: fixed t_end, soft Nd via sigmoid threshold
from pyrcel.integrator import nd_from_parcel
nd = nd_from_parcel(y0, args, t_end, epsilon=1e-8)
grad_nd = jax.grad(nd_from_parcel, argnums=(0, 1))(y0, args, t_end)
```

The `ts` array passed to `max_supersaturation` must span the supersaturation
peak; see [Differentiable $S_\text{max}$](#differentiable-s_textmax-hermite-cubic-interpolation)
below for guidance on sizing and scaling `ts`.

---

## Differentiable $S_\text{max}$: Hermite cubic interpolation

Computing a gradient-accurate $S_\text{max}$ requires more than returning
`jnp.max(sol.ys[:, 6])`. The adaptive solver takes slightly different internal
step structures for slightly different parameter values (especially updraft
speed $V$). With a discrete output grid, the argmax occasionally jumps by one
time step as $V$ varies, producing an aliased, non-monotone $S_\text{max}(V)$
curve with alternating-sign gradient kinks of $\sim 10^{-3}$ relative
amplitude.

`max_supersaturation` eliminates the aliasing via analytic peak-finding on a
**cubic Hermite interpolant** constructed from the discrete trajectory returned
by `_solve` (the same kernel as `integrate_parcel`). No second dense-output
solve is required.

### Stage 1 — coarse bracket

The argmax of the saved discrete supersaturation values $\{S_j\}$ locates the
approximate peak bin:

$$
i^\ast = \operatorname{argmax}_{j}\, S_j, \qquad j = 0,\dots,n-1
$$

$i^\ast$ is `stop_gradient`'d and used solely to define the two sub-intervals
$[t_{i^\ast - 1},\, t_{i^\ast}]$ and $[t_{i^\ast},\, t_{i^\ast + 1}]$.

### Stage 2 — endpoint derivatives

Three evaluations of the ODE right-hand side at the bracket points supply exact
endpoint slopes:

$$
\dot{S}_{i^\ast - 1},\quad \dot{S}_{i^\ast},\quad \dot{S}_{i^\ast + 1}
$$

where $\dot{S}_j = \mathbf{f}(t_j, \mathbf{y}_j, \boldsymbol{\theta})[6]$.
These are the same Hermite coefficients diffrax stores per step internally when
`SaveAt(dense=True)` is used.

### Stage 3 — analytic quadratic solve

On each sub-interval $[t_0, t_1]$ with $h = t_1 - t_0$, the Hermite cubic in
the normalised coordinate $u = (t - t_0)/h \in [0,1]$ is

$$
p(u) = (2u^3-3u^2+1)\,S_0 + (-2u^3+3u^2)\,S_1
       + h\,(u^3-2u^2+u)\,\dot{S}_0 + h\,(u^3-u^2)\,\dot{S}_1
$$

Setting $dp/du = 0$ yields a quadratic in $u$:

$$
A u^2 + B u + C = 0, \qquad
A = 3\!\left[\tfrac{2(S_0-S_1)}{h} + \dot{S}_0 + \dot{S}_1\right], \quad
B = \tfrac{6(S_1-S_0)}{h} - 2(2\dot{S}_0+\dot{S}_1), \quad
C = \dot{S}_0
$$

The root in $(0,1)$ at which $dp^2/du^2 < 0$ is the interior maximum. If both
sub-intervals have an interior maximum, the larger wins. All root-finding
operates under `jax.lax.stop_gradient` (including the discriminant and the
maximum-of-two selection).

### Why `stop_gradient` on $u^\ast$ is exact

Treating the peak parameter $u^\ast$ as a non-differentiable constant is not an
approximation. By the **envelope theorem**, at a smooth interior maximum

$$
\frac{dS_\text{max}}{d\theta} = \frac{\partial S}{\partial \theta}\bigg|_{u^\ast}
$$

the implicit dependence of $u^\ast$ on $\theta$ drops out identically.
Gradient flows through $p(u^\ast)$ via the endpoint values $S_{i^\ast-1},\,
S_{i^\ast},\, S_{i^\ast+1}$ (ODE adjoint) and the endpoint derivatives
$\dot{S}_{i^\ast-1},\, \dot{S}_{i^\ast},\, \dot{S}_{i^\ast+1}$ (ODE RHS at
those points).

### Output-grid sizing and $t_\text{end}$ scaling

`ts` is used directly for both the ODE output grid and the coarse bracket,
so it must be dense enough to contain the supersaturation peak.

**Rule 1 — ensure the peak lies in the interior.** `ts[-1]` must extend past
$S_\text{max}$ with at least one grid point beyond the peak. For a parcel
ascending from cloud base, the peak occurs at altitude
$z \approx z_0 + V \cdot t_\text{peak}$. A robust horizon choice is

$$
t_\text{end} \approx \frac{1500\,\text{m}}{V},
\qquad \text{clamped to } [200\,\text{s},\, 15000\,\text{s}]
$$

which guarantees the solver reaches $\sim 1500$ m above cloud base regardless
of updraft speed.

**Rule 2 — fix the array shape across $V$ for JIT reuse.** All calls with the
same `ts.shape[0]` reuse the same compiled XLA kernel. Scale only `t_end` with
$V$, keeping the number of output points constant:

```python
_N_TS = 600  # fixed; triggers one JIT compilation

def ts_for_V(V: float) -> jnp.ndarray:
    t_end = max(200.0, min(1500.0 / V, 15000.0))
    return jnp.linspace(0.0, t_end, _N_TS)
```

**Rule 3 — `t_end(V)` does not bias the adjoint.** Because $S_\text{max}$ is
at an interior maximum, the envelope theorem also gives
$\partial S_\text{max}/\partial t_\text{end} = 0$, so the $t_\text{end} \propto
1/V$ coupling introduces no gradient error.

---

## Differentiating through the ODE: `RecursiveCheckpointAdjoint`

Reverse-mode differentiation through an ODE (the adjoint method) requires
storing or recomputing intermediate trajectory values during the backward pass.
`diffrax` implements this via its `adjoint` argument to `diffeqsolve`; the
default in diffrax ≥ 0.4 is `RecursiveCheckpointAdjoint` [[Chen2018]](#Chen2018)
[[Kidger2021]](#Kidger2021).

`RecursiveCheckpointAdjoint` uses a recursive binary tree of checkpoints along
the forward trajectory: the trajectory is divided in half, and each half is
re-solved during the backward pass from the stored midpoint. This gives memory
complexity $O(\log N_\text{steps})$ in contrast to the $O(N_\text{steps})$
cost of storing the full forward trajectory. No configuration is needed — the
default is in force whenever `_solve` is called via
`max_supersaturation` or `nd_from_parcel`.

For the adjoint ODE the backward pass computes $\partial \mathcal{L}/\partial
\mathbf{y}_0$ and $\partial \mathcal{L}/\partial \boldsymbol{\theta}$ by
integrating the adjoint equations backward in time

$$
\frac{d\boldsymbol{\lambda}}{dt} = -\boldsymbol{\lambda}^\top \frac{\partial \mathbf{f}}{\partial \mathbf{y}}, \qquad
\frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \boldsymbol{\theta}} = -\boldsymbol{\lambda}^\top \frac{\partial \mathbf{f}}{\partial \boldsymbol{\theta}}
$$

where $\boldsymbol{\lambda}(T) = \partial \mathcal{L}/\partial \mathbf{y}(T)$
is the terminal condition.  Because `parcel_ode_sys` is fully JAX-traceable,
`diffrax` constructs the Jacobian-vector products $\boldsymbol{\lambda}^\top
(\partial \mathbf{f}/\partial \mathbf{y})$ with `jax.vjp` automatically —
no hand-coded adjoint is required.

### Adjoint conditioning failures

At extreme parameter combinations — particularly high updraft speeds
($V \gtrsim 2$ m/s) paired with very small particles ($\mu \lesssim 0.03\,\mu$m)
— the parcel ODE develops very rapid dynamics near activation. The implicit
solver must then take very stiff Newton steps whose Jacobians can be
nearly singular. `RecursiveCheckpointAdjoint` backpropagates through these
discrete steps, and the linear solve in the adjoint inherits that
ill-conditioning, eventually producing NaN or inf and raising:

```
EquinoxRuntimeError: A linear solver returned non-finite (NaN or inf) output.
```

This is a numerical issue with the discrete adjoint of a stiff implicit solver,
not a bug in the model. The forward pass succeeds at the same point; only the
gradient computation fails.

**Mitigations (roughly in order of invasiveness):**

1. **Tighten forward tolerances**: smaller `rtol`/`atol` force shorter steps,
   keeping each step's Jacobian better-conditioned. Try `rtol=1e-9` and
   radius `atol=1e-14`. The adjoint then inherits a smoother trajectory.

2. **Use `BacksolveAdjoint`**: replaces the discrete adjoint with a
   continuous adjoint ODE solved backwards in time (Pontryagin), which avoids
   differentiating through the Newton iterations entirely.  More numerically
   stable in this regime but can accumulate drift relative to the true gradient
   for very stiff problems.  Pass via `diffeqsolve(... adjoint=BacksolveAdjoint(...))`.

3. **Lower-order implicit solver**: `diffrax.Kvaerno3` or `diffrax.Euler`
   have simpler per-step Jacobian structure and can improve conditioning at
   the cost of accuracy.

4. **Numerical fallback**: for parameter-space sweeps where exact gradients at
   isolated grid points fail, `numpy.gradient` on the pre-computed $S_\text{max}$
   grid provides a reliable finite-difference approximation. The
   [Sensitivity Sweep](../examples/sensitivity_sweep.md) example demonstrates
   this pattern.

---

## MBN2014 activation: gradient via the Implicit Function Theorem

The Morales Betancourt & Nenes (2014) [[MBN2014]](#MBN2014) activation parameterization
finds $S_\text{max}$ as the root of a balance equation

$$
F(S_\text{max},\, \boldsymbol{\theta}) = 0
$$

where the zero-crossing of $F$ depends on the aerosol properties
$\boldsymbol{\theta} = (\mathbf{N}, \boldsymbol{\mu}, \boldsymbol{\sigma},
\boldsymbol{\kappa}, V, T, P, \alpha_c)$. The root is located by bisection,
which is not differentiable by default (the bisection loop contains a
data-dependent branch).

v2 implements the gradient via the **Implicit Function Theorem** (IFT): at the
root $F = 0$, total differentiation gives

$$
\frac{\partial F}{\partial S_\text{max}} \, dS_\text{max}
+ \frac{\partial F}{\partial \boldsymbol{\theta}} \cdot d\boldsymbol{\theta} = 0,
\qquad \Longrightarrow \qquad
\frac{\partial S_\text{max}}{\partial \boldsymbol{\theta}}
= -\frac{\partial F / \partial \boldsymbol{\theta}}{\partial F / \partial S_\text{max}}
$$

This is implemented using `jax.custom_vjp` in `pyrcel.activation._mbn2014`. The
forward pass runs the bisection to find $S_\text{max}$; the custom backward
pass evaluates the IFT formula above using `jax.grad(F)`, which is
$O(1)$ in memory and independent of the number of bisection iterations.

The two partial derivatives needed are

$$
\frac{\partial F}{\partial S_\text{max}} = \frac{dF}{dS_\text{max}}\bigg|_{\boldsymbol{\theta}\,\text{fixed}}, \qquad
\frac{\partial F}{\partial \boldsymbol{\theta}} = \nabla_{\boldsymbol{\theta}} F\big|_{S_\text{max}\,\text{fixed}}
$$

both evaluated at the converged root. Because $F$ is an analytic function of
all its arguments, both partial derivatives are computed exactly by JAX's
forward-mode AD (`jax.jacfwd`).

---

## Differentiable activated droplet number: sigmoid soft threshold

The hard-threshold Nd diagnostic (`model.run(mode="nd")`) counts bins where the
final wet radius exceeds the Köhler critical radius:

$$
N_d^{\text{hard}} = \sum_i N_i \cdot \mathbf{1}\!\left[r_i(t_\text{end}) \geq r_{\text{crit},i}\right]
$$

The Heaviside indicator $\mathbf{1}[\cdot]$ has zero gradient almost everywhere,
making this non-differentiable. The `nd_from_parcel` function replaces the
indicator with a logistic sigmoid:

$$
N_d^{\text{soft}} = \sum_i N_i \cdot \sigma\!\left(\frac{r_i(t_\text{end}) - r_{\text{crit},i}}{\varepsilon}\right),
\qquad
\sigma(x) = \frac{1}{1+e^{-x}}
$$

where $\varepsilon$ (default $10^{-8}$ m $\approx 10$ nm) controls the sharpness
of the threshold. In the limit $\varepsilon \to 0$, $N_d^\text{soft} \to
N_d^\text{hard}$ (using the approximate Köhler formula); in practice
$|N_d^\text{soft} - N_d^\text{hard}|/N_d^\text{hard} \lesssim 10\%$
at $\varepsilon = 10^{-8}$ m, with the residual discrepancy arising from the
approximate vs. exact Köhler critical radius.

The critical radius is computed with the approximate Köhler formula
[[Petters2007]](#Petters2007):

$$
r_{\text{crit},i} = \sqrt{\frac{3 \kappa_i r_{d,i}^3}{A}},
\qquad
A = \frac{2 M_w \sigma_w(T)}{R T \rho_w}
$$

This form is JAX-traceable and vectorises over all bins in a single operation,
unlike the exact `fminbound` solve used by `binned_activation`.

Gradients of $N_d^\text{soft}$ with respect to $(\mathbf{y}_0,
\boldsymbol{\theta})$ flow through both the ODE trajectory (via
`RecursiveCheckpointAdjoint`) and the Köhler computation:

```python
from pyrcel.integrator import nd_from_parcel

grad_fn = jax.grad(nd_from_parcel, argnums=(0, 1))
grad_y0, grad_args = grad_fn(y0, args, t_end, epsilon=1e-8)
dNd_dV = grad_args[4].V   # gradient w.r.t. updraft speed (ConstantV module)
```

---

## Low-level integrator API

`ParcelModel.run()` is the recommended entry point, but the underlying functions
are public and fully composable for advanced workflows — building a custom
optimization loop, wrapping the ODE in an external framework, or constructing
partial pipelines with `jax.vmap`.

### Equilibration

```python
from pyrcel.equilibrate import equilibrate_initial_state

y0 = equilibrate_initial_state(T0, S0, P0, r_drys, kappas, Nis)
# y0: shape (7 + nr,) — [z, P, T, wv, wc, wi, S, r_0, ..., r_{nr-1}]
```

`equilibrate_initial_state` computes the equilibrium wet-radius spectrum at the
initial conditions and assembles the full `(7 + nr)`-element state vector. It is
the direct equivalent of `ParcelModel._setup_run` and is itself JAX-traceable
(uses `optimistix` bisection).

### Core solve

| Function | Returns | When to use |
|---|---|---|
| `integrate_parcel(y0, args, ts)` | `diffrax.Solution` | Full solution object with `sol.ts`, `sol.ys`, `sol.result` |
| `integrate_parcel_arrays(y0, args, ts)` | `(ts, ys, success)` | Raw arrays; avoids unpacking the diffrax object |

```python
from pyrcel.integrator import integrate_parcel_arrays

ts_out = jnp.linspace(0.0, t_end, n_out)
ts, ys, ok = integrate_parcel_arrays(y0, args, ts_out)
# ys: shape (n_out, 7 + nr)
```

`args` is the tuple `(r_drys, Nis, kappas, accom, V)` where `V` is a
`ConstantV` or `InterpolatedUpdraft`.

### Differentiable diagnostics

| Function | Returns | Notes |
|---|---|---|
| `max_supersaturation(y0, args, ts)` | `float` | $S_\text{max}$ via Hermite cubic; `jax.grad`-able |
| `nd_from_parcel(y0, args, t_end)` | `float` | Soft-threshold $N_d$ (m⁻³); `jax.grad`-able |

Both are designed for the differentiable path and documented in detail in the
[adjoint](#differentiating-through-the-ode-recursivecheckpointadjoint) and
[sigmoid Nd](#differentiable-activated-droplet-number-sigmoid-soft-threshold)
sections above.

### Terminated-run pipeline

When `terminate=True` the interactive path uses a three-function pipeline:

```python
from pyrcel.integrator import (
    find_smax,
    terminate_cutoff_time,
    integrate_parcel_terminated,
)

# 1. Precisely localize S_max via a dS/dt = 0 event
t_smax, smax, y_smax, activated = find_smax(y0, args, t_end)

# 2. Compute the altitude-based cutoff time (terminate_depth m past S_max)
t_cut = terminate_cutoff_time(
    y0, args, t_end,
    terminate_depth=10.0, t_smax=t_smax, y_smax=y_smax, activated=activated,
)

# 3. Integrate to the cutoff and return the full trajectory
ts, ys = integrate_parcel_terminated(y0, args, t_cut, output_dt=1.0)
```

`find_smax` is **not** differentiable (event detection is discontinuous; see
[the event section](#supersaturation-maximum-event-detection)). Use
`max_supersaturation` on the differentiable path.

### Tolerance constants

The module-level constants `STATE_RTOL`, `STATE_ATOL`, and `RADIUS_ATOL` are
importable and match the CVode configuration from the legacy model:

```python
from pyrcel.integrator import STATE_RTOL, atol_vector

atol = atol_vector(nr)  # per-component vector matching legacy CVode setup
```

---

## Performance tips

### JIT compilation warm-up

The first call to any `jax.jit`-compiled function traces and compiles the
computation graph to XLA HLO. For a 200-bin parcel simulation this takes
$\sim 10$–$30$ s on a modern CPU. Subsequent calls with the same shapes are
fast (typically $< 1$ s). To amortise this cost:

- Call the function once with representative inputs before timing production runs.
- Use `jax.block_until_ready` to ensure the first call completes before starting
  a timer.

```python
import jax
from pyrcel.integrator import max_supersaturation

# Warm up the JIT cache
_ = jax.block_until_ready(max_supersaturation(y0, args, ts))

# Now time the production call
import time
t0 = time.perf_counter()
smax = float(jax.block_until_ready(max_supersaturation(y0, args, ts)))
print(f"S_max = {smax:.5f}  [{time.perf_counter()-t0:.3f}s]")
```

### `jax.vmap` for ensemble runs

For a batch of $n$ parcel simulations that share the same ODE structure but
differ in initial conditions or parameters (e.g., an updraft-speed ensemble),
`jax.vmap` compiles a single batched kernel that runs all $n$ members in parallel:

```python
import jax

# Batch over updraft speeds
V_samples = jnp.linspace(0.1, 3.0, 128)

def single_run(V):
    new_args = (r_drys, Nis, kappas, accom, pm.ConstantV(V))
    return max_supersaturation(y0, new_args, ts)

smax_batch = jax.jit(jax.vmap(single_run))(V_samples)  # shape (128,)
```

This is substantially faster than a Python `for` loop: the XLA compiler can fuse
the batch dimension with the ODE internal loop and exploit SIMD or GPU
parallelism. The `run_updraft_ensemble` convenience function does exactly this:

```python
result = pm.run_updraft_ensemble(
    [aerosol], T0=283.15, S0=0.0, P0=85000.0,
    mean=0.5, std=0.2, n=1024
)
```

### Gradient computation and `jax.jacfwd`

For scalar outputs (`max_supersaturation`, `nd_from_parcel`), `jax.grad` is the
most efficient gradient computation. For vector outputs or Jacobian matrices,
`jax.jacfwd` (forward-mode) is usually faster when the number of inputs is small
relative to outputs, while `jax.jacrev` (reverse-mode) is faster in the opposite
regime:

```python
# Gradient of S_max w.r.t. a vector of N_i (nr inputs, 1 output → jacrev)
grad_N = jax.jacrev(lambda N: max_supersaturation(
    y0, (r_drys, N, kappas, accom, V), ts
))(Nis)
```

### GPU dispatch

Pass `device="gpu"` to `ParcelModel` or wrap integrator calls in a
`jax.default_device` context:

```python
import jax

with jax.default_device(jax.devices("gpu")[0]):
    smax = max_supersaturation(y0, args, ts)
```

For ensemble runs the GPU backend provides the largest speedups: a 1024-member
ensemble that takes $\sim 60$ s on 8 CPU cores typically completes in $< 5$ s on
a single A100.

---

## Numerical accuracy

The v2 solver at default tolerances agrees with the v1 CVode baseline to within:

| Quantity | Typical discrepancy |
|---|---|
| $S_\text{max}$ | $< 10^{-5}$ (absolute) |
| Activated fraction | $\lesssim 0.1\%$ (relative) |
| Droplet number $N_d$ | $\lesssim 1\%$ (relative) |
| Trajectory $r_i(t)$ | $< 10^{-10}$ m (absolute) |

These values are confirmed by the golden-data regression test suite in
`tests/test_golden.py` and `tests/test_integration.py`, which compare v2 output
against frozen v1 reference trajectories across a range of aerosol configurations
and updraft speeds.

---

## References

<a id="Nenes2001"></a>**[[Nenes2001]](#Nenes2001)** Nenes, A., Ghan, S., Abdul-Razzak, H., Chuang, P. Y., & Seinfeld, J. H.
(2001). Kinetic limitations on cloud droplet formation and impact on cloud albedo.
*Tellus B*, **53**(2), 133–149. <https://doi.org/10.1034/j.1600-0889.2001.d01-12.x>

<a id="Ghan2011"></a>**[[Ghan2011]](#Ghan2011)** Ghan, S. J., et al. (2011). Droplet nucleation: Physically-based
parameterizations and comparative evaluation. *J. Adv. Model. Earth Syst.*, **3**,
M10001. <https://doi.org/10.1029/2011MS000074>

<a id="Seinfeld2006"></a>**[[Seinfeld2006]](#Seinfeld2006)** Seinfeld, J. H., & Pandis, S. N. (2006). *Atmospheric Chemistry and
Physics: From Air Pollution to Climate Change* (2nd ed.). Wiley.

<a id="Petters2007"></a>**[[Petters2007]](#Petters2007)** Petters, M. D., & Kreidenweis, S. M. (2007). A single parameter
representation of hygroscopic growth and cloud condensation nucleus activity.
*Atmos. Chem. Phys.*, **7**, 1961–1971. <https://doi.org/10.5194/acp-7-1961-2007>

<a id="Kvaerno2004"></a>**[[Kvaerno2004]](#Kvaerno2004)** Kvaernø, A. (2004). Singly diagonally implicit Runge–Kutta methods with
an explicit first stage. *BIT Numerical Mathematics*, **44**(3), 489–502.
<https://doi.org/10.1023/B:BITN.0000046811.70614.38>

<a id="Soderlind2003"></a>**[[Soderlind2003]](#Soderlind2003)** Söderlind, G. (2003). Digital filters in adaptive time-stepping.
*ACM Trans. Math. Softw.*, **29**(1), 1–26. <https://doi.org/10.1145/641876.641877>

<a id="Chen2018"></a>**[[Chen2018]](#Chen2018)** Chen, R. T. Q., Rubanova, Y., Bettencourt, J., & Duvenaud, D. (2018).
Neural ordinary differential equations. *NeurIPS*, 6571–6583.
<https://arxiv.org/abs/1806.07366>

<a id="Kidger2021"></a>**[[Kidger2021]](#Kidger2021)** Kidger, P. (2021). *On Neural Differential Equations* (Ph.D. thesis).
University of Oxford. <https://arxiv.org/abs/2202.02435>

<a id="MBN2014"></a>**[[MBN2014]](#MBN2014)** Morales Betancourt, R., & Nenes, A. (2014). Droplet activation
parameterization: the population-splitting concept revisited. *Geosci. Model Dev.*,
**7**, 2345–2357. <https://doi.org/10.5194/gmd-7-2345-2014>

<a id="Rothenberg2016"></a>**[[Rothenberg2016]](#Rothenberg2016)** Rothenberg, D., & Wang, C. (2016). Metamodeling of droplet
activation for global climate models. *J. Atmos. Sci.*, **73**(4), 1255–1272.
<https://doi.org/10.1175/JAS-D-15-0223.1>
