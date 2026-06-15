# Design Doc: Migrating `pyrcel` from numba/Assimulo to the JAX + `diffrax` ecosystem

- **Status:** Draft / for discussion
- **Author:** (design sketch)
- **Scope:** A **clean v2 rewrite**. Replace the `numba`-compiled right-hand-side (RHS)
  and the Assimulo/SUNDIALS `CVode` integrator with a JAX-native implementation built on
  `diffrax`, preserving the *physics and numerical results* of the current model but
  **not** its internal code paths or APIs where they get in the way.
- **Decisions locked in (see §9):**
  - This is **v2**: drop `numba` and Assimulo/`CVode` **outright**. No backwards-compat
    shims, no dual-backend dispatch, no legacy code paths.
  - Validation is **against the `master` branch** (current numba + CVode), used purely as
    an external oracle to generate frozen golden data — not as a runtime dependency of v2.
  - Use the **JAX ecosystem end-to-end**, including `optimistix` for the equilibration
    root-finds, so all numerics live in one ecosystem.
  - Differentiability requirement is scoped: we must support **gradients through the
    integration with respect to the (already-equilibrated) initial conditions** — i.e.
    `d(output)/d(y0)` and `d(output)/d(physical params)`. We do **not** require gradients
    to flow *through* the equilibration solve itself.
- **Non-goals:** Adding new physics (ice/freezing, time-varying `N`, etc.). Those are
  tracked separately; this migration should be physics-preserving.

---

## 1. Motivation

The current model has two compiled/native dependencies that create friction:

1. **`numba` AOT/JIT for the derivative.** `pyrcel/_parcel_aux_numba.py` uses
   `numba.njit` plus `numba.pycc.CC` (ahead-of-time compilation) and a parallel
   `prange` loop over aerosol bins. `numba.pycc` is **deprecated** and slated for
   removal, and the AOT path duplicates the thermodynamics that already live in
   `pyrcel/thermo.py` (drift risk — the two copies must be kept in sync by hand).
2. **Assimulo / SUNDIALS `CVode`.** This is the only supported solver
   (see `README.md`: "As of version 1.3, we no longer support any ODE solver backends
   other than `cvode`"). It is `conda`-only, hard to `pip install`, and the README
   itself flags the SUNDIALS dependency as the hardest part of installation.

Moving to JAX + `diffrax` would:

- **Simplify installation** — pure-Python `pip`/`conda` install, no SUNDIALS build.
- **Unblock the project's stated roadmap** — the README already advertises a planned
  "suite of ODE solvers from the `diffrax` toolkit," and there is prior art on the
  `jax_jacobian` branch ("Initial testing with autodifferentiation").
- **Enable differentiability.** This is the big one. The model author's research is
  literally *"Metamodeling of Droplet Activation for Global Climate Models"* (Rothenberg
  & Wang, 2016). With `jax.grad`/`jax.jacfwd` we can get exact gradients of, e.g.,
  `S_max` or activated fraction with respect to `V`, `T0`, `P0`, `kappa`, distribution
  parameters, and the accommodation coefficient `accom` — directly useful for sensitivity
  analysis, parameter inversion, and training metamodels/emulators.
- **Enable batched ensembles.** `jax.vmap` + `jit` (and GPU/TPU) make the existing Monte
  Carlo workflow (`monte_carlo_outputs.csv`) embarrassingly parallel without a Python
  loop over `ParcelModel` instances.
- **Single source of truth** for thermodynamics (one JAX implementation instead of
  numba-copy + numpy-copy).

This is not free: JAX has real gotchas (float32 default, bounded `max_steps` under `jit`,
compile latency) and the stiff solver is a *different* algorithm than CVode's BDF, so
results will not be bit-identical. That is precisely why this doc pairs the migration
with a regression harness (§7).

---

## 2. Current model structure (as-is)

### 2.1 State vector

The integrated state `y` has length `N_STATE_VARS + nr` where `N_STATE_VARS = 7`
(`pyrcel/constants.py`) and `nr` is the total number of aerosol size bins across all
species:

| index            | symbol | meaning                              | units   |
| ---------------- | ------ | ------------------------------------ | ------- |
| `y[0]`           | `z`    | altitude                             | m       |
| `y[1]`           | `P`    | pressure                             | Pa      |
| `y[2]`           | `T`    | temperature                          | K       |
| `y[3]`           | `wv`   | water-vapor mass mixing ratio        | kg/kg   |
| `y[4]`           | `wc`   | cloud liquid water mass mixing ratio | kg/kg   |
| `y[5]`           | `wi`   | cloud ice water mass mixing ratio    | kg/kg   |
| `y[6]`           | `S`    | supersaturation (0 = 100% RH)        | –       |
| `y[7:7+nr]`      | `r_i`  | wet radius of each aerosol bin       | m       |

### 2.2 RHS (`parcel_ode_sys`)

Defined in `pyrcel/_parcel_aux_numba.py`. Per call it:

1. Computes thermodynamic quantities (`es`, virtual temp `Tv`, densities `rho_air`,
   `rho_air_dry`, vapor pressure).
2. Loops over the `nr` bins (parallel `prange`), computing for each:
   - non-continuum diffusivity `dv(T, r, P, accom)` and conductivity `ka(T, r, rho)`,
   - the growth coefficient `G = 1/(G_a + G_b)`,
   - equilibrium supersaturation `Seq(r, r_dry, T, kappa)` (κ-Köhler),
   - droplet growth `dr/dt = (G/r)(S - Seq)` and its contribution to `dwc/dt`.
3. Forms the bulk tendencies: `dP/dt` (hydrostatic), `dwv/dt = -(dwc/dt + dwi/dt)`,
   `dT/dt` (adiabatic + latent), `dS/dt` via the Ghan (2011) form
   `dS/dt = α V − γ dwc/dt`, and `dz/dt = V`.

The per-bin loop is **trivially vectorizable** — there are no cross-bin dependencies
except the `dwc/dt` reduction (a sum). This maps cleanly onto array ops in JAX.

### 2.3 Equilibration / setup (`ParcelModel._setup_run`, `pyrcel/parcel.py`)

- Concatenates all species' `r_drys`, `kappas`, `Nis` into flat arrays.
- For each bin, finds the equilibrium wet radius `r0` by **`scipy.optimize.bisect`** on
  `Seq(r, r_dry, T0, kappa) - S0`, bracketed by `r_dry` and the Köhler critical radius
  `r_crit` (from `kohler_crit`, which uses **`scipy.optimize.fminbound`**).
- Builds `y0` and computes initial `wc0`, `wv0`.

### 2.4 Integration (`pyrcel/integrator.py`)

- `CVODEIntegrator` wraps Assimulo `CVode` (BDF, `maxord=5`, Newton iteration).
- Tolerances: `rtol = 1e-7`; `atol = [1e-4,1e-4,1e-4,1e-10,1e-10,1e-4,1e-8]` for the 7
  state vars plus `1e-12` for every radius.
- **Chunked loop**: integrates in `solver_dt` chunks, interpolating `n_out` points per
  chunk for output at `output_dt` cadence.
- **Event/termination**: `ExtendedProblem` watches the sign of `dS/dt`; once the
  supersaturation maximum is detected, it sets a cutoff time `t_smax + terminate_depth/V`
  and stops storing results past it (i.e. it integrates an extra `terminate_depth`
  meters after `S_max`).
- Output array `x` is converted to DataFrames/NetCDF in `pyrcel/output.py`.

### 2.5 Time-varying updraft

`V` may be a callable `V(t)` (`parcel.py` re-wraps the RHS, `hasattr(self.V, "__call__")`).
The JAX port must preserve this (cleanly, via a callable/interpolant passed as `args`).

---

## 3. Target architecture

```
                       +-------------------------------+
   AerosolSpecies  --> |  ParcelModel (v2)             |
   (setup/equilibrate) |  - set_initial_conditions     |
                       |  - run(...)   [diffrax only]   |
                       +---------------+---------------+
                                       |
        +------------------------------+------------------------------+
        |                              |                              |
+-------v--------+        +------------v-----------+      +-----------v----------+
| thermo_jax.py  |        | parcel_aux_jax.py      |      | integrator_diffrax.py|
| es, dv, ka,    | -----> | parcel_ode_sys (jit'd, |----> | diffeqsolve(         |
| sigma_w, Seq   |        | fully vectorized RHS)  |      |   Kvaerno5, PID,     |
| (jnp)          |        |                        |      |   Event, SaveAt)     |
+-------^--------+        +------------------------+      +----------------------+
        |                  +---------------------+
        +----------------- | equilibrate (optx)  | -- y0 -->  (initial conditions)
                           +---------------------+
```

Guiding principles:

- **One backend, no dispatch.** v2 integrates with `diffrax` only. There is no
  `solver=` switch and no Assimulo/numba code. The `master` branch stays around purely as
  an *offline oracle* to bake golden reference data (§7); v2 never imports it.
- **One thermodynamics implementation.** `thermo_jax.py` (in `jnp`) is the single source
  of truth; the numba copy in `_parcel_aux_numba.py` is **deleted**. Any numpy-facing
  helper is a thin wrapper over the same expressions (with `jax_enable_x64`), not a second
  hand-maintained copy.
- **Functional core.** The RHS is a pure function `f(t, y, args) -> dy/dt`. Static aerosol
  data (`r_drys`, `Nis`, `kappas`, `accom`, `V`) is carried in `args` — either as a plain
  pytree (tuple/dict) or an `equinox.Module` (decision pending, see §5.4).
- **Equilibration in `optimistix`.** Initial wet radii are solved with `optimistix` (not
  `scipy`) to keep all numerics in one ecosystem (§4.3). Gradients are required through
  the *integration* w.r.t. the resulting `y0`, but not through the equilibration solve.

---

## 4. Detailed migration plan, component by component

### 4.1 Thermodynamics → `thermo_jax.py`

Direct `np.` → `jnp.` translation of `sigma_w`, `ka`, `dv`, `es`, `Seq`. These are
elementwise and trivially differentiable. Notes:

- `Seq` contains `r**3 - r_dry**3` in both numerator and denominator — **catastrophic
  cancellation** when `r ≈ r_dry` (freshly equilibrated small particles). This is exactly
  why **float64 is mandatory** (see §6.1).
- Keep the numpy `thermo.py` functions; the JAX versions must match them to `~1e-12`
  (unit test, §7.1) reusing the existing `generate_data.py` reference fixture pattern.

### 4.2 RHS → `parcel_aux_jax.py`

Replace the `prange` loop with vectorized array ops:

```python
import jax.numpy as jnp
from functools import partial
import jax

@jax.jit  # nr inferred from array shape, not a static arg
def parcel_ode_sys(t, y, args):
    r_drys, Nis, kappas, accom, V = args   # V: scalar, or a callable/interpolant
    V_t = V(t) if callable(V) else V       # (with Equinox, V is a Module: V_t = V(t))
    P, T, S = y[1], y[2], y[6]
    wv = y[3]
    rs = y[7:]

    pv_sat = es(T - 273.15)
    Tv = (1.0 + 0.61 * wv) * T
    e = (1.0 + S) * pv_sat
    rho_air = P / c.Rd / Tv
    rho_air_dry = (P - e) / c.Rd / T

    dv_r = dv(T, rs, P, accom)          # vectorized over bins
    ka_r = ka(T, rs, rho_air)
    G_a = (c.rho_w * c.R * T) / (pv_sat * dv_r * c.Mw)
    G_b = (c.L * c.rho_w * ((c.L * c.Mw / (c.R * T)) - 1.0)) / (ka_r * T)
    G = 1.0 / (G_a + G_b)
    delta_S = S - Seq(rs, r_drys, T, kappas)
    drs_dt = (G / rs) * delta_S

    dwc_dt = 4.0 * jnp.pi * c.rho_w / rho_air_dry * jnp.sum(Nis * rs * rs * drs_dt)
    dwv_dt = -dwc_dt
    dP_dt = -rho_air * c.g * V_t
    dT_dt = -c.g * V_t / c.Cp - c.L * dwv_dt / c.Cp
    alpha = (c.g * c.Mw * c.L) / (c.Cp * c.R * T**2) - (c.g * c.Ma) / (c.R * T)
    gamma = (P * c.Ma) / (c.Mw * pv_sat) + (c.Mw * c.L * c.L) / (c.Cp * c.R * T * T)
    dS_dt = alpha * V_t - gamma * dwc_dt
    return jnp.concatenate([jnp.array([V_t, dP_dt, dT_dt, dwv_dt, dwc_dt, 0.0, dS_dt]), drs_dt])
```

**Subtle fidelity point — reduction order.** The numba version accumulates `dwc_dt`
in a `prange` loop (nondeterministic/parallel order); JAX uses a tree reduction
(`jnp.sum`). Both are mathematically the same sum but differ at the ULP level. This is
within tolerance but worth noting for the harness.

### 4.3 Equilibration / setup → `optimistix`

`_setup_run` finds each bin's equilibrium wet radius `r0` (root of
`Seq(r, r_dry, T0, kappa) - S0`) bracketed between `r_dry` and the Köhler critical radius
`r_crit`. v2 replaces the `scipy` calls with `optimistix`:

- `scipy.optimize.bisect` → **`optimistix.Bisection`** (robust, bracketing) or
  `optimistix.Newton` (faster, `Seq` is smooth and we have its derivative for free via
  `jax`). Bisection mirrors the current method most closely and is the safer first cut.
- `kohler_crit`'s `scipy.optimize.fminbound` → either the **closed-form approximate
  critical point** already implemented in `thermo.kohler_crit(approx=True)`, or an
  `optimistix` minimiser / a root-find on `d Seq/dr = 0` (the derivative is available
  analytically through `jax.grad`).

The whole per-bin equilibration is `vmap`-friendly (one independent root-find per bin).

**Differentiability scope.** Per the locked decision, we need `d(output)/d(y0)` through the
*integration*, not gradients *through* the equilibration solve. So equilibration can run
eagerly (even outside `jit`) and simply hand `y0` to the integrator as the differentiation
input. (As a free bonus, `optimistix` root-finds *are* differentiable via the implicit
function theorem, so we can opt into end-to-end gradients later without re-architecting.)

**Fidelity note for the harness.** Because we are also changing the root-finder, `y0` may
differ from `master` at the `~1e-12`–`1e-10` level. The harness (§7) should either (a)
assert the v2 `y0` matches `master`'s `y0` to a tight tolerance as its own test, or
(b) seed the v2 integration with `master`'s exact `y0` when isolating integrator-only
differences. Recommend doing both: a dedicated equilibration-equivalence test, plus
integrator runs seeded from a common `y0`.

### 4.4 Integrator → `integrator_diffrax.py` (`DiffraxIntegrator(Integrator)`)

Implements the existing `Integrator` ABC interface (`integrate(t_end) -> (x, t, success)`)
so it slots into `ParcelModel.run` via `Integrator.solver("diffrax")`.

```python
import diffrax as dfx
import optimistix as optx

term = dfx.ODETerm(parcel_ode_sys)
root_finder = dfx.with_stepsize_controller_tols(dfx.VeryChord)()   # inherits rtol/atol
solver = dfx.Kvaerno5(root_finder=root_finder)                     # A-L stable, stiffly accurate
controller = dfx.PIDController(
    rtol=1e-7,
    atol=jnp.array([1e-4,1e-4,1e-4,1e-10,1e-10,1e-4,1e-8] + [1e-12]*nr),
)
saveat = dfx.SaveAt(ts=output_times)        # dense output at output_dt cadence
sol = dfx.diffeqsolve(term, solver, t0, t_end, dt0=None, y0=y0, args=args,
                      stepsize_controller=controller, saveat=saveat,
                      event=event, max_steps=...)
```

Mapping the current solver knobs:

| Current (CVode)                 | `diffrax` equivalent                                       |
| ------------------------------- | ---------------------------------------------------------- |
| BDF, `maxord=5`, Newton         | `Kvaerno5` (ESDIRK) + Newton/`VeryChord` root finder       |
| `rtol`/per-state `atol`         | `PIDController(rtol=..., atol=<vector>)`                    |
| `maxh = min(0.1, output_dt)`    | `PIDController(dtmax=...)`                                  |
| `solver_dt` chunking for output | `SaveAt(ts=...)` with dense interpolation (no manual loop) |
| `max_steps`                     | `diffeqsolve(max_steps=...)` (**must be finite for jit**)  |

The `solver_dt`/`output_dt` chunking loop largely **disappears**: `diffrax` does adaptive
stepping internally and `SaveAt(ts=...)` gives output at any requested cadence via dense
interpolation. (We may still chunk later if/when ice time-splitting lands, per the comment
in `pm_test.py`.)

### 4.5 Termination at `S_max` (events) — the trickiest mapping

The current behavior: detect the supersaturation maximum (`dS/dt` changes sign +→−), then
integrate an **extra** `terminate_depth` meters and stop.

`diffrax` events terminate *at* the trigger, not after a delay. Options:

- **(A) Two-phase solve (most faithful).**
  1. Solve with `dfx.Event(cond_fn=lambda t, y, args, **kw: rhs(t, y, args)[6])` — a
     *continuous* condition returning `dS/dt`; with an `optimistix.Newton` root finder
     `diffrax` finds the exact `t_smax`.
  2. Resume a second bounded solve from that state for `terminate_depth/V` more seconds,
     with its own `SaveAt`.
- **(B) Single fixed-horizon solve + post-trim (simplest, great for the harness).**
  Because `z = ∫V dt` is monotonic, just integrate to `t_end`, save at `output_dt`, and
  trim the trajectory in post-processing to `z ≤ z_smax + terminate_depth`. No event
  needed. Best for first validation since it removes event-localization as a variable.
- **(C) `diffrax.steady_state_event()`** detects `‖f‖` small — not the same as the
  supersaturation *max*, so not a direct substitute. Mention but don't use for this.

**Recommendation:** start with (B) for validation, then implement (A) for production
parity and performance (early stop).

### 4.6 Output

`diffrax` returns `sol.ts`, `sol.ys`. Adapt to the existing `(x, t, success)` contract
expected by `ParcelModel.run`/`output.parcel_to_dataframes` — `x = np.asarray(sol.ys)`,
`t = np.asarray(sol.ts)`, `success = (sol.result == dfx.RESULTS.successful)`. No changes
needed downstream in `output.py`.

### 4.7 Integration structure: one adaptive solve vs. chunking

The `master` model integrates in `solver_dt` chunks (a Python `while` loop calling
Assimulo repeatedly). Three things motivated that design: (1) producing output at an
`output_dt` cadence, (2) live terminal feedback (the formatted `step | time | z T S`
table), and (3) a CVode/BDF-specific worry that a multistep method "predicts too far ahead"
on a stiff problem (see the comment in `pm_test.py`). It's worth re-deriving whether we
still want chunking in v2, because it interacts with autodiff and performance in
non-obvious ways.

**Key realization: in `diffrax`, the three original motivations are mostly handled
natively, without chunking.**

- **Output cadence** → `SaveAt(ts=output_times)` gives dense-interpolated output at any
  cadence from a *single* solve.
- **Live feedback** → `diffeqsolve(..., progress_meter=dfx.TextProgressMeter())` (or
  `dfx.TqdmProgressMeter()`) renders a progress bar/percentage during one big solve via
  host callbacks. (Confirmed present in `diffrax` 0.6.2.)
- **Stiffness "look-ahead"** → that was a property of multistep BDF. `Kvaerno5` is a
  *single-step* ESDIRK with an adaptive `PIDController`; it does not extrapolate across a
  long horizon, so the original rationale for capping the horizon does not carry over.

#### Trade-offs if we chunk anyway

| Concern | One big `diffeqsolve` | Python chunk loop (un-`jit`'d boundaries) | `lax.scan` chunks (compiled) |
| --- | --- | --- | --- |
| **Gradients (correctness)** | ✅ exact, via `RecursiveCheckpointAdjoint` | ✅ exact (backprop chains through chunks) | ✅ exact |
| **Gradient memory/perf** | ✅ optimal `O(log n_steps)` checkpointing built in | ⚠️ if the whole loop is wrapped in one `jit`/`grad`, the Python loop **unrolls** → larger graph, longer compile; chunk boundaries act as (suboptimal) manual checkpoints | ➖ `scan` keeps it rolled, but you now hand-manage checkpointing and adaptive state |
| **Early termination (`terminate=True`)** | ✅ native `Event` (data-dependent stop stays inside XLA) | ✅ trivial (check in Python between chunks) | ⚠️ awkward (data-dependent length needs `while_loop`) |
| **Per-chunk Python feedback (custom table)** | ❌ not possible mid-solve (must use host-callback progress meter instead) | ✅ natural — print between chunks | ❌ no Python between iterations |
| **`vmap` over an ensemble** | ✅ clean; host-callback meters get messy so disable them | ❌ Python loop doesn't `vmap`; host I/O per element is chaotic | ✅ clean |
| **Adaptive step continuity** | ✅ never reset | ⚠️ must thread `solver_state`/`controller_state`/`made_jump` across calls or the controller re-initialises each chunk (wasted steps at every boundary) | ⚠️ same threading needed |
| **Compile count** | 1 | 1 (same compiled chunk fn reused) provided shapes are constant | 1 |

**The core tension:** you cannot simultaneously (a) drop into Python between chunks to
print a custom table *and* (b) keep the entire multi-chunk integration inside a single XLA
graph for `jit`/`grad`/`vmap`. Per-chunk Python feedback forces the integration out of one
compiled unit, which is exactly the form you *don't* want for differentiable/batched runs.

#### Recommendation: decouple feedback from integration; default to a single solve

Adopt **two presentation modes over the same numerical core**:

1. **Core / differentiable / batched mode (default for `grad`, `vmap`, library use).**
   A single `jit`-ed `diffeqsolve` with `SaveAt(ts=...)`, `Event`-based early termination,
   `RecursiveCheckpointAdjoint`, and `progress_meter=NoProgressMeter()`. This is the
   cleanest for performance and autodiff and is what the differentiability requirement
   (gradients through the integration w.r.t. `y0`) needs.

2. **Interactive / CLI mode (`console=True`, `run_parcel`).** Same single solve, but with
   `progress_meter=TextProgressMeter()`/`TqdmProgressMeter()` for the live bar, and a
   nicely formatted **summary table computed from the dense output after the solve**
   (`S_max`, `t_smax`, activated fractions, etc.). This reproduces the *spirit* of the old
   console output without a Python chunk loop.

If a *live, per-interval* stats table (updating every ~10 s of sim time mid-run) turns out
to be a hard product requirement, implement it as a thin **interactive-only** chunk loop
that threads `solver_state`/`controller_state` across `diffeqsolve` calls — and explicitly
keep it **out of** the `jit`/`grad`/`vmap` path. In other words, chunking becomes a
*pure-UI affordance*, never part of the numerical/differentiable contract.

**Note on future physics.** When ice + operator/time-splitting lands (flagged in
`pm_test.py`), that introduces *physics-driven* sub-stepping — best expressed as a
`lax.scan` over splitting steps inside the compiled core. That is a separate concern from
the *output/feedback* chunking discussed here and does not change this recommendation.

---

## 5. Critical assessment of the JAX tool ecosystem

The user asked specifically: which JAX tools to use *alongside* `diffrax`, or *as a
replacement* for its core capabilities. Assessment:

### 5.1 Recommended companions to `diffrax`

> **Dependency cost is ~zero.** `diffrax==0.6.2` already declares `equinox>=0.11.10`,
> `optimistix>=0.0.7`, `lineax>=0.0.5`, and `jaxtyping>=0.2.24` as hard dependencies (and
> all are already installed in this repo's `.venv`). Adopting any of them adds **no new
> top-level dependency** — they ship with `diffrax` regardless.

| Library        | Role here                                                                                 | Verdict |
| -------------- | ----------------------------------------------------------------------------------------- | ------- |
| **Optimistix** | Root-finding for (a) implicit-solver Newton steps (internal to `diffrax`) and (b) equilibration, replacing `scipy.bisect`/`fminbound`. | **Use** |
| **Lineax**     | Linear solves inside the implicit Newton iteration; pulled in transitively by `diffrax`. Only touch directly if we need a custom `linear_solver` for the Jacobian. | Transitive |
| **Equinox**    | `eqx.Module` to model the parcel vector field as a typed pytree (holds `r_drys`, `Nis`, `kappas`, `V`), giving a clean callable RHS and natural `vmap`/`grad` boundaries. Same author as `diffrax`; already a dep. | **Pending — see §5.4** |
| **jaxtyping**  | Shape/dtype annotations (`Float[Array, "nr"]`) — cheap correctness wins; also a `diffrax` dep. | Optional/nice |
| **interpax**   | Smooth interpolation of a time-varying updraft `V(t)` (or tabulated profiles) in a `jit`/`grad`-friendly way, replacing the Python callable. | Use if `V(t)` |

### 5.2 Possible replacements for `diffrax`'s core — and why we still pick `diffrax`

| Alternative                          | Assessment                                                                                                                                  |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------- |
| **`jax.experimental.ode.odeint`**    | Only Dormand–Prince (RK45), **non-stiff explicit**. The parcel ODE is *stiff* (fast condensation relaxation near equilibrium). It would take impractically tiny steps or fail. **Reject as the core solver.** Useful only as an explicit-solver cross-check on non-stiff cases. |
| **Hand-rolled BDF in JAX**           | Reinventing exactly what `diffrax`'s ESDIRK + adaptive control already does, with worse testing. **Reject.**                                |
| **`probdiffeq`** (probabilistic ODE) | Interesting ODE-filter solvers with calibrated uncertainty; genuinely JAX-native. But it's a different numerical paradigm, less battle-tested for stiff microphysics, and harder to validate against CVode. **Not now**; revisit only if uncertainty quantification becomes a goal. |
| **`scipy`/`Assimulo` (status quo)**  | Retired in v2. The `master` branch's `cvode` path is used **offline** to bake golden reference data (§7); it is not a runtime dependency of v2.       |

**Conclusion:** `diffrax` (`Kvaerno5` + `PIDController` + `Optimistix` root finder) is the
right core; `Optimistix` is a confirmed companion; `jax.experimental.ode` is *not* a
viable core replacement for this stiff system.

### 5.3 The differentiability payoff (why this migration is more than a port)

With the RHS in JAX, `diffrax` gives gradients through the whole solve:

- `RecursiveCheckpointAdjoint` (default, discretise-then-optimise) for general
  `d S_max / d θ`.
- `ImplicitAdjoint` (implicit function theorem) if we adopt a steady-state/event-style
  terminal solve.

This directly serves emulator/metamodel training and sensitivity studies — capabilities
the numba/CVode stack cannot offer.

### 5.4 Equinox: concrete walkthrough and recommendation

**What Equinox actually is.** `equinox.Module` is a frozen `dataclass` that is registered
as a JAX *pytree*. That means an instance can be passed straight through `jit`, `grad`,
and `vmap`: its array-valued fields become traced/differentiable leaves, while fields you
mark `static` (ints, strings, Python callables) are treated as compile-time constants.
Equinox also offers `eqx.filter_jit` / `eqx.filter_grad` / `eqx.filter_vmap`, which
auto-partition a pytree into "arrays (trace these)" vs "everything else (hold static)",
so you don't hand-manage `static_argnums`. `diffrax`'s own examples define vector fields
as `eqx.Module`s.

**What it looks like here.**

*Plain-pytree style (no Equinox):*
```python
def parcel_ode_sys(t, y, args):
    r_drys, Nis, kappas, accom, V = args   # positional unpack
    ...
args = (r_drys, Nis, kappas, accom, V)     # a tuple is already a valid pytree
sol = dfx.diffeqsolve(dfx.ODETerm(parcel_ode_sys), ..., args=args)
```

*Equinox style:*
```python
class ParcelVectorField(eqx.Module):
    r_drys: jax.Array
    Nis: jax.Array
    kappas: jax.Array
    accom: float
    V: AbstractUpdraft        # itself an eqx.Module: ConstantV, or an interpolated profile

    def __call__(self, t, y, args):
        ...
        V_t = self.V(t)       # named, self-documenting; V can be const or V(t)
        ...

field = ParcelVectorField(r_drys, Nis, kappas, accom, ConstantV(1.0))
sol = dfx.diffeqsolve(dfx.ODETerm(field), ...)         # field is the pytree of params
dSmax_daccom = eqx.filter_grad(run_to_smax)(field)      # grad w.r.t. array leaves only
```

**Pros (concrete to `pyrcel`):**
- **Named fields kill an existing class of bug.** The current code already tripped over
  positional args — `parcel.py` rewraps the RHS to mutate `args[3]` for `V(t)`, and
  `integrator.py` has a comment about an `args` vs `rhs_args` indexing bug. Named fields
  on a Module remove this footgun.
- **Clean time-varying `V`.** `V` becomes a small `eqx.Module` (`ConstantV`, `StepV`, or
  an `interpax`-backed profile) instead of a closed-over Python callable, keeping the RHS
  pure and `jit`-able — directly solving the §6.7 concern.
- **Differentiation boundaries are explicit.** `eqx.filter_grad`/`partition` make it
  trivial to take gradients w.r.t. *some* parameters (say `accom`, `V`) while holding
  others fixed — exactly the kind of sensitivity study this migration unlocks.
- **Idiomatic + free.** It matches upstream `diffrax` docs/examples (easier maintenance)
  and adds no dependency (§5.1).

**Cons / costs:**
- **Conceptual surface.** Contributors must learn the pytree array/static partition model
  and that Modules are immutable — updates use `eqx.tree_at(...)` rather than attribute
  assignment. For a scientific-Python audience used to mutable objects this is a real (if
  small) ramp.
- **Interacts with the existing mutable API.** Today `set_initial_conditions` mutates a
  `ParcelModel` in place. If the *vector field* is an `eqx.Module`, that's fine (it's
  rebuilt per run anyway), but we should consciously keep the *user-facing* `ParcelModel`
  a normal mutable Python class and use Equinox only for the inner, functional vector-field
  / parameter bundle. Mixing the two paradigms carelessly is the main risk.
- **Marginal benefit if params stay a fixed 5-tuple.** A tuple is not hard to manage; the
  payoff scales with how much we lean on `grad`/`vmap`/`V(t)`.

**Decision (locked):** **Adopt Equinox for the inner vector field and the updraft `V`
abstraction only**, and keep `ParcelModel` itself a plain class. The naming/`V(t)`/grad
ergonomics are a genuine win for a differentiable rewrite, the dependency is already
present, and confining it to the functional core avoids paradigm-mixing in the public API.

---

## 6. Risks, gotchas, and numerical-fidelity concerns

### 6.1 float64 is mandatory (highest-risk gotcha)

JAX defaults to **float32**. Radii are `1e-8–1e-5` m, `S` is `~1e-4`, and `Seq` has a
`r³−r_dry³` cancellation. float32 will *destroy* accuracy and likely make the solver
diverge. We **must** set, at import time, before any array is created:

```python
import jax
jax.config.update("jax_enable_x64", True)
```

The harness must assert `x64` is enabled.

### 6.2 BDF vs ESDIRK ⇒ not bit-identical

CVode uses multistep BDF; `Kvaerno5` is a single-step ESDIRK. Both are A-/L-stable and
appropriate for stiff problems, but they take different internal steps, so trajectories
differ at the tolerance level. The acceptance criterion is *agreement within physically
meaningful tolerances* (e.g. `S_max` within a small relative tolerance, profiles
`allclose`), **not** bitwise equality. This is the central reason for §7.

### 6.3 `max_steps` must be finite under `jit`

`diffeqsolve` under `jit`/autodiff needs a bounded `max_steps`. Stiff activation runs can
take many adaptive steps; set generously (and surface a clear error when exceeded, mapping
to the existing `ParcelModelError`).

### 6.4 Compile latency & performance profile

- First call pays JIT compile cost (seconds); the README already warns users about numba
  warm-up, so this is not a regression in spirit — but recompilation triggers (changing
  `nr`, Python branches) must be controlled.
- For a *single* small run, JAX-on-CPU may not beat numba. The wins are **`vmap` ensembles
  + GPU** and **autodiff**. Set expectations accordingly; benchmark (§7.5).
- **Prefer one adaptive solve over a Python chunk loop (§4.7).** Wrapping a multi-chunk
  Python loop in `jit`/`grad` unrolls it into the graph (longer compile, larger adjoint
  tape); a single `diffeqsolve` keeps the step loop rolled inside XLA with optimal
  checkpointing.
- **Progress meters use host callbacks.** Fine for a single interactive run, but disable
  them (`NoProgressMeter`) under `vmap`/large batches and in the differentiable core to
  avoid per-element host syncs.

### 6.5 `nr` and recompilation

`nr` is encoded as an array shape, not a Python int. Different bin counts ⇒ different
shapes ⇒ recompiles. Acceptable (the existing model also rebuilds state per config), but
avoid making `nr` a `static_argnum` churn point.

### 6.6 Event localization differences

If we use event-based termination (§4.5 option A), the exact `t_smax` is found by a root
find and may land a hair differently than Assimulo's event detection. Validate `S_max`
and `t_smax` explicitly.

### 6.7 Time-varying `V`

Must remain supported. Pass `V` as part of `args` — either a constant or an
`interpax`/`diffrax` interpolation — rather than closing over Python state, so the RHS
stays pure and `jit`-able.

---

## 7. Testing & validation harness

Goal: **prove there are no significant differences** between the `master` model
(numba + CVode) and the v2 JAX + `diffrax` model. Because v2 *deletes* numba/Assimulo, the
`master` model is used **offline** to produce **frozen fixtures**; the v2 test suite then
runs with no dependency on `master`, numba, or Assimulo. Layered from cheapest/strongest
to most integrated.

**Agreed acceptance tolerances:** `S_max` relative tolerance `≤ 1e-3`; activated-fraction
absolute tolerance `≤ 1e-3`. (Tighten opportunistically if the implementation supports it.)

### 7.0 Fixture generation from `master` (oracle, run once per refresh)

A small generator script, run in an environment with `master` checked out (numba +
Assimulo available), writes versioned fixtures into the repo (`tests/fixtures/`):

- **RHS fixtures:** a set of `(y, args)` inputs (sampled physical states + randomized
  states) and the numba `parcel_ode_sys` outputs `dy/dt` for each → `npz`.
- **Equilibration fixtures:** for each config, the `y0` produced by `_setup_run`.
- **Trajectory fixtures:** full CVode trajectories + derived scalars (`S_max`, `t_smax`,
  per-mode activated fraction from `pyrcel/activation.py`) → NetCDF/`npz`.

These fixtures are committed so CI never needs SUNDIALS/numba.

### 7.1 Thermodynamics unit tests (cheap, exact)

Reuse the existing `generate_data.py` + `results.dict` reference pattern (regenerated from
`master`). Assert `thermo_jax.*` matches the frozen numpy/numba values to `~1e-12` over the
existing parameter grids (`dv`, `dv_cont`, `ka`, `es`, `sigma_w`, `Seq`).

### 7.2 RHS pointwise equivalence (strongest single test)

Evaluate the v2 JAX `parcel_ode_sys` at the **frozen `(y, args)` inputs from §7.0** and
assert `allclose` against the frozen numba outputs (tight rtol, e.g. `1e-9`). This isolates
the RHS from the solver and catches transcription bugs immediately — the cornerstone test,
and it runs without numba installed.

### 7.2b Equilibration equivalence

Assert the v2 `optimistix` equilibration reproduces `master`'s frozen `y0` to a tight
tolerance (`~1e-10`). Where we want to isolate *integrator-only* differences (§7.3), seed
the v2 integration with the frozen `y0` so equilibration differences don't leak in.

### 7.3 Golden trajectory regression

1. **Reference data** = the frozen CVode trajectories from §7.0.
2. **Run the v2 model** on the same configs, seeded with the frozen `y0`.
3. **Compare** with the agreed tolerances:
   - `S_max`: relative tolerance `≤ 1e-3`.
   - activated fraction: absolute tolerance `≤ 1e-3`.
   - profiles interpolated to common `output_dt`: `assert_allclose` with per-variable
     tolerances mirroring the solver `atol`.
   Seed configs: the three `examples/*.yml`, the `monte_carlo_outputs.csv` cases, and
   `pm_test.py`.

### 7.4 Config matrix (parametrize with `pytest`)

- bins `nr`: small (mono, 1) and large (250+).
- aerosol modes: single vs multiple species (`simple.yml` has sulfate + mixed).
- updraft `V`: slow/medium/fast (e.g. 0.1, 0.5, 5 m/s) and a `V(t)` case.
- `T0`, `P0`, `S0`: boundary-layer + mid-troposphere points.
- `accom`: 1.0 and 0.1 (as in `pm_test.py`).

### 7.5 Cross-solver & physics invariants

- **Cross-solver bound:** compare `Kvaerno5` vs `Kvaerno3` vs (non-stiff) `Tsit5` on
  non-stiff configs to quantify solver-induced spread — contextualizes the
  CVode↔diffrax differences captured in the frozen fixtures.
- **Invariants:** total water `wv + wc (+ wi)` conservation, `S` single-peaked,
  monotonic `z`. Cheap, solver-agnostic sanity checks that need no oracle.
- **Benchmarks:** wall-clock for a single v2 run (and, for context, the `master` numba
  run captured separately) and `vmap` ensemble (diffrax CPU/GPU) to document the
  performance trade-off honestly.

### 7.6 Differentiability smoke tests

`jax.grad` of `S_max` w.r.t. `y0` / `V` / `accom` runs and is finite; finite-difference
check agrees to a loose tolerance. This exercises the *integration* differentiability that
is the locked requirement (gradients through the solve from the equilibrated `y0`).

### 7.7 CI

The v2 `pixi` environment needs only `jax`, `diffrax` (which pulls `optimistix`,
`equinox`, `lineax`, `jaxtyping`), `numpy`, `pandas`, `xarray`, `pyyaml`. Run §7.1–7.2b +
a small §7.3 subset on every push (fast, no SUNDIALS); run the full matrix nightly/on
demand. Fixture *regeneration* (§7.0) uses a separate, `master`-checked-out environment
with numba + Assimulo and is run manually when the oracle needs refreshing.

---

## 8. Phased rollout

| Phase | Deliverable                                                                                                       | Exit criterion                                  |
| ----- | ----------------------------------------------------------------------------------------------------------------- | ----------------------------------------------- |
| 0     | `v2` feature branch in a private worktree; new `pixi`/deps (`jax`, `diffrax`); enable `x64`; **fixture generator run on `master`** (§7.0). | env builds; `import` works; fixtures committed  |
| 1     | `thermo_jax.py` + `parcel_aux_jax.py` (vectorized RHS); delete numba module.                                      | §7.1 + §7.2 pass (vs frozen fixtures)           |
| 2     | `optimistix` equilibration replacing `scipy`.                                                                     | §7.2b passes (`y0` matches `master`)            |
| 3     | `DiffraxIntegrator` (`Kvaerno5`+`PIDController`), **single adaptive solve** with `SaveAt(ts=...)` (§4.7); rip out Assimulo + `solver_dt` chunking. | §7.3 passes on `examples/*.yml`                 |
| 4     | Event-based termination (§4.5 A) for `S_max` + `terminate_depth`.                                                 | `S_max`/`t_smax` parity; full §7.4 matrix green |
| 5     | `jit`/`vmap`/`grad` ergonomics via Equinox vector field + `V(t)`; differentiability through integration w.r.t. `y0`. | §7.6 + ensemble `vmap` works                    |
| 6     | CLI/interactive mode (`progress_meter` + post-solve summary table); docs, benchmarks, `README`/`pyproject` updates, v2 release prep. | docs updated; perf documented                   |

The branch/worktree, fixture generation, and harness scaffolding (Phase 0) are the agreed
next step once this doc is approved.

---

## 9. Decisions (locked) and the one remaining question

**Locked in:**

1. **v2 / no backwards-compat.** Drop `numba` and Assimulo/`CVode` outright; no
   dual-backend dispatch or legacy paths.
2. **Validate against `master`.** The current numba+CVode model is the oracle, used
   offline to bake frozen fixtures (§7.0); it is not a v2 runtime dependency.
3. **`optimistix` for equilibration**, to keep all numerics in the JAX ecosystem.
   Differentiability is required **through the integration w.r.t. the equilibrated
   `y0`** (and physical params), *not* through the equilibration solve itself.
4. **Acceptance tolerances:** `S_max` rel-tol `≤ 1e-3`, activated-fraction abs-tol
   `≤ 1e-3`.

5. **Equinox — adopted.** Use `eqx.Module` for the inner vector field and the updraft `V`
   abstraction; keep `ParcelModel` a plain class (§5.4). Adds no new dependency.
6. **float64 / device.** Run in float64 (`jax_enable_x64`) for numerical stability;
   **CPU-float64 is the primary correctness target**, with GPU as a `vmap`-ensemble bonus.
7. **Integration structure — single adaptive solve** for the numerical/differentiable
   core, with feedback/output decoupled from integration structure (§4.7). A thin
   interactive-only chunk loop is allowed for the CLI's live stats table, kept out of the
   `jit`/`grad`/`vmap` path.
```
