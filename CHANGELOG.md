# Changelog

All notable changes to pyrcel are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Version numbers follow [Semantic Versioning](https://semver.org/).

---

## [2.0.0] — unreleased

This is a major release. The numerical core has been rewritten from scratch
on [JAX](https://github.com/google/jax) +
[diffrax](https://github.com/patrick-kidger/diffrax), replacing the
numba + Assimulo/SUNDIALS stack from v1. The public interface is largely
unchanged; see the [migration guide](https://pyrcel.readthedocs.io/en/latest/user_guide/migration/)
for the handful of breaking changes.

### Added

- **Differentiable integration** — the full parcel model is now
  differentiable end-to-end via `jax.grad`. Gradients flow through the
  ODE solve, the equilibration step, and the activation diagnostics.
- **Batched ensemble runs** — `run_updraft_ensemble` and
  `smax_nact_ensemble` run an array of independent parcel solves in a
  single compiled `jax.vmap` kernel with no Python loop overhead.
- **GPU acceleration** — pass `device="gpu"` to `ParcelModel` or set
  `JAX_PLATFORM_NAME=gpu`; install with `pip install "pyrcel[gpu]"`.
- **`ModelOutput`** — `ParcelModel.run()` returns a structured
  `ModelOutput` object with `.to_pandas()`, `.to_polars()`,
  `.to_xarray()`, `.to_netcdf()`, `.to_csv()`, and `.to_parquet()`.
- **Short-circuit run modes** — `run(mode="smax")` returns the peak
  supersaturation as a plain float (ideal for the differentiable path);
  `run(mode="nd")` returns the activated droplet number concentration.
- **`InterpolatedUpdraft`** — a JAX-traceable, `vmap`-safe time-varying
  updraft specified by `(ts, Vs)` arrays. Replaces bare callables from v1.
- **JAX-native activation parameterizations** — `ARG2000` and `MBN2014`
  are now differentiable callable classes in `pyrcel.activation`; the
  underlying implementations are fully `jax.jit`- and `jax.grad`-compatible.
- **Hermite cubic S_max refinement** — S_max is now located via a
  single-kernel Hermite cubic interpolant of the terminal supersaturation
  trajectory, giving a grid-independent peak estimate without a second ODE
  solve.
- **`console=True` mode** — pretty-prints setup tables, phase timers
  (with first-compile detection), trajectory sample, and activation summary
  without adding any overhead to the compiled solve.
- **`live=True` mode** — streams a z/T/S step table during integration,
  replicating the interactive output of the v1 CVode integrator.
- **`pyrcel.legacy` subpackage** — the original NumPy/SciPy
  implementations of thermodynamics and activation functions are preserved
  here for cross-validation; they are not part of the v2 numerical path.
- **MkDocs/Material documentation site** — full rebuild with API
  reference (mkdocstrings), migration guide, numerical methods description,
  worked examples, and GPU setup guide.
- **Python 3.11, 3.12, and 3.13 support.**
- **PEP 561 type annotations** throughout the public API, validated with
  pyrefly.

### Changed

- **ODE solver**: SUNDIALS CVode (BDF, via Assimulo) → `diffrax.Kvaerno5`
  (ESDIRK-5/4). Results agree within solver tolerances (rtol = 1e-7); see
  [numerical differences](https://pyrcel.readthedocs.io/en/latest/user_guide/migration/#known-numerical-differences).
- **`ParcelModel.run()` return value**: was `(parcel_df, aer_dfs)` tuple;
  is now a `ModelOutput` object. Call `.to_pandas()` to recover the v1 tuple.
- **Activation parameterizations**: `pyrcel.activation.arg2000` and
  `pyrcel.activation.mbn2014` are now thin wrappers; prefer the new
  `ARG2000()` and `MBN2014()` callable classes for new code.
- **`ParcelModel` constructor**: two new optional kwargs (`console`,
  `device`); all existing call-sites are unaffected.
- **Equilibration**: `scipy.optimize.bisect` → `optimistix.Bisection` (same
  algorithm, tighter tolerance, gradients supported).
- **Build / install tooling**: pixi and CircleCI replaced by uv and GitHub
  Actions.
- **Kinetic activation gate**: `binned_activation` now restricts the
  "grown past critical radius" search to equilibrium-activated bins only,
  preventing false positives from the approximate Köhler formula for tiny
  particles near the r_crit ≈ r_dry transition.

### Removed

- **Assimulo / SUNDIALS dependency** — gone entirely. No Fortran, no native
  extension compilation, no conda-forge workarounds. `pip install pyrcel`
  is sufficient.
- **numba dependency** — the JIT layer is now `jax.jit` / XLA.
- **`solver=` parameter on `ParcelModel.run()`** — there is one solver;
  the argument is silently ignored if passed (will become an error in a
  future release).
- **`ParcelModelJAX` class name** — still importable as an alias for
  `ParcelModel`; will be removed in v2.1.
- **Python < 3.11 support.**

### Fixed

- Mean/median mismatch in aerosol mode summary table (`#26`).
- Post-hoc Hermite cubic S_max finder returned NaN for runs without
  termination or with `live=True` (`#80`).
- `configure_logging()` set `log.propagate = False`, silently breaking
  pytest `caplog` capture of pyrcel log records.

---

## [1.3.5] — 2022-10-17

Last release on the v1 (numba + CVode) branch. See the
[v1.3.5 tag](https://github.com/darothen/pyrcel/releases/tag/v1.3.5) for
full details. To install: `pip install "pyrcel==1.3.5"`.

---

[2.0.0]: https://github.com/darothen/pyrcel/compare/v1.3.5...v2.0.0
[1.3.5]: https://github.com/darothen/pyrcel/releases/tag/v1.3.5
