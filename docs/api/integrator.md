# Integrator

Low-level ODE integration functions. These are the building blocks used by
`ParcelModel.run()` and are public for advanced workflows: custom optimization
loops, partial pipelines, or composing with external JAX transforms.

All functions are JAX-traceable and compatible with `jax.grad` / `jax.vmap`
unless noted otherwise.

!!! tip "When to use these directly"
    Use the integrator functions directly when you need:

    - Gradients through the ODE (use `max_supersaturation` or `nd_from_parcel`)
    - A batched `vmap` kernel without the `ParcelModel` overhead
    - The raw `diffrax.Solution` object
    - Fine-grained control over `ts`, tolerances, or step limits

## Core solve

::: pyrcel.integrator.integrate_parcel

---

::: pyrcel.integrator.integrate_parcel_arrays

## Differentiable diagnostics

::: pyrcel.integrator.max_supersaturation

---

::: pyrcel.integrator.nd_from_parcel

## Terminated-run pipeline

::: pyrcel.integrator.find_smax

!!! warning "Not differentiable"
    `find_smax` uses event detection (`dS/dt = 0`) which is discontinuous and
    therefore not compatible with `jax.grad`. Use `max_supersaturation` on the
    differentiable path.

---

::: pyrcel.integrator.terminate_cutoff_time

---

::: pyrcel.integrator.integrate_parcel_terminated

## Tolerance constants

::: pyrcel.integrator.atol_vector

The module also exports `STATE_RTOL`, `STATE_ATOL`, and `RADIUS_ATOL` as
importable constants that match the CVode configuration from the legacy model.
