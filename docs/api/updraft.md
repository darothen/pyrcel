# Updraft

Updraft speed specification. `ParcelModel` accepts a scalar (constant speed), a
`ConstantV` instance, or an `InterpolatedUpdraft` for time-varying profiles.
All types are JAX pytrees and are compatible with `jax.grad` and `jax.vmap`.

::: pyrcel.updraft.AbstractUpdraft

---

::: pyrcel.updraft.ConstantV

---

::: pyrcel.updraft.InterpolatedUpdraft

---

::: pyrcel.updraft.as_updraft
