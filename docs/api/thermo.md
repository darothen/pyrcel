# Thermodynamics

JAX implementations of the aerosol/atmospheric thermodynamics functions. All
functions are elementwise, broadcast naturally over arrays, and are differentiable
via `jax.grad`. They are faithful translations of the NumPy reference
implementations in `pyrcel.legacy.thermo`.

!!! note "Legacy reference"
    The original NumPy implementations remain available in `pyrcel.legacy.thermo`
    for cross-checking. They are not part of the v2 numerical path and are not
    differentiable.

::: pyrcel.thermo.sigma_w

---

::: pyrcel.thermo.es

---

::: pyrcel.thermo.Seq

---

::: pyrcel.thermo.rho_air

---

::: pyrcel.thermo.ka

---

::: pyrcel.thermo.dv
