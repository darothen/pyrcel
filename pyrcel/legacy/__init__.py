"""Reference NumPy/numba implementations of the pyrcel parcel model.

This subpackage preserves the original, non-JAX implementations of the
parcel model and its supporting thermodynamics and activation routines.
It exists for cross-checking and historical reference — not for production
use. The supported API is in the parent ``pyrcel`` package.

Submodules
----------
thermo
    NumPy thermodynamics functions (reference for ``pyrcel.thermo``).
activation
    NumPy/SciPy activation parameterizations (reference for ``pyrcel.activation``).
parcel_aux
    numba-JIT ODE right-hand side (reference for ``pyrcel.parcel_aux``).
parcel
    Legacy ``ParcelModel`` class backed by numba and Assimulo/CVode.
integrator
    Legacy ODE integrator interface (scipy / Assimulo wrappers).
driver
    Legacy ``run_model`` / ``iterate_runs`` convenience wrappers.
"""
