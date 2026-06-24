#!/usr/bin/env python
"""Sensitivity of S_max to aerosol and updraft properties via automatic differentiation.

Scientific motivation
---------------------
The peak supersaturation S_max achieved in an adiabatically rising parcel is the
central quantity linking aerosol properties to cloud droplet number. Classical
activation parameterizations (Abdul-Razzak & Ghan 2000; Morales Betancourt &
Nenes 2014) derive closed-form approximations to S_max as a function of the aerosol
size distribution, hygroscopicity, and updraft speed. With a differentiable parcel
model we can instead compute the *exact* sensitivities

    ∂S_max/∂θ   for θ ∈ {V, N, μ, κ}

by propagating gradients through the full ODE integration. These sensitivities are
directly useful for:

  - validating and tuning activation parameterizations,
  - aerosol retrieval (inverting S_max observations to constrain N or κ), and
  - training neural-network metamodels where gradient information is available
    as training signal.

This example computes the sensitivity table at a single base state, then uses the
gradient to verify a finite-difference perturbation experiment.

Usage
-----
    python examples/differentiable_smax.py
    python examples/differentiable_smax.py --V 0.5 --N 500 --bins 50
"""

from __future__ import annotations

import argparse
import time

import jax
import jax.numpy as jnp
import numpy as np

import pyrcel as pm
from pyrcel.integrator import max_supersaturation

jax.config.update("jax_enable_x64", True)


def differentiable_smax() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--V", type=float, default=1.0, help="base updraft speed (m/s)")
    p.add_argument("--T0", type=float, default=283.15, help="initial temperature (K)")
    p.add_argument("--P0", type=float, default=85000.0, help="initial pressure (Pa)")
    p.add_argument("--S0", type=float, default=-0.02, help="initial supersaturation")
    p.add_argument("--mu", type=float, default=0.05, help="lognormal mode radius (µm)")
    p.add_argument("--sigma", type=float, default=2.0, help="lognormal geometric std dev")
    p.add_argument("--N", type=float, default=1000.0, help="aerosol number (cm⁻³)")
    p.add_argument("--kappa", type=float, default=0.54, help="hygroscopicity κ")
    p.add_argument("--bins", type=int, default=100, help="number of size bins")
    p.add_argument("--t-end", type=float, default=300.0, help="integration time horizon (s)")
    p.add_argument("--n-out", type=int, default=600, help="number of output time points")
    a = p.parse_args()

    # ------------------------------------------------------------------
    # Build the parcel and extract the internal (y0, args) representation
    # ------------------------------------------------------------------
    aerosol = pm.AerosolSpecies(
        "sulfate", pm.Lognorm(mu=a.mu, sigma=a.sigma, N=a.N), kappa=a.kappa, bins=a.bins
    )
    model = pm.ParcelModel([aerosol], V=a.V, T0=a.T0, S0=a.S0, P0=a.P0)
    y0 = jnp.asarray(model.y0)
    ts = jnp.linspace(0.0, a.t_end, a.n_out)

    # Verify the base-state S_max
    r_drys, Nis, kappas, accom, V_obj = model.args
    base_args = (r_drys, Nis, kappas, accom, V_obj)

    print(f"\nBase state: V={a.V} m/s, N={a.N} cm⁻³, μ={a.mu} µm, κ={a.kappa}")

    # ------------------------------------------------------------------
    # Warm up the JIT cache
    # The first call traces and compiles the XLA computation; subsequent
    # calls reuse the cache and are much faster.
    # ------------------------------------------------------------------
    print("\nWarming up JIT cache (first call includes compilation)...", end=" ", flush=True)
    t_start = time.perf_counter()
    _ = jax.block_until_ready(max_supersaturation(y0, base_args, ts))
    print(f"done in {time.perf_counter() - t_start:.1f}s")

    smax_base = float(max_supersaturation(y0, base_args, ts))
    print(f"S_max = {smax_base * 100:.4f}%")

    # ------------------------------------------------------------------
    # Gradient w.r.t. updraft speed V
    # ------------------------------------------------------------------
    # We define a scalar function of V by wrapping max_supersaturation.
    # ConstantV is an equinox Module, so jax.grad treats its field V as
    # a differentiable leaf.
    def smax_of_V(V_val: jax.Array) -> jax.Array:
        args = (r_drys, Nis, kappas, accom, pm.ConstantV(V_val))
        return max_supersaturation(y0, args, ts)

    V_base = jnp.float64(a.V)
    print("\nComputing ∂S_max/∂V via jax.grad...", end=" ", flush=True)
    t0 = time.perf_counter()
    dsmax_dV = float(jax.grad(smax_of_V)(V_base))
    print(f"done in {time.perf_counter() - t0:.2f}s")

    # ------------------------------------------------------------------
    # Gradient w.r.t. total aerosol number N (distributed over bins)
    # ------------------------------------------------------------------
    # Scale all Nis by a scalar factor α; ∂S_max/∂α at α=1 gives the
    # directional derivative along the "uniform scaling" direction, which
    # is proportional to ∂S_max/∂N when N is the total number.
    Nis_base = jnp.asarray(Nis)

    def smax_of_alpha(alpha: jax.Array) -> jax.Array:
        args = (r_drys, alpha * Nis_base, kappas, accom, V_obj)
        return max_supersaturation(y0, args, ts)

    print("Computing ∂S_max/∂N (via uniform scaling) via jax.grad...", end=" ", flush=True)
    t0 = time.perf_counter()
    dsmax_dalpha = float(jax.grad(smax_of_alpha)(jnp.float64(1.0)))
    # ∂S_max/∂N = (∂S_max/∂α) / N_total  [units: %/(cm⁻³)]
    N_total = float(np.sum(np.asarray(Nis))) * 1e-6  # m⁻³ → cm⁻³
    dsmax_dN = dsmax_dalpha / N_total
    print(f"done in {time.perf_counter() - t0:.2f}s")

    # ------------------------------------------------------------------
    # Gradient w.r.t. hygroscopicity κ (uniform shift across all bins)
    # ------------------------------------------------------------------
    kappas_base = jnp.asarray(kappas)

    def smax_of_kappa(kappa_val: jax.Array) -> jax.Array:
        args = (r_drys, Nis_base, jnp.full_like(kappas_base, kappa_val), accom, V_obj)
        return max_supersaturation(y0, args, ts)

    print("Computing ∂S_max/∂κ via jax.grad...", end=" ", flush=True)
    t0 = time.perf_counter()
    dsmax_dkappa = float(jax.grad(smax_of_kappa)(jnp.float64(a.kappa)))
    print(f"done in {time.perf_counter() - t0:.2f}s")

    # ------------------------------------------------------------------
    # Print sensitivity table
    # ------------------------------------------------------------------
    # Convert decimal gradients to % per unit for display
    pct = 100.0
    print("\n" + "=" * 66)
    print("S_max sensitivity table  (∂S_max [%] / ∂θ)")
    print("=" * 66)
    print(f"  {'Parameter':<12}  {'Base value':<14}  {'∂S_max/∂θ  [% / unit]':<24}  Unit of θ")
    print(f"  {'-' * 12}  {'-' * 14}  {'-' * 24}  {'-' * 14}")
    print(f"  {'V':<12}  {a.V:<14.4f}  {dsmax_dV * pct:<24.4e}  m/s")
    print(f"  {'N':<12}  {a.N:<14.1f}  {dsmax_dN * pct:<24.4e}  cm⁻³")
    print(f"  {'κ':<12}  {a.kappa:<14.4f}  {dsmax_dkappa * pct:<24.4e}  —")
    print("=" * 66)

    # ------------------------------------------------------------------
    # Finite-difference verification for V
    # ------------------------------------------------------------------
    dV = a.V * 1e-3
    smax_p = float(smax_of_V(jnp.float64(a.V + dV)))
    smax_m = float(smax_of_V(jnp.float64(a.V - dV)))
    dsmax_dV_fd = (smax_p - smax_m) / (2 * dV)

    rel_err = abs(dsmax_dV - dsmax_dV_fd) / abs(dsmax_dV_fd)
    print("\nFinite-difference check (∂S_max/∂V, decimal units):")
    print(f"  JAX grad:  {dsmax_dV:.6e} /m/s")
    print(f"  FD (±{dV:.4f} m/s): {dsmax_dV_fd:.6e} /m/s")
    print(f"  Relative error: {rel_err:.2e}  {'✓' if rel_err < 0.01 else '✗'}")

    return 0


if __name__ == "__main__":
    differentiable_smax()
