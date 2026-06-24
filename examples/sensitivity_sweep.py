#!/usr/bin/env python
"""Smax sensitivity sweep: parcel model vs ARG2000 vs MBN2014.

Computes Smax and its sensitivities ∂Smax/∂V and ∂Smax/∂μ on a dense 2D
grid of updraft speed V and lognormal median radius μ. Three methods are
compared:

  * Parcel model — exact ODE integration (JAX/diffrax), differentiable via
    automatic differentiation.
  * ARG2000 — Abdul-Razzak & Ghan (2000) closed-form parameterization.
  * MBN2014 — Morales Betancourt & Nenes (2014) parameterization.

∂Smax/∂V is computed via jax.grad for all three methods.
∂Smax/∂μ is computed via jax.grad for ARG2000 and MBN2014 (analytical);
for the parcel model it is approximated by numerical differentiation on
the grid using numpy.gradient.

Results are cached in a compressed .npz file so the expensive parcel
integrations are not repeated on subsequent calls. Use --recompute to
force a fresh sweep.

Usage
-----
    python examples/sensitivity_sweep.py            # compute, cache, plot
    python examples/sensitivity_sweep.py --recompute  # force fresh run
    python examples/sensitivity_sweep.py --no-plot    # skip figure
    python examples/sensitivity_sweep.py --cache PATH # custom cache
    python examples/sensitivity_sweep.py --plot PATH  # custom figure path
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

import pyrcel as pm
from pyrcel.activation import arg2000, mbn2014
from pyrcel.integrator import max_supersaturation
from pyrcel.updraft import ConstantV

jax.config.update("jax_enable_x64", True)

# ── Grid extents ──────────────────────────────────────────────────────────────
_V_MIN,  _V_MAX  = 0.1, 3.0    # updraft speed (m/s)
_MU_MIN, _MU_MAX = 0.02, 0.30  # median dry radius (µm)

# ── Fixed aerosol / parcel parameters ────────────────────────────────────────
_SIGMA = 2.0
_N     = 1000.0   # cm⁻³
_KAPPA = 0.54
_T0    = 283.0    # K
_P0    = 85000.0  # Pa
_S0    = -0.02
_BINS  = 50

_DEFAULT_CACHE = "output/sensitivity_sweep_cache.npz"


def sensitivity_sweep() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-V",  type=int,   default=25,                help="grid points along V")
    p.add_argument("--n-mu", type=int,   default=25,                help="grid points along μ")
    p.add_argument("--sigma",  type=float, default=_SIGMA)
    p.add_argument("--N",      type=float, default=_N,              help="aerosol number (cm⁻³)")
    p.add_argument("--kappa",  type=float, default=_KAPPA)
    p.add_argument("--T0",     type=float, default=_T0)
    p.add_argument("--P0",     type=float, default=_P0)
    p.add_argument("--cache",  type=str,   default=_DEFAULT_CACHE,  metavar="PATH")
    p.add_argument("--recompute", action="store_true",              help="ignore existing cache")
    p.add_argument("--no-plot",   action="store_true",              help="skip figure generation")
    p.add_argument("--plot",   type=str,   default=None,            metavar="PATH",
                   help="figure output path (default: cache stem + .png)")
    a = p.parse_args()

    V_grid  = np.logspace(np.log10(_V_MIN),  np.log10(_V_MAX),  a.n_V)
    mu_grid = np.logspace(np.log10(_MU_MIN), np.log10(_MU_MAX), a.n_mu)

    cache_path = Path(a.cache)
    data: dict | None = None

    if not a.recompute and cache_path.exists():
        cached = np.load(cache_path, allow_pickle=False)
        if _cache_valid(cached, V_grid, mu_grid, a):
            print(f"Loaded from cache: {cache_path}")
            data = dict(cached)
        else:
            print("Cache exists but parameters differ — recomputing.")

    if data is None:
        data = _compute(V_grid, mu_grid, a)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache_path, **data)
        print(f"\nCache saved to {cache_path}")

    if not a.no_plot:
        fig_path = a.plot or str(cache_path.with_suffix(".png"))
        _plot(data, fig_path)

    return 0


# ── Cache validation ──────────────────────────────────────────────────────────

_REQUIRED_KEYS = {
    "smax_parcel", "dsmax_dV_parcel", "dsmax_dmu_parcel",
    "smax_arg",    "dsmax_dV_arg",    "dsmax_dmu_arg",
    "smax_mbn",    "dsmax_dV_mbn",    "dsmax_dmu_mbn",
    "V_grid", "mu_grid",
}
_SCALAR_PARAMS = ("sigma", "N", "kappa", "T0", "P0")


def _cache_valid(cached: np.lib.npyio.NpzFile, V_grid, mu_grid, a) -> bool:
    try:
        return (
            _REQUIRED_KEYS.issubset(set(cached.files))
            and np.allclose(cached["V_grid"], V_grid)
            and np.allclose(cached["mu_grid"], mu_grid)
            and all(np.isclose(float(cached[k]), getattr(a, k)) for k in _SCALAR_PARAMS)
        )
    except (KeyError, ValueError):
        return False


# ── Computation ───────────────────────────────────────────────────────────────

def _compute(V_grid: np.ndarray, mu_grid: np.ndarray, a) -> dict:
    n_V, n_mu = len(V_grid), len(mu_grid)
    V_jax  = jnp.asarray(V_grid,  dtype=jnp.float64)
    mu_jax = jnp.asarray(mu_grid, dtype=jnp.float64)

    print(f"\nSweep:  {n_V} × {n_mu} grid  "
          f"(V ∈ [{_V_MIN}, {_V_MAX}] m/s,  μ ∈ [{_MU_MIN}, {_MU_MAX}] µm)")
    print(f"Aerosol: N = {a.N:.0f} cm⁻³, σ = {a.sigma}, κ = {a.kappa}\n")

    # ── Analytical parameterizations — vectorised over full grid ──────────────
    # Helper: double-vmap fn(V_scalar, mu_scalar) → scalar over (n_V, n_mu)
    def _grid2(fn):
        return jax.jit(
            jax.vmap(jax.vmap(fn, in_axes=(None, 0)), in_axes=(0, None))
        )

    def _grads2(fn):
        """Grid of (∂fn/∂V, ∂fn/∂μ)."""
        return _grid2(lambda V, mu: jax.grad(fn, argnums=(0, 1))(V, mu))

    def _smax_arg(V, mu):
        s, _, _ = arg2000(
            V, a.T0, a.P0,
            jnp.array([mu]), jnp.array([a.sigma]),
            jnp.array([a.N]), jnp.array([a.kappa]),
        )
        return s

    def _smax_mbn(V, mu):
        s, _, _ = mbn2014(
            V, a.T0, a.P0,
            jnp.array([mu]), jnp.array([a.sigma]),
            jnp.array([a.N]), jnp.array([a.kappa]),
        )
        return s

    print("ARG2000: forward pass … ", end="", flush=True)
    t0 = time.perf_counter()
    smax_arg = np.array(_grid2(_smax_arg)(V_jax, mu_jax))
    print(f"{time.perf_counter()-t0:.1f}s", end="   |   ", flush=True)
    print("gradients … ", end="", flush=True)
    t0 = time.perf_counter()
    dV_a, dmu_a = _grads2(_smax_arg)(V_jax, mu_jax)
    dsmax_dV_arg, dsmax_dmu_arg = np.array(dV_a), np.array(dmu_a)
    print(f"{time.perf_counter()-t0:.1f}s")

    print("MBN2014: forward pass … ", end="", flush=True)
    t0 = time.perf_counter()
    smax_mbn = np.array(_grid2(_smax_mbn)(V_jax, mu_jax))
    print(f"{time.perf_counter()-t0:.1f}s", end="   |   ", flush=True)
    print("gradients … ", end="", flush=True)
    t0 = time.perf_counter()
    dV_m, dmu_m = _grads2(_smax_mbn)(V_jax, mu_jax)
    dsmax_dV_mbn, dsmax_dmu_mbn = np.array(dV_m), np.array(dmu_m)
    print(f"{time.perf_counter()-t0:.1f}s")

    # ── Parcel model — vmap over V for each μ value ───────────────────────────
    # t_end generous enough for the slowest updraft to reach S_max
    t_end = 5.0 * 300.0 / _V_MIN
    ts_global = jnp.linspace(0.0, t_end, max(600, int(t_end)))

    # Get accommodation coefficient from a representative model build
    _aero0 = pm.AerosolSpecies(
        "sulfate", pm.Lognorm(mu=mu_grid[0], sigma=a.sigma, N=a.N),
        kappa=a.kappa, bins=_BINS,
    )
    _, _, _, accom_const, _ = pm.ParcelModel([_aero0], V=1.0, T0=a.T0, S0=_S0, P0=a.P0).args

    # JIT-compiled once; shapes are fixed (V: (n_V,), y0: (7+bins,),
    # r_drys/Nis/kappas: (bins,)) so cache is reused across all μ iterations.
    @jax.jit
    def _batch(V_arr, y0, r_drys, Nis, kappas_arr):
        def smax_fn(V_val):
            return max_supersaturation(
                y0,
                (r_drys, Nis, kappas_arr, accom_const, ConstantV(V_val)),
                ts_global,
            )
        return (
            jax.vmap(smax_fn)(V_arr),
            jax.vmap(jax.grad(smax_fn))(V_arr),
        )

    print(f"\nParcel model: iterating over {n_mu} μ values (vmapped over {n_V} V values each)")
    smax_parcel    = np.zeros((n_V, n_mu))
    dsmax_dV_parcel = np.zeros((n_V, n_mu))

    for i_mu, mu_val in enumerate(mu_grid):
        aerosol = pm.AerosolSpecies(
            "sulfate",
            pm.Lognorm(mu=mu_val, sigma=a.sigma, N=a.N),
            kappa=a.kappa,
            bins=_BINS,
        )
        mdl = pm.ParcelModel([aerosol], V=1.0, T0=a.T0, S0=_S0, P0=a.P0)
        y0_jax       = jnp.asarray(mdl.y0,      dtype=jnp.float64)
        r_drys_jax   = jnp.asarray(mdl.args[0], dtype=jnp.float64)
        Nis_jax      = jnp.asarray(mdl.args[1], dtype=jnp.float64)
        kappas_jax   = jnp.asarray(mdl.args[2], dtype=jnp.float64)

        tag = "cold" if i_mu == 0 else "warm"
        print(f"  [{i_mu+1:02d}/{n_mu}] μ = {mu_val:.4f} µm ({tag}) … ", end="", flush=True)
        t0 = time.perf_counter()

        smax_row, dsmax_dV_row = jax.block_until_ready(
            _batch(V_jax, y0_jax, r_drys_jax, Nis_jax, kappas_jax)
        )
        print(f"{time.perf_counter()-t0:.1f}s")

        smax_parcel[:, i_mu]    = np.asarray(smax_row)
        dsmax_dV_parcel[:, i_mu] = np.asarray(dsmax_dV_row)

    # ∂Smax/∂μ for parcel: numerical differentiation across the μ axis.
    # numpy.gradient handles non-uniform (log-spaced) grids correctly.
    dsmax_dmu_parcel = np.gradient(smax_parcel, mu_grid, axis=1)

    return {
        "V_grid":  V_grid,
        "mu_grid": mu_grid,
        "sigma":   np.float64(a.sigma),
        "N":       np.float64(a.N),
        "kappa":   np.float64(a.kappa),
        "T0":      np.float64(a.T0),
        "P0":      np.float64(a.P0),
        "smax_parcel":      smax_parcel,
        "dsmax_dV_parcel":  dsmax_dV_parcel,
        "dsmax_dmu_parcel": dsmax_dmu_parcel,
        "smax_arg":         smax_arg,
        "dsmax_dV_arg":     dsmax_dV_arg,
        "dsmax_dmu_arg":    dsmax_dmu_arg,
        "smax_mbn":         smax_mbn,
        "dsmax_dV_mbn":     dsmax_dV_mbn,
        "dsmax_dmu_mbn":    dsmax_dmu_mbn,
    }


# ── Figure ────────────────────────────────────────────────────────────────────

def _plot(data: dict, path: str) -> None:
    try:
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
    except ImportError:
        print("matplotlib not available; skipping plot")
        return

    V_grid  = data["V_grid"]
    mu_grid = data["mu_grid"]
    pct = 100.0  # fraction → %

    row_defs = [
        {
            "label": "$S_{\\mathrm{max}}$ (%)",
            "keys":  ("smax_parcel",     "smax_arg",     "smax_mbn"),
            "scale": pct,
            "cmap":  "viridis",
            "sym":   False,
        },
        {
            "label": r"$\partial S_{\mathrm{max}}/\partial V$" "\n"
                     r"(% (m s$^{-1}$)$^{-1}$)",
            "keys":  ("dsmax_dV_parcel",  "dsmax_dV_arg",  "dsmax_dV_mbn"),
            "scale": pct,
            "cmap":  "plasma",
            "sym":   False,
        },
        {
            "label": r"$\partial S_{\mathrm{max}}/\partial \mu$" "\n"
                     r"(% µm$^{-1}$)",
            "keys":  ("dsmax_dmu_parcel", "dsmax_dmu_arg", "dsmax_dmu_mbn"),
            "scale": pct,
            "cmap":  "RdBu_r",
            "sym":   True,
        },
    ]
    col_titles = ["Parcel model", "ARG2000", "MBN2014"]

    fig, axes = plt.subplots(3, 3, figsize=(11, 9), squeeze=False)
    fig.subplots_adjust(hspace=0.35, wspace=0.07)

    for i_row, rdef in enumerate(row_defs):
        arrays = [data[k] * rdef["scale"] for k in rdef["keys"]]
        all_vals = np.concatenate([arr.ravel() for arr in arrays])
        vmin, vmax = np.nanpercentile(all_vals, [2, 98])
        if rdef["sym"]:
            lim = max(abs(vmin), abs(vmax))
            vmin, vmax = -lim, lim
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        for i_col, (arr, ctitle) in enumerate(zip(arrays, col_titles)):
            ax = axes[i_row, i_col]
            im = ax.pcolormesh(
                V_grid, mu_grid, arr.T,
                cmap=rdef["cmap"], norm=norm, shading="auto",
            )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.tick_params(labelsize=8)

            if i_row == 0:
                ax.set_title(ctitle, fontsize=11, fontweight="bold", pad=7)
            if i_row < 2:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel("$V$ (m s$^{-1}$)", labelpad=4)
            if i_col == 0:
                ax.set_ylabel("$\\mu$ (µm)", labelpad=4)
            else:
                ax.set_yticklabels([])

        cbar = fig.colorbar(
            im,
            ax=axes[i_row, :].tolist(),
            location="right",
            shrink=0.85,
            pad=0.02,
        )
        cbar.set_label(rdef["label"], fontsize=9, labelpad=8)
        cbar.ax.tick_params(labelsize=8)

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved figure to {path}")
    plt.close(fig)


if __name__ == "__main__":
    sensitivity_sweep()
