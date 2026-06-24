#!/usr/bin/env python
"""Exact vs numerical ∂Smax/∂V: parcel model, ARG2000, MBN2014.

For each of three methods, computes ∂Smax/∂V two ways on a 10×10 grid of
updraft speed V and lognormal median radius μ:

  * Exact adjoint  — jax.grad propagated through the full computation.
  * Numerical      — central finite differences via numpy.gradient on the
                     pre-computed Smax grid.

For the parameterizations (ARG2000, MBN2014) both paths are cheap; for the
parcel model the adjoint is obtained by calling jax.grad at each grid point
serially, reusing a single JIT-compiled gradient function across all 100
evaluations.

Results are cached to a compressed .npz file; use --recompute to force a
fresh sweep.

Usage
-----
    python examples/sensitivity_sweep.py            # compute, cache, plot
    python examples/sensitivity_sweep.py --recompute
    python examples/sensitivity_sweep.py --no-plot
    python examples/sensitivity_sweep.py --cache PATH --plot PATH
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
_V_MIN, _V_MAX = 0.1, 3.0  # updraft speed (m/s)
_MU_MIN, _MU_MAX = 0.02, 0.30  # median dry radius (µm)

# ── Fixed aerosol / parcel parameters ────────────────────────────────────────
_SIGMA = 2.0
_N = 1000.0  # cm⁻³
_KAPPA = 0.54
_T0 = 283.0  # K
_P0 = 85000.0  # Pa
_S0 = -0.02
_BINS = 50

_DEFAULT_CACHE = "output/sensitivity_sweep_cache.npz"

_REQUIRED_KEYS = {
    "smax_parcel",
    "dsmax_dV_parcel_exact",
    "dsmax_dV_parcel_num",
    "smax_arg",
    "dsmax_dV_arg_exact",
    "dsmax_dV_arg_num",
    "smax_mbn",
    "dsmax_dV_mbn_exact",
    "dsmax_dV_mbn_num",
    "V_grid",
    "mu_grid",
}
_SCALAR_PARAMS = ("sigma", "N", "kappa", "T0", "P0")


def sensitivity_sweep() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-V", type=int, default=10)
    p.add_argument("--n-mu", type=int, default=10)
    p.add_argument("--sigma", type=float, default=_SIGMA)
    p.add_argument("--N", type=float, default=_N, help="aerosol number (cm⁻³)")
    p.add_argument("--kappa", type=float, default=_KAPPA)
    p.add_argument("--T0", type=float, default=_T0)
    p.add_argument("--P0", type=float, default=_P0)
    p.add_argument("--cache", type=str, default=_DEFAULT_CACHE, metavar="PATH")
    p.add_argument("--recompute", action="store_true")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("--plot", type=str, default=None, metavar="PATH")
    a = p.parse_args()

    V_grid = np.logspace(np.log10(_V_MIN), np.log10(_V_MAX), a.n_V)
    mu_grid = np.logspace(np.log10(_MU_MIN), np.log10(_MU_MAX), a.n_mu)

    cache_path = Path(a.cache)
    data: dict | None = None

    if not a.recompute and cache_path.exists():
        cached = np.load(cache_path, allow_pickle=False)
        if _cache_valid(cached, V_grid, mu_grid, a):
            print(f"Loaded from cache: {cache_path}")
            data = dict(cached)
        else:
            print("Cache parameters differ — recomputing.")

    if data is None:
        data = _compute(V_grid, mu_grid, a)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache_path, **data)
        print(f"\nCache saved to {cache_path}")

    if not a.no_plot:
        fig_path = a.plot or str(cache_path.with_suffix(".png"))
        _plot(data, fig_path)

    return 0


def _cache_valid(cached, V_grid, mu_grid, a) -> bool:
    try:
        return (
            _REQUIRED_KEYS.issubset(set(cached.files))
            and np.allclose(cached["V_grid"], V_grid)
            and np.allclose(cached["mu_grid"], mu_grid)
            and all(np.isclose(float(cached[k]), getattr(a, k)) for k in _SCALAR_PARAMS)
        )
    except (KeyError, ValueError):
        return False


def _compute(V_grid: np.ndarray, mu_grid: np.ndarray, a) -> dict:
    n_V, n_mu = len(V_grid), len(mu_grid)
    V_jax = jnp.asarray(V_grid, dtype=jnp.float64)
    mu_jax = jnp.asarray(mu_grid, dtype=jnp.float64)

    print(f"\nGrid: {n_V} × {n_mu}  (V ∈ [{_V_MIN}, {_V_MAX}] m/s,  μ ∈ [{_MU_MIN}, {_MU_MAX}] µm)")
    print(f"Aerosol: N = {a.N:.0f} cm⁻³, σ = {a.sigma}, κ = {a.kappa}\n")

    # ── ARG2000 and MBN2014 — vectorised over full grid ───────────────────────
    def _smax_arg(V, mu):
        s, _, _ = arg2000(
            V,
            a.T0,
            a.P0,
            jnp.array([mu]),
            jnp.array([a.sigma]),
            jnp.array([a.N]),
            jnp.array([a.kappa]),
        )
        return s

    def _smax_mbn(V, mu):
        s, _, _ = mbn2014(
            V,
            a.T0,
            a.P0,
            jnp.array([mu]),
            jnp.array([a.sigma]),
            jnp.array([a.N]),
            jnp.array([a.kappa]),
        )
        return s

    def _grid2(fn):
        return jax.jit(jax.vmap(jax.vmap(fn, in_axes=(None, 0)), in_axes=(0, None)))

    def _dV_grid2(fn):
        return _grid2(lambda V, mu: jax.grad(fn, argnums=0)(V, mu))

    print("ARG2000 … ", end="", flush=True)
    t0 = time.perf_counter()
    smax_arg = np.array(_grid2(_smax_arg)(V_jax, mu_jax))
    dsmax_dV_arg_ex = np.array(_dV_grid2(_smax_arg)(V_jax, mu_jax))
    print(f"{time.perf_counter() - t0:.1f}s")

    print("MBN2014 … ", end="", flush=True)
    t0 = time.perf_counter()
    smax_mbn = np.array(_grid2(_smax_mbn)(V_jax, mu_jax))
    dsmax_dV_mbn_ex = np.array(_dV_grid2(_smax_mbn)(V_jax, mu_jax))
    print(f"{time.perf_counter() - t0:.1f}s")

    dsmax_dV_arg_num = np.gradient(smax_arg, V_grid, axis=0)
    dsmax_dV_mbn_num = np.gradient(smax_mbn, V_grid, axis=0)

    # ── Parcel model ──────────────────────────────────────────────────────────
    # Use generous t_end; the event-based solver stops at S_max automatically.
    t_end = 5.0 * 300.0 / _V_MIN
    ts_global = jnp.linspace(0.0, t_end, max(600, int(t_end)))

    # accom is a physical constant, same for every model build.
    from pyrcel import constants as _c

    _accom = _c.ac

    # Explicit-arg JIT: shapes are fixed across all (V, μ) calls so the kernel
    # is compiled once on the first call and reused for all 100 evaluations.
    def _smax_parcel(V_val, y0, r_drys, Nis, kappas_arr):
        return max_supersaturation(
            y0, (r_drys, Nis, kappas_arr, _accom, ConstantV(V_val)), ts_global
        )

    _fwd = jax.jit(_smax_parcel)
    _grad = jax.jit(jax.grad(_smax_parcel, argnums=0))

    smax_parcel = np.zeros((n_V, n_mu))
    dsmax_dV_parcel_ex = np.full((n_V, n_mu), np.nan)

    total = n_V * n_mu
    done = 0
    print(f"Parcel model: {total} forward + adjoint evaluations")

    for i_mu, mu_val in enumerate(mu_grid):
        aerosol = pm.AerosolSpecies(
            "sulfate",
            pm.Lognorm(mu=mu_val, sigma=a.sigma, N=a.N),
            kappa=a.kappa,
            bins=_BINS,
        )
        mdl = pm.ParcelModel([aerosol], V=1.0, T0=a.T0, S0=_S0, P0=a.P0)
        y0_j = jnp.asarray(mdl.y0, dtype=jnp.float64)
        rd_j = jnp.asarray(mdl.args[0], dtype=jnp.float64)
        ni_j = jnp.asarray(mdl.args[1], dtype=jnp.float64)
        ka_j = jnp.asarray(mdl.args[2], dtype=jnp.float64)

        for i_V, V_val in enumerate(V_grid):
            V_j = jnp.float64(V_val)
            cold = done == 0

            label = f"  [{done + 1:03d}/{total}] V={V_val:.2f} m/s  μ={mu_val:.4f} µm"
            if cold:
                label += "  (cold JIT — compiling forward + adjoint …)"
            print(label, end="", flush=True)
            t0 = time.perf_counter()

            try:
                s = float(jax.block_until_ready(_fwd(V_j, y0_j, rd_j, ni_j, ka_j)))
                g = float(jax.block_until_ready(_grad(V_j, y0_j, rd_j, ni_j, ka_j)))
            except Exception as exc:
                print(f"  ✗ ({exc.__class__.__name__})")
                s, g = np.nan, np.nan
            else:
                print(f"  {time.perf_counter() - t0:.1f}s")

            smax_parcel[i_V, i_mu] = s
            dsmax_dV_parcel_ex[i_V, i_mu] = g
            done += 1

    dsmax_dV_parcel_num = np.gradient(smax_parcel, V_grid, axis=0)

    return {
        "V_grid": V_grid,
        "mu_grid": mu_grid,
        "sigma": np.float64(a.sigma),
        "N": np.float64(a.N),
        "kappa": np.float64(a.kappa),
        "T0": np.float64(a.T0),
        "P0": np.float64(a.P0),
        "smax_parcel": smax_parcel,
        "dsmax_dV_parcel_exact": dsmax_dV_parcel_ex,
        "dsmax_dV_parcel_num": dsmax_dV_parcel_num,
        "smax_arg": smax_arg,
        "dsmax_dV_arg_exact": dsmax_dV_arg_ex,
        "dsmax_dV_arg_num": dsmax_dV_arg_num,
        "smax_mbn": smax_mbn,
        "dsmax_dV_mbn_exact": dsmax_dV_mbn_ex,
        "dsmax_dV_mbn_num": dsmax_dV_mbn_num,
    }


def _plot(data: dict, path: str) -> None:
    try:
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot")
        return

    V_grid = data["V_grid"]
    mu_grid = data["mu_grid"]
    pct = 100.0  # fraction → %

    row_defs = [
        {
            "label": "Parcel model",
            "exact": data["dsmax_dV_parcel_exact"] * pct,
            "num": data["dsmax_dV_parcel_num"] * pct,
        },
        {
            "label": "ARG2000",
            "exact": data["dsmax_dV_arg_exact"] * pct,
            "num": data["dsmax_dV_arg_num"] * pct,
        },
        {
            "label": "MBN2014",
            "exact": data["dsmax_dV_mbn_exact"] * pct,
            "num": data["dsmax_dV_mbn_num"] * pct,
        },
    ]
    col_titles = [
        "Exact adjoint  (jax.grad)",
        r"Numerical  (np.gradient on grid)",
    ]

    fig, axes = plt.subplots(3, 2, figsize=(9, 9), squeeze=False)
    fig.subplots_adjust(hspace=0.30, wspace=0.08)

    for i_row, rdef in enumerate(row_defs):
        exact = rdef["exact"]
        num = rdef["num"]

        # Shared colour limits: symmetric around 0, using finite values only
        all_vals = np.concatenate([exact[np.isfinite(exact)], num[np.isfinite(num)]])
        lim = np.nanpercentile(np.abs(all_vals), 98)
        norm = mcolors.Normalize(vmin=0, vmax=lim)

        for i_col, arr in enumerate([exact, num]):
            ax = axes[i_row, i_col]
            im = ax.pcolormesh(
                V_grid,
                mu_grid,
                arr.T,
                cmap="plasma",
                norm=norm,
                shading="auto",
            )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.tick_params(labelsize=8)

            if i_row == 0:
                ax.set_title(col_titles[i_col], fontsize=10, pad=7)
            if i_row < 2:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel("$V$ (m s$^{-1}$)", labelpad=4)
            if i_col == 0:
                ax.set_ylabel(rdef["label"] + "\n$\\mu$ (µm)", labelpad=4)
            else:
                ax.set_yticklabels([])

        cbar = fig.colorbar(
            im,
            ax=axes[i_row, :].tolist(),
            location="right",
            shrink=0.85,
            pad=0.02,
        )
        cbar.set_label(
            r"$\partial S_{\mathrm{max}}/\partial V$  (% (m s$^{-1}$)$^{-1}$)",
            fontsize=9,
            labelpad=8,
        )
        cbar.ax.tick_params(labelsize=8)

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved figure to {path}")
    plt.close(fig)


if __name__ == "__main__":
    sensitivity_sweep()
