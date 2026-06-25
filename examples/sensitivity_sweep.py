#!/usr/bin/env python
"""Exact vs numerical ∂Smax/∂V: parcel model, ARG2000, MBN2014.

For each of three methods, computes ∂Smax/∂V two ways on a 10×10 grid of
updraft speed V and lognormal median radius μ:

  * Exact adjoint  — jax.grad propagated through the full computation.
  * Numerical      — central differences on an extended (V, μ) grid (one extra
                     log-spaced node in each direction) so every original node
                     uses a central difference; values are interpolated back via
                     bilinear interpolation in log(V) × log(μ) space.

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


def _extend_grid_log(grid: np.ndarray) -> np.ndarray:
    """Extend a log-spaced 1-D grid by one node at each end, preserving spacing."""
    log_g = np.log10(grid)
    d = log_g[1] - log_g[0]
    return 10 ** np.concatenate([[log_g[0] - d], log_g, [log_g[-1] + d]])


def _num_grad_log_interp(
    smax_ext: np.ndarray,
    V_ext: np.ndarray,
    mu_ext: np.ndarray,
    V_grid: np.ndarray,
    mu_grid: np.ndarray,
) -> np.ndarray:
    """d(Smax)/d(V) via central FD on extended grid, log-bilinear interp to original nodes.

    Central differences on the (n_V+2)×(n_μ+2) extended grid give second-order
    accuracy at every original node (none sit on a boundary any more).
    Interpolation is bilinear in log(V)×log(μ) space, matching the log-spaced grids.
    """
    from scipy.interpolate import RegularGridInterpolator

    dsmax_dV_ext = np.gradient(smax_ext, V_ext, axis=0)
    interp = RegularGridInterpolator(
        (np.log10(V_ext), np.log10(mu_ext)),
        dsmax_dV_ext,
        method="linear",
        bounds_error=False,
        fill_value=None,
    )
    coords = np.array([[np.log10(V), np.log10(mu)] for V in V_grid for mu in mu_grid])
    return interp(coords).reshape(len(V_grid), len(mu_grid))


def _compute(V_grid: np.ndarray, mu_grid: np.ndarray, a) -> dict:
    n_V, n_mu = len(V_grid), len(mu_grid)
    V_jax = jnp.asarray(V_grid, dtype=jnp.float64)
    mu_jax = jnp.asarray(mu_grid, dtype=jnp.float64)

    # Extended grids (one extra node in each direction) for numerical FD.
    V_ext = _extend_grid_log(V_grid)
    mu_ext = _extend_grid_log(mu_grid)
    V_ext_jax = jnp.asarray(V_ext, dtype=jnp.float64)
    mu_ext_jax = jnp.asarray(mu_ext, dtype=jnp.float64)

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
    smax_arg_ext = np.array(_grid2(_smax_arg)(V_ext_jax, mu_ext_jax))
    print(f"{time.perf_counter() - t0:.1f}s")

    print("MBN2014 … ", end="", flush=True)
    t0 = time.perf_counter()
    smax_mbn = np.array(_grid2(_smax_mbn)(V_jax, mu_jax))
    dsmax_dV_mbn_ex = np.array(_dV_grid2(_smax_mbn)(V_jax, mu_jax))
    smax_mbn_ext = np.array(_grid2(_smax_mbn)(V_ext_jax, mu_ext_jax))
    print(f"{time.perf_counter() - t0:.1f}s")

    dsmax_dV_arg_num = _num_grad_log_interp(smax_arg_ext, V_ext, mu_ext, V_grid, mu_grid)
    dsmax_dV_mbn_num = _num_grad_log_interp(smax_mbn_ext, V_ext, mu_ext, V_grid, mu_grid)

    # ── Parcel model ──────────────────────────────────────────────────────────
    # accom is a physical constant, same for every model build.
    from pyrcel import constants as _c

    _accom = _c.ac

    # Scale t_end with V so the solver reaches ~1500 m altitude for every grid
    # point without running on past S_max. A fixed global t_end (e.g. 15000 s)
    # forces high-V runs to integrate for tens of kilometres of ascent, exhausting
    # max_steps and producing NaN before the output grid is complete.
    # N_TS is fixed so all ts arrays share the same shape → JIT compiles once.
    #
    # Note: t_end(V) ∝ 1/V, so ts varies with V. This does NOT bias the adjoint.
    # By the envelope theorem, d(S_max)/d(t_end) = 0 whenever the supersaturation
    # peak is strictly interior to [t0, t_end] — which _ts_for_V guarantees — so
    # the ∂S_max/∂t_end · dt_end/dV coupling term vanishes identically.
    _N_TS = 600

    def _ts_for_V(V_val: float) -> jnp.ndarray:
        t_end = max(200.0, min(1500.0 / float(V_val), 15000.0))
        return jnp.linspace(0.0, t_end, _N_TS)

    # ts is passed as an explicit arg so it can vary per grid point while the
    # JIT-compiled kernel is reused (shape is always (_N_TS,)).
    def _smax_parcel(V_val, y0, r_drys, Nis, kappas_arr, ts):
        return max_supersaturation(y0, (r_drys, Nis, kappas_arr, _accom, ConstantV(V_val)), ts)

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
            ts_j = _ts_for_V(V_val)
            cold = done == 0

            label = f"  [{done + 1:03d}/{total}] V={V_val:.2f} m/s  μ={mu_val:.4f} µm"
            if cold:
                label += "  (cold JIT — compiling forward + adjoint …)"
            print(label, end="", flush=True)
            t0 = time.perf_counter()

            try:
                s = float(jax.block_until_ready(_fwd(V_j, y0_j, rd_j, ni_j, ka_j, ts_j)))
                g = float(jax.block_until_ready(_grad(V_j, y0_j, rd_j, ni_j, ka_j, ts_j)))
            except Exception as exc:
                print(f"  ✗ ({exc.__class__.__name__})")
                s, g = np.nan, np.nan
            else:
                print(f"  {time.perf_counter() - t0:.1f}s")

            smax_parcel[i_V, i_mu] = s
            dsmax_dV_parcel_ex[i_V, i_mu] = g
            done += 1

    # ── Extended-grid forward pass for numerical FD ───────────────────────────
    # The original n_V×n_mu nodes sit exactly at [1:-1, 1:-1] in the extended
    # grid. We only need to evaluate the parcel model on the outer ring (44 pts)
    # and reuse smax_parcel for the interior, giving central differences at all
    # original nodes (no boundary one-sided differences).
    smax_parcel_ext = np.full((n_V + 2, n_mu + 2), np.nan)
    smax_parcel_ext[1:-1, 1:-1] = smax_parcel

    # Build list of outer-ring (ext_i, ext_j, V_val, mu_val) to evaluate.
    outer: list[tuple[int, int, float, float]] = []
    for j, mu_val in enumerate(mu_ext):  # top + bottom rows
        outer.append((0, j, float(V_ext[0]), float(mu_val)))
        outer.append((n_V + 1, j, float(V_ext[-1]), float(mu_val)))
    for i, V_val in enumerate(V_grid):  # left + right columns (no corners)
        outer.append((i + 1, 0, float(V_val), float(mu_ext[0])))
        outer.append((i + 1, n_mu + 1, float(V_val), float(mu_ext[-1])))

    print(f"\nExtended grid: {len(outer)} additional forward evaluations (no adjoint)")
    for k, (ext_i, ext_j, V_val, mu_val) in enumerate(outer):
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
        V_j = jnp.float64(V_val)
        ts_j = _ts_for_V(V_val)
        print(f"  ext [{k + 1:02d}/{len(outer)}] V={V_val:.3f}  μ={mu_val:.5f}", end="", flush=True)
        t0 = time.perf_counter()
        try:
            s = float(jax.block_until_ready(_fwd(V_j, y0_j, rd_j, ni_j, ka_j, ts_j)))
        except Exception as exc:
            print(f"  ✗ ({exc.__class__.__name__})")
            s = np.nan
        else:
            print(f"  {time.perf_counter() - t0:.1f}s")
        smax_parcel_ext[ext_i, ext_j] = s

    dsmax_dV_parcel_num = _num_grad_log_interp(smax_parcel_ext, V_ext, mu_ext, V_grid, mu_grid)

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
        from matplotlib.gridspec import GridSpec
        from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
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

    # --- Global colour limits -------------------------------------------------

    all_grad = np.concatenate(
        [arr[np.isfinite(arr)] for rdef in row_defs for arr in (rdef["exact"], rdef["num"])]
    )
    grad_vmax = float(np.nanpercentile(np.abs(all_grad), 98))
    grad_norm = mcolors.Normalize(vmin=0, vmax=grad_vmax)

    all_err = np.concatenate(
        [
            (rdef["exact"] - rdef["num"])[np.isfinite(rdef["exact"]) & np.isfinite(rdef["num"])]
            for rdef in row_defs
        ]
    )
    err_vlim = max(float(np.nanpercentile(np.abs(all_err), 98)), 1e-10)
    err_norm = mcolors.TwoSlopeNorm(vmin=-err_vlim, vcenter=0.0, vmax=err_vlim)

    # --- Font sizes -----------------------------------------------------------
    FS_TICK = 11
    FS_LABEL = 13
    FS_TITLE = 11

    # --- Layout: 3×4 [exact | num | gap | error]; colorbars below as h-bars --
    fig = plt.figure(figsize=(13, 10))
    gs = GridSpec(
        3,
        4,
        figure=fig,
        hspace=0.28,
        wspace=0.12,
        width_ratios=[10, 10, 1.5, 8],
        left=0.09,
        right=0.95,
        top=0.93,
        bottom=0.22,
    )
    plot_axes = np.array([[fig.add_subplot(gs[i, j]) for j in (0, 1, 3)] for i in range(3)])

    # Horizontal colorbars: derive position from the numerical column bbox.
    pos_num = gs[0, 1].get_position(fig)
    cbar_cx = (pos_num.x0 + pos_num.x1) / 2
    cbar_w = (pos_num.x1 - pos_num.x0) * 1.4  # column width + 40%
    cbar_x0 = cbar_cx - cbar_w / 2
    cbar_h = 0.028

    cax_grad = fig.add_axes([cbar_x0, 0.13, cbar_w, cbar_h])  # above
    cax_err = fig.add_axes([cbar_x0, 0.05, cbar_w, cbar_h])  # below

    col_titles = [
        "Exact adjoint  (jax.grad)",
        "Numerical  (central diff, extended grid)",
        r"Error  (exact $-$ numerical)",
    ]

    # Tick positions: explicit decimals, minor ticks suppressed.
    x_ticks = [
        t for t in [0.1, 0.2, 0.5, 1, 2, 3] if V_grid.min() * 0.99 <= t <= V_grid.max() * 1.01
    ]
    y_ticks = [
        t
        for t in [0.02, 0.03, 0.05, 0.1, 0.2, 0.3]
        if mu_grid.min() * 0.99 <= t <= mu_grid.max() * 1.01
    ]
    dec_fmt = FuncFormatter(lambda x, _: f"{x:g}")

    for i_row, rdef in enumerate(row_defs):
        exact = rdef["exact"]
        num = rdef["num"]
        err = np.where(np.isfinite(exact) & np.isfinite(num), exact - num, np.nan)

        for i_col, (arr, norm, cmap) in enumerate(
            [(exact, grad_norm, "plasma"), (num, grad_norm, "plasma"), (err, err_norm, "RdBu_r")]
        ):
            ax = plot_axes[i_row, i_col]
            ax.pcolormesh(V_grid, mu_grid, arr.T, cmap=cmap, norm=norm, shading="auto")
            ax.set_xscale("log")
            ax.set_yscale("log")

            ax.xaxis.set_major_locator(FixedLocator(x_ticks))
            ax.xaxis.set_minor_locator(NullLocator())
            ax.xaxis.set_major_formatter(dec_fmt)
            ax.yaxis.set_minor_locator(NullLocator())
            ax.tick_params(labelsize=FS_TICK)

            if i_row == 0:
                ax.set_title(col_titles[i_col], fontsize=FS_TITLE, pad=6)
            if i_row < 2:
                ax.tick_params(axis="x", labelbottom=False)
            else:
                ax.set_xlabel("$V$ (m s$^{-1}$)", fontsize=FS_LABEL, labelpad=4)

            if i_col == 0:
                ax.yaxis.set_major_locator(FixedLocator(y_ticks))
                ax.yaxis.set_major_formatter(dec_fmt)
                ax.set_ylabel(rdef["label"] + "\n$\\mu$ (µm)", fontsize=FS_LABEL, labelpad=4)
            else:
                ax.tick_params(axis="y", labelleft=False)

    # --- Horizontal shared colorbars at bottom --------------------------------
    sm_grad = plt.cm.ScalarMappable(norm=grad_norm, cmap="plasma")
    cbar_grad = fig.colorbar(sm_grad, cax=cax_grad, orientation="horizontal")
    cbar_grad.set_label(
        r"$\partial S_{\mathrm{max}}/\partial V$  (% (m s$^{-1}$)$^{-1}$)",
        fontsize=FS_LABEL - 1,
        labelpad=5,
    )
    cbar_grad.ax.tick_params(labelsize=FS_TICK)

    sm_err = plt.cm.ScalarMappable(norm=err_norm, cmap="RdBu_r")
    cbar_err = fig.colorbar(sm_err, cax=cax_err, orientation="horizontal")
    cbar_err.set_label(
        r"$\Delta\,\partial S_{\mathrm{max}}/\partial V$  (% (m s$^{-1}$)$^{-1}$)",
        fontsize=FS_LABEL - 1,
        labelpad=5,
    )
    cbar_err.ax.tick_params(labelsize=FS_TICK)

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved figure to {path}")
    plt.close(fig)


if __name__ == "__main__":
    sensitivity_sweep()
