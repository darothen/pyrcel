#!/usr/bin/env python
"""Propagate a Gaussian updraft-velocity distribution through the parcel model.

Samples ``n`` vertical velocities from ``Normal(mean, std)`` and runs a single
``jax.vmap``-ed ensemble of parcel-model simulations (one compiled, batched call) to
estimate the resulting distributions of peak supersaturation ``S_max`` and activated
droplet number concentration ``N_act``.

This showcases the headline benefit of the JAX/diffrax backend: vectorized ensembles for
uncertainty propagation / sensitivity studies that the numba+CVode stack cannot do in one
call. Each ensemble member integrates to its *own* supersaturation maximum via the
``dS/dt`` event, so widely varying velocities are handled correctly.

Usage
-----
    python examples/jax/updraft_ensemble.py --mean 0.5 --std 0.2 --n 512
    python examples/jax/updraft_ensemble.py --mean 1.0 --std 0.4 --n 1024 \
        --out output/ensemble.nc --plot output/ensemble.png
"""

from __future__ import annotations

import argparse
import os
import time

import numpy as np

import pyrcel as pm


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mean", type=float, default=0.5, help="mean updraft, m/s")
    p.add_argument("--std", type=float, default=0.2, help="updraft std dev, m/s")
    p.add_argument("--n", type=int, default=512, help="ensemble size")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--v-min", type=float, default=0.02, help="clip floor for sampled V, m/s")
    p.add_argument("--mu", type=float, default=0.05, help="lognormal geometric mean radius, micron")
    p.add_argument("--sigma", type=float, default=2.0, help="lognormal geometric std dev")
    p.add_argument("--N", type=float, default=1000.0, help="aerosol number, cm^-3")
    p.add_argument("--kappa", type=float, default=0.54, help="hygroscopicity")
    p.add_argument("--bins", type=int, default=50, help="number of size bins")
    p.add_argument("--T0", type=float, default=283.15)
    p.add_argument("--P0", type=float, default=85000.0)
    p.add_argument("--S0", type=float, default=-0.02)
    p.add_argument(
        "--z-cap", type=float, default=400.0, help="target height for the t_end heuristic, m"
    )
    p.add_argument("--out", type=str, default=None, help="optional NetCDF output path")
    p.add_argument("--plot", type=str, default=None, help="optional PNG histogram path")
    a = p.parse_args()

    aerosol = pm.AerosolSpecies(
        "sulfate", pm.Lognorm(mu=a.mu, sigma=a.sigma, N=a.N), kappa=a.kappa, bins=a.bins
    )

    print(
        f"Running {a.n}-member updraft ensemble: V ~ N({a.mean}, {a.std}^2) m/s "
        f"(clipped at {a.v_min})"
    )
    t0 = time.perf_counter()
    res = pm.run_updraft_ensemble(
        [aerosol],
        T0=a.T0,
        S0=a.S0,
        P0=a.P0,
        mean=a.mean,
        std=a.std,
        n=a.n,
        seed=a.seed,
        v_min=a.v_min,
        z_cap=a.z_cap,
    )
    elapsed = time.perf_counter() - t0

    V = res["V"]
    smax_pct = res["S_max"] * 100.0
    nact_cm3 = res["N_act"] * 1e-6
    activated = res["activated"]
    n_reached = int(np.sum(activated))

    print(
        f"  done in {elapsed:.1f}s (compile + solve) | t_end={res['t_end']:.0f}s | "
        f"{n_reached}/{a.n} members reached an interior S_max"
    )
    _print_stats("V [m/s]", V)
    _print_stats("S_max [%]", smax_pct[activated])
    _print_stats("N_act [cm^-3]", nact_cm3[activated])

    if a.out is not None:
        _save_netcdf(res, a, a.out)
    if a.plot is not None:
        _plot(V[activated], smax_pct[activated], nact_cm3[activated], a.plot)

    return 0


def _print_stats(label: str, x: np.ndarray) -> None:
    if x.size == 0:
        print(f"  {label:>16}: (no samples)")
        return
    pct = np.percentile(x, [5, 25, 50, 75, 95])
    print(
        f"  {label:>16}: mean={x.mean():.4g}  std={x.std():.4g}  "
        f"p05={pct[0]:.4g}  p50={pct[2]:.4g}  p95={pct[4]:.4g}"
    )


def _save_netcdf(res, a, path: str) -> None:
    import xarray as xr

    ds = xr.Dataset(
        {
            "V": (("member",), res["V"], {"units": "m/s", "long_name": "updraft speed"}),
            "S_max": (
                ("member",),
                res["S_max"] * 100.0,
                {"units": "%", "long_name": "peak supersaturation"},
            ),
            "N_act": (
                ("member",),
                res["N_act"] * 1e-6,
                {"units": "cm-3", "long_name": "activated droplet number"},
            ),
            "T_smax": (
                ("member",),
                res["T_smax"],
                {"units": "K", "long_name": "temperature at S_max"},
            ),
            "activated": (
                ("member",),
                res["activated"],
                {"long_name": "reached interior S_max before t_end"},
            ),
        },
        attrs={
            "source": "pyrcel JAX/diffrax updraft ensemble",
            "V_mean": a.mean,
            "V_std": a.std,
            "v_min": a.v_min,
            "T0": a.T0,
            "S0": a.S0,
            "P0": a.P0,
            "kappa": a.kappa,
            "t_end": res["t_end"],
        },
    )
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    ds.to_netcdf(path)
    print(f"  saved ensemble to {path}")


def _plot(V, smax_pct, nact_cm3, path: str) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot")
        return

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    axes[0].hist(smax_pct, bins=40, color="C0", alpha=0.85)
    axes[0].set_xlabel("S_max [%]")
    axes[0].set_ylabel("count")
    axes[1].hist(nact_cm3, bins=40, color="C2", alpha=0.85)
    axes[1].set_xlabel("N_act [cm$^{-3}$]")
    axes[2].scatter(V, nact_cm3, s=8, alpha=0.4, color="C3")
    axes[2].set_xlabel("V [m/s]")
    axes[2].set_ylabel("N_act [cm$^{-3}$]")
    fig.suptitle("Updraft-velocity ensemble: output distributions")
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=120)
    print(f"  saved figure to {path}")


if __name__ == "__main__":
    raise SystemExit(main())
