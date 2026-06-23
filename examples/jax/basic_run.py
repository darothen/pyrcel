#!/usr/bin/env python
"""Basic single parcel-model run on the JAX/diffrax backend.

Equilibrates and integrates one parcel, prints the summary, writes a NetCDF file, and
(if matplotlib is available) saves a supersaturation/temperature-vs-height figure.

Usage
-----
    python examples/jax/basic_run.py
    python examples/jax/basic_run.py --V 0.5 --N 2000 --mu 0.05 --kappa 0.54 \
        --out output/basic_run.nc --plot output/basic_run.png
"""

from __future__ import annotations

import argparse
import os

import pyrcel as pm


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--V", type=float, default=1.0, help="updraft speed, m/s")
    p.add_argument("--T0", type=float, default=283.15, help="initial temperature, K")
    p.add_argument("--P0", type=float, default=85000.0, help="initial pressure, Pa")
    p.add_argument("--S0", type=float, default=-0.02, help="initial supersaturation (0 = 100% RH)")
    p.add_argument("--mu", type=float, default=0.05, help="lognormal geometric mean radius, micron")
    p.add_argument("--sigma", type=float, default=2.0, help="lognormal geometric std dev")
    p.add_argument("--N", type=float, default=1000.0, help="aerosol number, cm^-3")
    p.add_argument("--kappa", type=float, default=0.54, help="hygroscopicity")
    p.add_argument("--bins", type=int, default=100, help="number of size bins")
    p.add_argument("--t-end", type=float, default=300.0, help="max integration time, s")
    p.add_argument("--out", type=str, default="output/jax_basic_run.nc", help="NetCDF output path")
    p.add_argument("--plot", type=str, default=None, help="optional PNG figure path")
    a = p.parse_args()

    aerosol = pm.AerosolSpecies(
        "sulfate", pm.Lognorm(mu=a.mu, sigma=a.sigma, N=a.N), kappa=a.kappa, bins=a.bins
    )
    model = pm.ParcelModelJAX(
        [aerosol], V=a.V, T0=a.T0, S0=a.S0, P0=a.P0, console=True
    )
    model.run(a.t_end, output_dt=1.0, terminate=True, terminate_depth=10.0, progress=True)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    model.save_netcdf(a.out)

    if a.plot is not None:
        _plot(model, a.plot)

    return 0


def _plot(model, path: str) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot")
        return

    ds = model.to_dataset()
    fig, (ax_s, ax_t) = plt.subplots(1, 2, figsize=(9, 5), sharey=True)
    ax_s.plot(ds["S"], ds["height"], color="C0")
    ax_s.set_xlabel("Supersaturation [%]")
    ax_s.set_ylabel("Height [m]")
    ax_s.axvline(0.0, color="0.7", lw=0.8)
    ax_t.plot(ds["T"], ds["height"], color="C3")
    ax_t.set_xlabel("Temperature [K]")
    fig.suptitle(f"Parcel run: V={model.V} m/s, S_max={model.summary()['S_max'] * 100:.3f}%")
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=120)
    print(f"Saved figure to {path}")


if __name__ == "__main__":
    raise SystemExit(main())
