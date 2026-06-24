#!/usr/bin/env python
"""Activation comparison: parcel model vs ARG2000 vs MBN2014.

Computes peak supersaturation, activated droplet number, and activated
fraction for a single-mode sulfate aerosol population using the full parcel
model (ground truth) and two analytical parameterizations.

Usage
-----
    python examples/activation_comparison.py
    python examples/activation_comparison.py --V 2.0 --mu 0.05 --N 500
    python examples/activation_comparison.py --plot output/activation_comparison.png
"""

from __future__ import annotations

import argparse

import jax
import jax.numpy as jnp

import pyrcel as pm
from pyrcel.activation import arg2000, mbn2014

jax.config.update("jax_enable_x64", True)


def activation_comparison() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--V", type=float, default=1.0, help="updraft speed (m/s)")
    p.add_argument("--T0", type=float, default=283.0, help="initial temperature (K)")
    p.add_argument("--P0", type=float, default=85000.0, help="initial pressure (Pa)")
    p.add_argument("--S0", type=float, default=-0.02, help="initial supersaturation")
    p.add_argument("--mu", type=float, default=0.05, help="median dry radius (µm)")
    p.add_argument("--sigma", type=float, default=2.0, help="geometric std dev")
    p.add_argument("--N", type=float, default=1000.0, help="number concentration (cm⁻³)")
    p.add_argument("--kappa", type=float, default=0.54, help="hygroscopicity κ")
    p.add_argument("--bins", type=int, default=100, help="size bins for parcel model")
    p.add_argument("--t-end", type=float, default=300.0, help="max integration time (s)")
    p.add_argument("--plot", type=str, default=None, metavar="PATH",
                   help="save figure to PATH (requires matplotlib)")
    a = p.parse_args()

    print(f"\nConditions:  V = {a.V} m/s | T = {a.T0} K | P = {a.P0:.0f} Pa")
    print(f"Aerosol:     μ = {a.mu} µm | σ = {a.sigma} | N = {a.N:.0f} cm⁻³ | κ = {a.kappa}")

    # ── Parcel model (ground truth) ──────────────────────────────────────────
    aerosol = pm.AerosolSpecies(
        "sulfate",
        pm.Lognorm(mu=a.mu, sigma=a.sigma, N=a.N),
        kappa=a.kappa,
        bins=a.bins,
    )
    model = pm.ParcelModel([aerosol], V=a.V, T0=a.T0, S0=a.S0, P0=a.P0)
    model.run(a.t_end, output_dt=1.0, terminate=True, progress=False)
    s = model.summary()

    sp = s["per_species"][0]
    smax_pm  = s["S_max"] * 100       # fraction → %
    N_act_pm = sp["N_act"] * 1e-6     # m⁻³ → cm⁻³
    frac_pm  = sp["eq_act_frac"]

    # ── Analytical parameterizations ─────────────────────────────────────────
    mus  = jnp.array([a.mu])
    sigs = jnp.array([a.sigma])
    Ns   = jnp.array([a.N])
    kaps = jnp.array([a.kappa])

    smax_a, Nact_a, frac_a = arg2000(a.V, a.T0, a.P0, mus, sigs, Ns, kaps)
    smax_m, Nact_m, frac_m = mbn2014(a.V, a.T0, a.P0, mus, sigs, Ns, kaps)

    results = [
        ("Parcel model", smax_pm,              N_act_pm,           frac_pm),
        ("ARG2000",      float(smax_a) * 100,  float(Nact_a[0]),   float(frac_a[0])),
        ("MBN2014",      float(smax_m) * 100,  float(Nact_m[0]),   float(frac_m[0])),
    ]

    # ── Table ────────────────────────────────────────────────────────────────
    W = 14
    sep = "─" * (18 + W * 3 + 6)
    print(f"\n{sep}")
    print(f"  {'Method':<16}  {'Smax (%)':>{W}}  {'N_act (cm⁻³)':>{W}}  {'f_act':>{W}}")
    print(f"  {'─'*16}  {'─'*W}  {'─'*W}  {'─'*W}")
    for label, smax, nact, frac in results:
        print(f"  {label:<16}  {smax:>{W}.4f}  {nact:>{W}.1f}  {frac:>{W}.4f}")
    print(sep)

    print("\n  Relative error vs parcel model:")
    for label, smax, nact, frac in results[1:]:
        ds = (smax - smax_pm) / smax_pm * 100 if smax_pm else float("nan")
        dn = (nact - N_act_pm) / N_act_pm * 100 if N_act_pm else float("nan")
        print(f"    {label}:  ΔSmax = {ds:+.1f}%,  ΔN_act = {dn:+.1f}%")

    if a.plot:
        _plot(results, a.plot)

    return 0


def _plot(results: list, path: str) -> None:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib not available; skipping plot")
        return

    labels = [r[0].replace(" ", "\n") for r in results]
    colors = ["#2A9D8F", "#E76F51", "#C9A227"]
    x = np.arange(len(labels))
    w = 0.55

    panels = [
        ([r[1] for r in results], "$S_{\\mathrm{max}}$ (%)"),
        ([r[2] for r in results], "$N_{\\mathrm{act}}$ (cm$^{-3}$)"),
        ([r[3] for r in results], "$f_{\\mathrm{act}}$"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.5))
    fig.subplots_adjust(wspace=0.45)

    for ax, (vals, ylabel) in zip(axes, panels):
        ax.bar(x, vals, width=w, color=colors, zorder=3)
        ax.axhline(vals[0], color="0.45", lw=1.0, ls="--", zorder=2,
                   label="Parcel model")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", color="0.92", zorder=1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    import os
    os.makedirs(os.path.dirname(path) if os.path.dirname(path) else ".", exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved figure to {path}")
    plt.close(fig)


if __name__ == "__main__":
    activation_comparison()
