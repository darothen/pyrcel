#!/usr/bin/env python
"""Basic single parcel-model run on the JAX/diffrax backend.

Two aerosol modes (accumulation + Aitken) illustrate how different size populations
partially activate or remain as haze droplets.

Usage
-----
    python examples/basic_run.py
    python examples/basic_run.py --V 0.5 --out output/basic_run.nc --plot output/basic_run.png
"""

from __future__ import annotations

import argparse
import os

import pyrcel as pm


def basic_run() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--V", type=float, default=1.0, help="updraft speed (m/s)")
    p.add_argument("--T0", type=float, default=283.0, help="initial temperature (K)")
    p.add_argument("--P0", type=float, default=85000.0, help="initial pressure (Pa)")
    p.add_argument("--S0", type=float, default=-0.02, help="initial supersaturation")
    p.add_argument("--t-end", type=float, default=100.0, help="max integration time (s)")
    p.add_argument("--out", type=str, default="output/jax_basic_run.nc")
    p.add_argument("--plot", type=str, default=None, help="optional PNG figure path")
    a = p.parse_args()

    accumulation = pm.AerosolSpecies(
        "accumulation",
        pm.Lognorm(mu=0.05, sigma=2.0, N=1000.0),
        kappa=0.54,
        bins=50,
    )
    aitken = pm.AerosolSpecies(
        "aitken",
        pm.Lognorm(mu=0.01, sigma=1.6, N=2000.0),
        kappa=0.3,
        bins=50,
    )

    model = pm.ParcelModel([accumulation, aitken], V=a.V, T0=a.T0, S0=a.S0, P0=a.P0, console=True)
    out = model.run(a.t_end, output_dt=1.0, terminate=False, progress=True)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    out.to_netcdf(a.out)
    print(f"Saved trajectory to {a.out}")

    if a.plot is not None:
        _plot(out, model, a.plot)

    return 0


def _plot(out: pm.ModelOutput, model: pm.ParcelModel, path: str) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot")
        return

    import numpy as np
    from matplotlib.lines import Line2D

    from pyrcel.equilibrate import kohler_crit_approx

    ds = out.to_xarray()
    s = model.summary()
    height = ds["height"].values
    t_smax = float(ds["t_smax"])
    z_smax = float(ds["height"].sel(time=t_smax, method="nearest"))

    S_COLOR = "#2271B5"
    T_COLOR = "#C44E52"
    SMAX_KW = dict(color="#C9A227", lw=1.0, ls="--", zorder=2)

    # Species colors — distinct from left-panel blue/red
    ACC_COLOR = "#2A9D8F"  # teal
    AIT_COLOR = "#E76F51"  # burnt orange

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(11, 3.7))
    fig.subplots_adjust(wspace=0.40)

    # ── left panel: S (bottom axis) + T (top axis), shared height y-axis ──────
    ax_T = ax_left.twiny()

    ax_left.plot(ds["S"].values, height, color=S_COLOR, lw=1.8)
    ax_left.axvline(0.0, color="0.78", lw=0.8, ls=":")
    ax_left.axhline(z_smax, **SMAX_KW)
    ax_left.set_xlim(0.0, 0.40)
    ax_left.set_xticks(np.arange(0.0, 0.41, 0.05))
    ax_left.set_ylim(0, height.max() + 2)
    ax_left.set_xlabel("Supersaturation (%)", color=S_COLOR, labelpad=6)
    ax_left.tick_params(axis="x", colors=S_COLOR)
    ax_left.spines["bottom"].set_edgecolor(S_COLOR)
    ax_left.spines["top"].set_visible(False)
    ax_left.spines["right"].set_visible(False)
    ax_left.set_ylabel("Height (m)")
    ax_left.text(
        0.03,
        z_smax,
        f"$S_{{\\max}} = {s['S_max'] * 100:.3f}\\%$",
        transform=ax_left.get_yaxis_transform(),
        fontsize=8,
        color="#C9A227",
        va="bottom",
    )

    ax_T.plot(ds["T"].values, height, color=T_COLOR, lw=1.8)
    ax_T.set_xlim(282.0, 283.0)
    ax_T.set_xticks([282.0, 282.2, 282.4, 282.6, 282.8, 283.0])
    ax_T.set_xlabel("Temperature (K)", color=T_COLOR, labelpad=6)
    ax_T.tick_params(axis="x", colors=T_COLOR)
    ax_T.spines["top"].set_edgecolor(T_COLOR)
    ax_T.spines["bottom"].set_visible(False)
    ax_T.spines["right"].set_visible(False)
    ax_T.spines["left"].set_visible(False)

    # ── right panel: per-bin radius traces ─────────────────────────────────────
    ax_right.axhline(z_smax, **SMAX_KW)
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.set_ylim(0, height.max() + 2)

    N_TRACES = 22
    legend_handles: list = []
    per_sp = s["per_species"]

    for sp_idx, aer in enumerate(out.aerosols):
        sp = aer.species
        rdry_m = np.asarray(aer.r_drys)
        rwet = ds[f"{sp}_size"].values  # µm, shape (time, nr)
        nr = aer.nr
        color = ACC_COLOR if sp_idx == 0 else AIT_COLOR

        # Per-bin activation: a bin activates if its critical supersaturation
        # is below the actual S_max.  Comparing s_crit < S_max is correct for
        # all sizes; comparing final wet radius to r_crit fails for large dry
        # particles whose r_crit >> any wet radius reached in the simulation.
        T_final = float(ds["T"].values[-1])
        _, s_crit = kohler_crit_approx(T_final, rdry_m, aer.kappa)
        activated = np.asarray(s_crit) < s["S_max"]

        sample = np.unique(np.round(np.linspace(0, nr - 1, N_TRACES)).astype(int))

        for i in sample:
            if activated[i]:
                ax_right.plot(rwet[:, i], height, color=color, lw=2.0, solid_capstyle="round")
            else:
                ax_right.plot(rwet[:, i], height, color=color, lw=0.9, ls="--", alpha=0.60)

        frac = per_sp[sp_idx]["eq_act_frac"]
        frac_str = "< 0.01" if 0 < frac < 0.01 else f"{frac:.2f}"
        legend_handles.append(
            Line2D(
                [0], [0], color=color, lw=2.0, label=f"{sp}  ($f_{{\\mathrm{{act}}}} = {frac_str}$)"
            )
        )

    # Style/activation proxies — interleaved so ncol=2 gives species left, style right
    legend_handles.insert(
        1,
        Line2D([0], [0], color="0.35", lw=2.0, label="activated"),
    )
    legend_handles.append(
        Line2D([0], [0], color="0.35", lw=0.9, ls="--", alpha=0.60, label="haze"),
    )

    ax_right.legend(
        handles=legend_handles,
        ncol=2,
        loc="lower left",
        bbox_to_anchor=(0, 1.01),
        fontsize=8.5,
        frameon=False,
        columnspacing=1.2,
        handlelength=1.8,
    )

    ax_right.set_xscale("log")
    ax_right.set_xlim(0.01, 10.0)
    ax_right.set_xticks([0.01, 0.1, 1, 10])
    ax_right.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:g}"))
    ax_right.set_xlabel("Droplet radius, µm", labelpad=6)

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved figure to {path}")


if __name__ == "__main__":
    basic_run()
