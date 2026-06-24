#!/usr/bin/env python
"""Parcel run with master-style live integration loop output.

Prints z / T / S after each integration chunk (``live=True``). Combine with
``console=True`` for setup tables and the post-solve activation summary.

Usage
-----
    python examples/live_run.py
    python examples/live_run.py --chunk-dt 5 --no-console
"""

from __future__ import annotations

import argparse

import pyrcel as pm


def live_run() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--chunk-dt",
        type=float,
        default=10.0,
        help="simulation-time chunk length for live rows (s)",
    )
    p.add_argument(
        "--no-console",
        action="store_true",
        help="live loop only (skip setup/summary tables)",
    )
    a = p.parse_args()

    aerosol = pm.AerosolSpecies(
        "NaCl",
        {"r_drys": [0.1, 0.25, 0.5], "Nis": [500.0, 200.0, 50.0]},
        kappa=1.2,
    )
    model = pm.ParcelModel(
        [aerosol],
        V=1.0,
        T0=283.15,
        S0=-0.02,
        P0=85000.0,
        console=not a.no_console,
    )
    model.run(
        100.0,
        output_dt=1.0,
        terminate=True,
        terminate_depth=10.0,
        live=True,
        live_chunk_dt=a.chunk_dt,
    )
    return 0


if __name__ == "__main__":
    live_run()
