"""CLI interface for running a parcel-model simulation with the v2 JAX/diffrax backend.

Usage
-----
.. code-block:: bash

    run_parcel config.yml                    # write output/<name>.nc
    run_parcel config.yml -o results/run1   # explicit output path (no .nc extension needed)
    run_parcel config.yml --no-console      # suppress the setup/summary table

YAML Namelist Format
--------------------
The namelist is a YAML file with three required top-level keys:

.. code-block:: yaml

    experiment_control:
        name: "simple"          # base name for the output file
        output_dir: "output/"   # directory to write the NetCDF file into

    model_control:
        output_dt: 1.0          # output cadence (s)
        t_end: 9999.0           # maximum integration time (s)
        terminate: true         # stop shortly after S_max (recommended)
        terminate_depth: 10.0   # extra height (m) to integrate past S_max

    initial_aerosol:
        - name: sulfate
          distribution: lognormal
          distribution_args: { mu: 0.15, N: 1000.0, sigma: 1.2 }
          kappa: 0.54
          bins: 250

    initial_conditions:
        temperature: 283.15       # K
        relative_humidity: 0.95   # decimal fraction (0–1)
        pressure: 85000.0         # Pa
        updraft_speed: 0.44       # m/s

Notes
-----
- Only ``lognormal`` distributions are supported; ``distribution_args`` maps to
  :class:`~pyrcel.distributions.Lognorm` keyword arguments (``mu``, ``sigma``, ``N``).
- ``relative_humidity`` is converted to initial supersaturation via
  ``S0 = RH - 1.0`` (negative for sub-saturated parcels).
- ``solver_dt`` from the old namelist is silently ignored; the v2 diffrax backend
  uses an adaptive step controller.
- Output is written as NetCDF4 via :meth:`~pyrcel.model_output.ModelOutput.to_netcdf`.
"""

from __future__ import annotations

import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pathlib import Path

import yaml

_DIST_MAP = {
    "lognormal": "Lognorm",
}

parser = ArgumentParser(
    prog="run_parcel",
    description=__doc__,
    formatter_class=RawDescriptionHelpFormatter,
)
parser.add_argument(
    "namelist",
    type=str,
    metavar="config.yml",
    help="YAML namelist controlling simulation configuration",
)
parser.add_argument(
    "-o",
    "--output",
    type=str,
    default=None,
    metavar="PATH",
    help="Override output path (default: <output_dir>/<name>.nc from namelist)",
)
parser.add_argument(
    "--no-console",
    action="store_true",
    default=False,
    help="Suppress the setup/summary console output",
)


def run_parcel() -> None:
    args = parser.parse_args()

    # --- load namelist -------------------------------------------------------
    try:
        with open(args.namelist, "rb") as f:
            cfg = yaml.safe_load(f)
    except OSError as exc:
        print(f"error: cannot read namelist '{args.namelist}': {exc}", file=sys.stderr)
        sys.exit(1)

    # --- build aerosol modes -------------------------------------------------
    import pyrcel as pm

    aerosol_modes = []
    for ap in cfg["initial_aerosol"]:
        dist_key = ap["distribution"]
        if dist_key not in _DIST_MAP:
            print(
                f"error: unsupported distribution '{dist_key}'. "
                f"Supported: {list(_DIST_MAP)}",
                file=sys.stderr,
            )
            sys.exit(1)
        dist_cls = getattr(pm, _DIST_MAP[dist_key])
        dist = dist_cls(**ap["distribution_args"])
        aer = pm.AerosolSpecies(ap["name"], dist, kappa=ap["kappa"], bins=ap["bins"])
        aerosol_modes.append(aer)

    # --- initialise model ----------------------------------------------------
    ic = cfg["initial_conditions"]
    S0 = ic["relative_humidity"] - 1.0  # RH=0.95 → S0=-0.05

    model = pm.ParcelModel(
        aerosol_modes,
        V=ic["updraft_speed"],
        T0=ic["temperature"],
        S0=S0,
        P0=ic["pressure"],
        console=not args.no_console,
    )

    # --- run -----------------------------------------------------------------
    mc = cfg["model_control"]
    run_kwargs = {
        "t_end": mc["t_end"],
        "output_dt": mc.get("output_dt", 1.0),
        "terminate": mc.get("terminate", True),
        "terminate_depth": mc.get("terminate_depth", 10.0),
    }
    output = model.run(**run_kwargs, mode="full")

    # --- write output --------------------------------------------------------
    if args.output:
        out_path = Path(args.output)
        if out_path.suffix != ".nc":
            out_path = out_path.with_suffix(".nc")
    else:
        ec = cfg["experiment_control"]
        out_dir = Path(ec.get("output_dir", "output"))
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{ec['name']}.nc"

    output.to_netcdf(out_path)
    if not args.no_console:
        print(f"Output written to {out_path}")


if __name__ == "__main__":
    run_parcel()
