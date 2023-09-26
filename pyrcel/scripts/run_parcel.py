#!/usr/bin/env python
"""
CLI interface to run parcel model simulation.

"""
import os
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import yaml

import pyrcel as pm
import pyrcel.util

parser = ArgumentParser(
    description=__doc__, formatter_class=RawDescriptionHelpFormatter
)
parser.add_argument(
    "namelist",
    type=str,
    metavar="config.yml",
    help="YAML namelist controlling simulation configuration",
)

DIST_MAP = {
    "lognormal": pm.Lognorm,
}


def run_parcel():
    # Read command-line arguments
    args = parser.parse_args()

    # Convert the namelist file into an in-memory dictionary
    try:
        print("Attempting to read simulation namelist {}".format(args.namelist))
        with open(args.namelist, "rb") as f:
            y = yaml.safe_load(f)
    except IOError:
        print("Couldn't read file {}".format(args.namelist))
        sys.exit(0)

    # Create the aerosol
    aerosol_modes = []
    print("Constructing aerosol modes")
    for i, aerosol_params in enumerate(y["initial_aerosol"], start=1):
        ap = aerosol_params

        dist = DIST_MAP[ap["distribution"]](**ap["distribution_args"])

        aer = pm.AerosolSpecies(ap["name"], dist, kappa=ap["kappa"], bins=ap["bins"])
        print("   {:2d})".format(i), aer)

        aerosol_modes.append(aer)

    # Set up the model
    ic = y["initial_conditions"]
    print("Initializing model")
    try:
        model = pm.ParcelModel(
            aerosol_modes,
            V=ic["updraft_speed"],
            T0=ic["temperature"],
            S0=-1.0 * (1.0 - ic["relative_humidity"]),
            P0=ic["pressure"],
            console=True,
            truncate_aerosols=True,
        )
    except pyrcel.util.ParcelModelError:
        print("Something went wrong setting up the model")
        sys.exit(0)

    # Run the model
    mc = y["model_control"]
    print("Beginning simulation")
    try:
        par_out, aer_out = model.run(
            max_steps=2000,
            solver="cvode",
            output_fmt="dataframes",
            # terminate=True,
            # terminate_depth=10.,
            **mc,
        )
    except pyrcel.util.ParcelModelError:
        print("Something went wrong during model run")
        sys.exit(0)

    # Output
    ec = y["experiment_control"]

    Smax = par_out["S"].max()
    T_fin = par_out["T"].iloc[-1]

    # Make output directory if it doesn't exist
    if not os.path.exists(ec["output_dir"]):
        os.makedirs(ec["output_dir"])

    out_file = os.path.join(ec["output_dir"], ec["name"]) + ".nc"
    try:
        print("Trying to save output to {}".format(out_file))
        pm.output.write_parcel_output(out_file, parcel=model)
    except (IOError, RuntimeError):
        print("Something went wrong saving to {}".format(out_file))
        sys.exit(0)

    # Succesful completion
    print("Done!")


if __name__ == "__main__":
    run_parcel()
