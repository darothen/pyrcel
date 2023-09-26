import os.path
import pickle
from datetime import datetime as ddt

import numpy as np
import pandas as pd
import xarray as xr

from . import __version__ as ver
from . import constants as c
from .thermo import rho_air
from .util import ParcelModelError

#: Acceptable output formats
OUTPUT_FORMATS = ["nc", "obj", "csv"]


def get_timestamp(fmt="%m%d%y_%H%M%S"):
    """Get current timestamp in MMDDYY_hhmmss format."""

    current_time = ddt.now()
    timestamp = current_time.strftime("%m%d%y_%H%M%S")

    return timestamp


def write_parcel_output(
    filename=None,
    format=None,
    parcel=None,
    parcel_df=None,
    aerosol_dfs=None,
    other_dfs=None,
):
    """Write model output to disk.

    Wrapper for methods to write parcel model output to disk.

    Parameters
    ----------
    filename : str
        Full filename to write output; if not supplied, will default
        to the current timestamp
    format : str
        Format to use from ``OUTPUT_FORMATS``; must be supplied if no
        filename is provided
    parcel : ParcelModel
        A ParcelModel which has already been integrated at least once
    parcel_df : DataFrame
        Model thermodynamic history
    aerosol_dfs : Panel
        Aerosol size history
    other_dfs : list of DataFrames
        Additional DataFrames to include in output; must have the same index
        as the parcel's results when transformed to a DataFrame!

    """

    if not filename:
        if not format:
            raise ParcelModelError("Must supply either a filename or format.")
        if format not in OUTPUT_FORMATS:
            raise ParcelModelError("Please supply a format from %r" % OUTPUT_FORMATS)
        basename = get_timestamp()
        extension = format
    else:
        basename, extension = os.path.splitext(filename)
        extension = extension[1:]  # strip '.'
        if extension not in OUTPUT_FORMATS:
            extension = format = "obj"
        else:
            format = extension

    if parcel.console:
        print()
        print(
            "Saving output to %s format with base filename %s" % (extension, basename)
        )
        print()

    # filename = "%s.%s" % (basename, extension)

    # Sanity - check, either we need the dataframes themselves or
    # we need the model
    if (parcel_df is None) and (aerosol_dfs is None):
        if parcel is None:
            raise ValueError("Need to supply either dataframes or model")
        else:
            parcel_df, aerosol_dfs = parcel_to_dataframes(parcel)
    # Concatenate on the additional dataframes supplied by the user
    if other_dfs is not None:
        for df in other_dfs:
            parcel_df = pd.concat([parcel_df, df], axis=1)

    # 1) csv
    if format == "csv":
        # Write parcel data
        parcel_df.to_csv("%s_%s.%s" % (basename, "parcel", extension))

        # Write aerosol data
        for species, data in list(aerosol_dfs.items()):
            data.to_csv("%s_%s.%s" % (basename, species, extension))

    # 2) nc
    elif format == "nc":
        ## Construct xarray datastructure to write to netCDF
        ds = xr.Dataset(attrs={"Conventions": "CF-1.0", "source": "pyrcel v%s" % ver})

        ds.coords["time"] = (
            "time",
            parcel.time,
            {"units": "seconds", "long_name": "simulation time"},
        )

        ## Aerosol coordinates and basic data
        for aerosol in parcel.aerosols:
            if parcel.console:
                print(aerosol)

            nr = aerosol.nr
            r_drys = aerosol.r_drys * 1e6
            kappas = [aerosol.kappa] * nr
            Nis = aerosol.Nis * 1e-6
            species = aerosol.species

            aer_coord = "%s_bins" % species

            ds.coords[aer_coord] = (
                aer_coord,
                np.array(list(range(1, aerosol.nr + 1)), dtype=np.int32),
                {"long_name": "%s size bin number" % species},
            )
            ds["%s_rdry" % species] = (
                (aer_coord,),
                r_drys,
                {"units": "micron", "long_name": "%s bin dry radii" % species},
            )
            ds["%s_kappas" % species] = (
                (aer_coord,),
                kappas,
                {"long_name": "%s bin kappa-kohler hygroscopicity" % species},
            )
            ds["%s_Nis" % species] = (
                (aer_coord,),
                Nis,
                {
                    "units": "cm-3",
                    "long_name": "%s bin number concentration" % species,
                },
            )

            size_data = aerosol_dfs[species].values * 1e6
            ds["%s_size" % species] = (
                ("time", aer_coord),
                size_data,
                {"units": "micron", "long_name": "%s bin wet radii" % species},
            )

        ## Parcel data
        ds["S"] = (
            ("time",),
            parcel_df["S"] * 100.0,
            {"units": "%", "long_name": "Supersaturation"},
        )
        ds["T"] = (
            ("time",),
            parcel_df["T"],
            {"units": "K", "long_name": "Temperature"},
        )
        ds["P"] = (
            ("time",),
            parcel_df["P"],
            {"units": "Pa", "long_name": "Pressure"},
        )
        ds["wv"] = (
            ("time",),
            parcel_df["wv"],
            {"units": "kg/kg", "long_name": "Water vapor mixing ratio"},
        )
        ds["wc"] = (
            ("time",),
            parcel_df["wc"],
            {"units": "kg/kg", "long_name": "Liquid water mixing ratio"},
        )
        ds["wi"] = (
            ("time",),
            parcel_df["wi"],
            {"units": "kg/kg", "long_name": "Ice water mixing ratio"},
        )
        ds["height"] = (
            ("time",),
            parcel_df["z"],
            {"units": "meters", "long_name": "Parcel height above start"},
        )
        ds["rho"] = (
            ("time",),
            parcel_df["rho"],
            {"units": "kg/m3", "long_name": "Air density"},
        )

        ds["wtot"] = (
            ("time",),
            parcel_df["wv"] + parcel_df["wc"],
            {"units": "kg/kg", "long_name": "Total water mixing ratio"},
        )

        if "alpha" in parcel_df:
            ds["alpha"] = (
                ("time",),
                parcel_df["alpha"],
                {"long_name": "ratio of Nkn/Neq"},
            )
        if "phi" in parcel_df:
            ds["phi"] = (
                ("time",),
                parcel_df["phi"],
                {"long_name": "fraction of not-strictly activated drops in Nkn"},
            )
        if "eq" in parcel_df:
            ds["eq"] = (
                ("time",),
                parcel_df["eq"],
                {"long_name": "Equilibrium Kohler-activated fraction"},
            )
        if "kn" in parcel_df:
            ds["kn"] = (
                ("time",),
                parcel_df["kn"],
                {"long_name": "Kinetic activated fraction"},
            )

        ## Save to disk
        ds.to_netcdf(basename + ".nc")

    # 3) obj (pickle)
    else:
        assert parcel

        with open(basename + ".obj", "w") as f:
            pickle.dump(parcel, f)


def parcel_to_dataframes(parcel):
    """Convert model simulation to dataframe.

    Parameters
    ----------
    parcel : ParcelModel

    Returns
    -------
    DataFrame and list of DataFrames
        The first returned DataFrame will be the parcel model thermo/dynamical
        output with keys ``T``, ``wv``, ``wc``, ``S``, ``z``, and indexed by
        ``time``. The list of DataFrames shares the same ``time`` index, but
        is divided up so that the radii of each aerosol species are tracked in
        different containers.

    Raises
    ------
    ParcelModelError
        If the argument does not have any output saved

    Notes
    -----
    Will double-check that the parcel model passed into the function actually
    has output by seeing if it as the attribute ``x``, which is where such output
    would be saved.

    """

    x = parcel.x
    heights = parcel.heights
    time = parcel.time

    parcel_out = pd.DataFrame(
        {var: x[:, i] for i, var in enumerate(c.STATE_VARS)}, index=time
    )

    ## Add some thermodynamic output to the parcel model dataframe
    parcel_out["rho"] = rho_air(parcel_out["T"], parcel_out["P"], parcel_out["S"] + 1.0)

    aerosol_dfs = {}
    species_shift = 0  # increment by nr to select the next aerosol's radii
    for aerosol in parcel.aerosols:
        nr = aerosol.nr
        species = aerosol.species

        labels = ["r%03d" % i for i in range(nr)]
        radii_dict = dict()
        for i, label in enumerate(labels):
            radii_dict[label] = x[:, c.N_STATE_VARS + species_shift + i]

        aerosol_dfs[species] = pd.DataFrame(radii_dict, index=time)
        species_shift += nr

    return parcel_out, aerosol_dfs
