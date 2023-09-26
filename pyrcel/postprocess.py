""" Collection of output post-processing routines.
"""
import numpy as np
import pandas as pd

from .activation import binned_activation


def simulation_activation(model, parcel_df, aerosols_panel):
    """Given the DataFrame output from a parcel model simulation, compute
    activation kinetic limitation diagnostics.

    Parameters
    ----------
    model : ParcelModel
        The ParcelModel
    parcel_df : DataFrame used to generate the results to be analyzed
        The DataFrame containing the parcel's thermodynamic trajectory
    aerosols_panel : Panel
        A Panel collection of DataFrames containing the aerosol size evolution

    Returns
    -------
    act_stats : DataFrame
        A DataFrame containing the activation statistics

    """

    initial_row = parcel_df.iloc[0]
    Smax_i, T_i = initial_row["S"], initial_row["T"]

    acts = {"eq": [], "kn": [], "alpha": [], "phi": []}

    initial_aerosols = model.aerosols
    N_all_modes = np.sum([aer.total_N for aer in initial_aerosols])
    N_fracs = {aer.species: aer.total_N / N_all_modes for aer in initial_aerosols}
    for i in range(len(parcel_df)):
        row_par = parcel_df.iloc[i]
        rows_aer = {key: aerosols_panel[key].iloc[i] for key in aerosols_panel}

        # Update thermo
        T_i = row_par["T"]
        if row_par["S"] > Smax_i:
            Smax_i = row_par["S"]

        eq_tot, kn_tot, alpha_tot, phi_tot = 0.0, 0.0, 0.0, 0.0
        for aerosol in initial_aerosols:
            N_frac = N_fracs[aerosol.species]
            rs = rows_aer[aerosol.species]

            eq, kn, alpha, phi = binned_activation(Smax_i, T_i, rs, aerosol)
            eq_tot += eq * N_frac
            kn_tot += kn * N_frac
            alpha_tot += alpha * N_frac
            phi_tot += phi * N_frac

        acts["kn"].append(kn_tot)
        acts["eq"].append(eq_tot)
        acts["alpha"].append(alpha_tot)
        acts["phi"].append(phi_tot)
    acts_total = pd.DataFrame(acts, index=parcel_df.index)

    return acts_total
