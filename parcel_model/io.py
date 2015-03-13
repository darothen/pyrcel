
import pandas as pd
import xray

from thermo import es

def write_parcel_output(parcel, format='xray'):
    """ Write model output to disk.

    Wrapper for methods to write parcel model output to disk.

    """
    
    pass

def parcel_to_dataframes(parcel):
    """ Convert model simulation to dataframe. 

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

    parcel_out = pd.DataFrame( {'T':x[:,1], 'wv':x[:,2],
                               'wc':x[:,5], 'S':x[:,4], 'z':heights},
                                index=time) )
    ## Add some thermodynamic output to the parcel model dataframe
    ess = es(parcel_out['T'] - 273.15)
    parcel_out['P'] = ess*(1. + 0.622*(parcel_out['S'] + 1.)/parcel_out['wv'])

    aerosol_dfs = {}
    species_shift = 0 # increment by nr to select the next aerosol's radii
    for aerosol in parcel_out.aerosols:
        nr = aerosol.nr
        species = aerosol.species

        labels = ["r%03d" % i for i in xrange(nr)]
        radii_dict = dict()
        for i, label in enumerate(labels):
            radii_dict[label] = x[:,5+species_shift+i]

        aerosol_dfs[species] = pd.DataFrame( radii_dict,
                                             index=time  ))
        species_shift += nr

    return parcel_out, aerosol_dfs