"""
.. module:: parcel
    :synopsis: Input/Output utilities for parcel run data.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""

__docformat__ = "reStructuredText"

import os

def write_csv(parcel_data, aerosol_data, output_dir=None):
    """Write output to CSV files.

    Utilize pandas fast output procedures to write the model run output to a set of CSV files.

    :param parcel_data:
        pandas DataFrame of the parcel thermodynamic properties
    :type parcel_data: pandas.DataFrame

    :param aerosol_data:
        dictionary of pandas DataFrames with the aerosol radii at each model step
    :type aerosol_data: dictionary
    """

    if not output_dir:
        output_dir = os.getcwd()

    # Write parcel data
    parcel_data.to_csv(os.path.join(output_dir, "parcel.csv"), header=None)

    # Write aerosol data
    for species, data in aerosol_data.iteritems():
        data.to_csv(os.path.join(output_dir, "%s.csv" % species), header=None)