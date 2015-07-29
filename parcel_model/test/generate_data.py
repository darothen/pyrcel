""" Generate test case data for future reference.
"""
from __future__ import absolute_import, print_function

from parcel_model import thermo

from itertools import product

import os, pickle
import numpy as np

REFERENCE_FN = "results.dict"

def generate_reference(overwrite=False):

    results = dict()
    results['temperatures'] = np.linspace(233, 333, 10) # K
    results['pressures'] = np.linspace(100, 1050, 10) # hPa
    results['radii'] = np.logspace(-3, 1, 10) # microns

    print("dv_cont", end=", ")
    results['dv_cont'] = \
        [ thermo.dv_cont(T, P*100) for T, P in \
            product(results['temperatures'], results['pressures']) ]
    print("(%d cases)" % len(results['dv_cont']))

    print("dv", end=", ")
    results['dv'] = \
        [ thermo.dv(T, r*1e-6, P*100) for T, r, P in \
            product(results['temperatures'], results['radii'],
                    results['pressures']) ]
    print("(%d cases)" % len(results['dv']))

    if (not os.path.exists(REFERENCE_FN)) or overwrite:
        with open(REFERENCE_FN, 'wb') as f:
            pickle.dump(results, f)

if __name__ == "__main__":

    generate_reference(True)