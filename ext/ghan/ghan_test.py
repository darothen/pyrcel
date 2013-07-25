from parcel_model import Lognorm
from parcel_model import ParcelModel, AerosolSpecies

import pandas
from activation import explicit, parameterization

import numpy as np

from collections import OrderedDict

P0 = 80000. # Pressure, Pa
T0 = 283.15 # Temperature, K
S0 = -0.00 # Supersaturation. 1-RH from wv term
#V = 0.7076312079 # m/s
V = 2.5

aerosol1 = AerosolSpecies('(NH4)2SO4', Lognorm(mu=0.02494149518, sigma=1.301854178, N=2299.298741 ),
                          bins=200, kappa=0.4983795899)
aerosol1.rho = 1760.

initial_aerosols = [aerosol1, ]

aer_species = [a.species for a in initial_aerosols]
aer_dict = OrderedDict()
for aerosol in initial_aerosols:
    aer_dict[aerosol.species] = aerosol

parameterization(initial_aerosols, V, T0, S0, P0)
explicit(initial_aerosols, V, T0, S0, P0)
