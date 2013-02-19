from parcel_model.lognorm import Lognorm
from parcel_model.parcel import ParcelModel, AerosolSpecies
from parcel_model.micro import kohler_crit, Rd, r_eff, activation

import pandas
import ghan_act

import numpy as np

from collections import OrderedDict

P0 = 100000. # Pressure, Pa
T0 = 279.0 # Temperature, K
S0 = -0.00 # Supersaturation. 1-RH from wv term
V = 0.1 # m/s
explicit = False

aerosol1 = AerosolSpecies('(NH4)2SO4', Lognorm(mu=0.05, sigma=2.0, N=1000.),
                          bins=200, kappa=0.7)
aerosol1.rho = 1760.

initial_aerosols = [aerosol1, ]

aer_species = [a.species for a in initial_aerosols]
aer_dict = OrderedDict()
for aerosol in initial_aerosols:
    aer_dict[aerosol.species] = aerosol

###
## Unpack the aerosol data for the Ghan code
na = [a.N for a in initial_aerosols]
sig = [a.sigma for a in initial_aerosols]
rhodry = [a.rho*1e-3 for a in initial_aerosols]
hygro = [a.kappa for a in initial_aerosols]
rad = [a.mu for a in initial_aerosols]

## Call Ghan code
fn, smax = ghan_act.activation(T0, P0/100., V*100., explicit,
                               na, sig, rhodry, hygro, rad)

print fn, smax

