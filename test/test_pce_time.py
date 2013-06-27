from lognorm import Lognorm
from parcel import AerosolSpecies
from pce_params import pce_deg1, pce_deg2, pce_deg3, pce_deg4
from activation import arg2000, fn2005

P0 = 80000. # Pressure, Pa
T0 = 283.15 # Temperature, K
S0 = -0.00 # Supersaturation. 1-RH from wv term
V = 0.7076312079 # m/s

aerosol1 = AerosolSpecies('(NH4)2SO4', Lognorm(mu=0.02494149518, sigma=1.301854178, N=2299.298741 ),
                          bins=200, kappa=0.4983795899)

pce1_Smax, pce1_ratio = pce_deg1(V, T0, P0, aerosol1)
pce2_Smax, pce2_ratio = pce_deg2(V, T0, P0, aerosol1)
pce3_Smax, pce3_ratio = pce_deg3(V, T0, P0, aerosol1)
pce4_Smax, pce4_ratio = pce_deg4(V, T0, P0, aerosol1)

ghan_Smax, ghan_ratio = arg2000(V, T0, P0, [aerosol1, ])

nenes_Smax, nenes_ratio = fn2005(V, T0, P0, [aerosol1, ])