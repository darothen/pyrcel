from parcel import ParcelModel
from aerosol import AerosolSpecies
from lognorm import Lognorm
from activation import arg2000, fn2005, act_fraction
from driver import iterate_runs
import time

P0 = 80000. # Pressure, Pa
T0 = 283. # Temperature, K
S0 = 0.0 # Supersaturation. 1-RH from wv term
V = 1.693813019475383e+00 # m/s

mu = 1.453787108519219e-02
sigma =  1.200005759431596e+00
kappa = 1.179201424192869e+00
N = 2.871871822641781e+02

## Add aerosols to the parcel
aerosol1 = AerosolSpecies('(NH4)2SO4', Lognorm(mu=mu, sigma=sigma, N=N),
                          bins=200, kappa=kappa)
initial_aerosols = [aerosol1, ]

aer_species = [a.species for a in initial_aerosols]
aer_dict = dict()
for aerosol in initial_aerosols:
    aer_dict[aerosol.species] = aerosol

t_end = 500.
dt = 0.01 # seconds

## Vanilla parcel model run
'''
pm = ParcelModel(initial_aerosols, V, T0, S0, P0, console=True)
start = time.time()
parcel, aerosols = pm.run(t_end, dt=dt, max_steps=2000, solver='cvode')
end = time.time()

ghan_Smax, ghan_ratio = arg2000(V, T0, P0, initial_aerosols)
nenes_Smax, nenes_ratio = fn2005(V, T0, P0, initial_aerosols)
'''
## Iteration straetegy
Smax, ghan_Smax, nenes_Smax = iterate_runs(V, initial_aerosols, T0, P0)


## print activation stuff


print "       Smax"
print " Nenes", nenes_Smax
print "  Ghan", ghan_Smax
#print "  Ming", ming_Smax
print "Parcel", Smax