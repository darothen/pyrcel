from aerosol import AerosolSpecies
from lognorm import Lognorm
from activation import act_fraction
from driver import iterate_runs
from ext.ghan.activation import explicit as ghan_explicit
from ext.ghan.activation import parameterization as ghan_param

import numpy as np

P0 = 85000. # Pressure, Pa
T0 = 283. # Temperature, K
S0 = 0.0 # Supersaturation. 1-RH from wv term
V = 2.5 # m/s

mu = 0.21
sigma =  1.5
kappa = 0.54
N = 500

## Add aerosols to the parcel
aerosol1 = AerosolSpecies('(NH4)2SO4', Lognorm(mu=mu, sigma=sigma, N=N),
                          bins=200, kappa=kappa, rho=1760.)
initial_aerosols = [aerosol1, ]

aer_species = [a.species for a in initial_aerosols]
aer_dict = dict()
for aerosol in initial_aerosols:
    aer_dict[aerosol.species] = aerosol

dt = 0.05 # seconds
t_end = 500./V

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
(p, a), ghan_Smax, nenes_Smax = iterate_runs(V, initial_aerosols, T0, P0,
                                             output="dataframes", 
                                             dt=dt, t_end=t_end)

Smax = p.S.max()
tmax = p.S.argmax()
a = a[aerosol1.species]
rs = np.array(a.ix[tmax].tolist())
eq, kn = act_fraction(Smax, T0, rs, kappa, aerosol1.r_drys, aerosol1.Nis)

_, ghan_exp_Smax = ghan_explicit(initial_aerosols, V, T0, S0, P0)
_, ghan_param_Smax = ghan_param(initial_aerosols, V, T0, S0, P0)

## print activation stuff

print "               Smax"
print "       Nenes", nenes_Smax
print "        Ghan", ghan_Smax
#print "  Ming", ming_Smax
print "      Parcel", Smax
print "--"*40
print " Ghan parcel", ghan_exp_Smax
print "  Ghan param", ghan_param_Smax
print "--"*40
print eq, kn

import matplotlib.pyplot as plt
def quick_plot():
    plt.figure(1)
    plt_meta = p.plot("z", "S", zorder=3)
    c = plt_meta.lines[-1].get_color()
    plt.vlines(p.z.ix[tmax], 0, Smax, color=c, linestyle='dashed')
    plt.hlines(Smax, 0, p.z.ix[tmax], color=c, linestyle='dashed')
    plt.ylim(0)

def activation_plot():
    eqs, kns = [], []

    rising = a.ix[:tmax]
    print "len = ",len(rising)

    if len(rising) > 20:
        step = len(rising)/20
        print "step = ", step
    else:
        step = 1
    #step = 1

    for t in a.index[::step]:
        print t,
        if t > tmax: break

        rs = np.array(a.ix[t].tolist())
        Smax = p.S.ix[t]
        eq, kn = act_fraction(Smax, T0, rs, kappa, 
                              aerosol1.r_drys, aerosol1.Nis)
        eqs.append(eq)
        kns.append(kn)
        print eq, kn
    print p.z.ix[:tmax:step].shape, len(eqs)

    plt.figure(2)
    plt.clf()
    plt.plot(p.z.ix[:tmax:step], eqs, color='k', label="Eq")
    plt.plot(p.z.ix[:tmax:step], kns, color='r', label="Kn")
    plt.legend(loc="upper left")
    plt.ylim(0, 1)

quick_plot()
plt.show(block=True)


