from lognorm import Lognorm
from parcel import ParcelModel, AerosolSpecies
from micro import kohler_crit

from pylab import *
ion()
import pandas

import numpy as np

P0 = 95000. # Pressure, Pa
T0 = 285.2 # Temperature, K
S0 = -0.02 # Supersaturation. 1-RH from wv term
V = 0.5 # m/s

z_top = 100.0 # meters
dt = 0.01 # seconds

## Iterate over number concentration of biomass burn aerosol
#aer2Nis = np.logspace(2, 4, 15)
aer2rs = np.logspace(np.log10(0.01), np.log10(0.2), 15)
results = {}
N = 5000.0

#for N in aer2Nis:
#    print N
for r in aer2rs:
    print r
    initial_aerosols = [AerosolSpecies('sulf', Lognorm(mu=0.05, sigma=2.0, N=100.),
                              bins=200, kappa=0.6),
                        AerosolSpecies('carbon', {'r_drys': [r, ], 'Nis': [N, ]}, kappa=0.01)]
    aer_species = [a.species for a in initial_aerosols]
    aer_dict = dict()
    for aerosol in initial_aerosols:
        aer_dict[aerosol.species] = aerosol

    pm = ParcelModel(initial_aerosols, V, T0, S0, P0, console=True)
    parcel, aerosols = pm.run(z_top, dt, max_steps=2000)

    parcel = parcel.ix[parcel.index % 1 == 0]
    aerosols2 = {}
    for species in aer_species:
        aerosol = aerosols[species]
        aerosols2[species] = aerosol.ix[aerosol.index % 1 == 0]

    ##

    for species in aer_species:
        if species == 'sulf': continue
        aerosol = aerosols2[species]
        aer_meta = aer_dict[species]
        Nis = aer_meta.Nis

        Neq = []
        Nkn = []
        Nunact = []
        S_max = S0

        for S, T, i in zip(parcel.S, parcel['T'], xrange(len(parcel.S))):

            r_crits, s_crits = zip(*[kohler_crit(T, r_dry, aer_meta.kappa) for r_dry in aer_meta.r_drys])
            s_crits = np.array(s_crits)
            r_crits = np.array(r_crits)
            if S > S_max: S_max = S

            big_s =  S_max >= s_crits
            Neq.append(np.sum(Nis[big_s]))

            rstep = np.array(aerosol.ix[i])
            #active_radii = (S > s_crits) | (rstep > r_crits)
            active_radii = (rstep > r_crits)
            #sar = np.min(active_radii) if len(active_radii) > 0 else 1e99
            if len(active_radii) > 0:
                Nkn.append(np.sum(Nis[active_radii]))
                Nunact.append(np.sum(Nis[(rstep < r_crits)]))
            else:
                Nkn.append(0.0)
                Nunact.append(np.sum(Nis))

            print parcel.index[i], Neq[i], Nkn[i], Nunact[i], S_max, S

        Neq = np.array(Neq)
        Nkn = np.array(Nkn)
        Nunact = np.array(Nunact)

        parcel[species+'_Neq'] = Neq
        parcel[species+'_Nkn'] = Nkn
        parcel[species+'_Nunact'] = Nunact

        alphaz = Nkn/Neq
        alphaz[isnan(alphaz)] = 0.
        phiz = Nunact/Nkn
        phiz[phiz == inf] = 1.

        parcel[species+'_alpha'] = alphaz
        parcel[species+'_phi'] = phiz
    ##
    result_dict = {}
    result_dict['parcel'] = parcel
    result_dict['aerosols'] = aerosols

    results[r] = result_dict

figure(1)
s_maxes = [results[r]['parcel'].S.max() for r in aer2rs]
p = plot(aer2rs, s_maxes, label="%1.2f" % aer_meta.kappa)
for r, s in zip(aer2rs, s_maxes):
    act = np.any(results[r]['parcel']['carbon_Neq']) > 0.
    if act:
        plot(r, s, 'or')
for r in aer2rs:
    r_crit, s_crit = kohler_crit(T0, r*1e-6, aer_meta.kappa)
    print r, s_crit, aer_meta.kappa, T0
    plot(r, s_crit, marker="x", color=p[0].get_color())
ylim(0, 0.01)
xlabel("dry radius - carbon ($\mu$m)")
ylabel("Smax - 1.0")
legend()

