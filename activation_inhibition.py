from lognorm import Lognorm
from parcel import ParcelModel, AerosolSpecies
from micro import kohler_crit
import pandas

import numpy as np

P0 = 95000. # Pressure, Pa
T0 = 285.2 # Temperature, K
S0 = -0.02 # Supersaturation. 1-RH from wv term
V = 0.5 # m/s

z_top = 100.0 # meters
dt = 0.01 # seconds

## Iterate over number concentration of biomass burn aerosol
aer2Nis = np.logspace(2, 4, 15)
results = {}

for N in aer2Nis:
    print N
    aerosol1 = AerosolSpecies('sulf', Lognorm(mu=0.05, sigma=2.0, N=500.),
                              bins=50, kappa=0.6)
    aerosol2 = AerosolSpecies('carbon', {'r_drys': [0.200, ], 'Nis': [N, ]}, kappa=0.075)

    initial_aerosols = [aerosol1, aerosol2, ]
    aer_species = [a.species for a in initial_aerosols]
    aer_dict = dict()
    for aerosol in initial_aerosols:
        aer_dict[aerosol.species] = aerosol

    pm = ParcelModel(initial_aerosols, V, T0, S0, P0, console=False)
    parcel, aerosols = pm.run(z_top, dt)

    result_dict = {}
    result_dict['parcel'] = parcel
    result_dict['aerosols'] = aerosols

    xs = np.arange(501)
    parcel = parcel.ix[parcel.index % 1 == 0]
    aero_subset = {}
    for key in aerosols:
        aerosol = aerosols[key]
        subset = aerosol.ix[aerosol.index % 1 == 0]
        aero_subset[key] = subset
    aerosols2 = aero_subset

    for species in aer_species:
        print species
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
        alphaz[np.isnan(alphaz)] = 0.
        phiz = Nunact/Nkn
        phiz[phiz == np.inf] = 1.

        parcel[species+'_alpha'] = alphaz
        parcel[species+'_phi'] = phiz

        print "=="*35
        print species + " Summary - "
        print "Max activated fraction"
        print "   Eq: ", Neq.max()/np.sum(aer_meta.Nis)
        print "  Kin: ", Nkn.max()/np.sum(aer_meta.Nis)
        print ""
        print "Alpha maximum: %2.2f" % alphaz.max()
        print "  Phi maximum: %2.2f" % phiz.max()
        print "=="*35

        result_dict[species] = {'Eq_max': Neq.max()/np.sum(aer_meta.Nis),
                                'Kn_max': Nkn.max()/np.sum(aer_meta.Nis)}

    results[N] = result_dict