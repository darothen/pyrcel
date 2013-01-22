from lognorm import Lognorm
from parcel import ParcelModel, AerosolSpecies
from micro import kohler_crit, Rd, r_eff, activation

import numpy as np
from pylab import *
ion()

rc('text', usetex=True)
rc('font', family='serif')
rc('font', size=16)
rc('legend', fontsize=12)

P0 = 100000. # Pressure, Pa
T0 = 294.0 # Temperature, K
S0 = 0.00 # Supersaturation. 1-RH from wv term
V = 0.5 # m/s

z_top = 20.0 # meters
dt = 0.001 # seconds

mode2_Ns = np.logspace(np.log10(100), np.log10(5000), 6)
#mode2_Ns = [500., ]

aerosol1 = AerosolSpecies('Mode 1', Lognorm(mu=0.05, sigma=2.0, N=100.),
                          bins=200, kappa=0.6)

fig, axes = subplots(2, 1, sharex=True, num=2)

for N in mode2_Ns:
    print N

    aerosol2 = AerosolSpecies('Mode 2', Lognorm(mu=0.05, sigma=2.0, N=N),
                              bins=200, kappa=0.05)

    initial_aerosols = [aerosol1, aerosol2, ]
    aer_species = [a.species for a in initial_aerosols]
    aer_dict = dict()
    for aerosol in initial_aerosols:
        aer_dict[aerosol.species] = aerosol

    print "   ... model run",
    pm = ParcelModel(initial_aerosols, V, T0, S0, P0, console=False)
    parcel, aerosols = pm.run(z_top, dt)

    xs = np.arange(501)
    parcel = parcel.ix[parcel.index % 1 == 0]
    aero_subset = {}
    for key in aerosols:
        aerosol = aerosols[key]
        subset = aerosol.ix[aerosol.index % 1 == 0]
        aero_subset[key] = subset
    aerosols2 = aero_subset
    print " done"

    print "   ... activation",
    ## Compute Activation stats
    for species in aer_species:
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

            print species, parcel.index[i], Neq[i], Nkn[i], Nunact[i], S_max, S

        Neq = np.array(Neq)
        Nkn = np.array(Nkn)
        Nunact = np.array(Nunact)

        kn_frac = Nkn.max()/np.sum(aer_meta.Nis)
        eq_frac = Neq.max()/np.sum(aer_meta.Nis)

        ax = axes[0] if species == "Mode 1" else axes[1]
        lbl_kn, lbl_eq = ("Kinetic", "Equilibrium") if N == mode2_Ns[0] else (None, None)
        ax.plot(N, kn_frac, color="r", marker="D", linestyle="None", label=lbl_kn)
        ax.plot(N, eq_frac, color="b", marker="D", linestyle="None", label=lbl_eq)
    print " done"

Ns = np.logspace(np.log10(100), np.log10(5000), 100)
actfracs = []
for N in Ns:
    aerosol2 = AerosolSpecies('Mode 2', Lognorm(mu=0.05, sigma=2.0, N=N), bins=200, kappa=0.05)
    initial_aerosols = [aerosol1, aerosol2, ]
    Smax, act = activation(V, T0, P0, initial_aerosols)
    actfracs.append(act)
actfracs = np.array(actfracs)
mode1_act, mode2_act = actfracs[:, 0], actfracs[:, 1]

legend_props = {'size':12}

ax1, ax2 = axes
ax1.plot(Ns, mode1_act, 'k-', label="Parameterized")
ax1.set_ylim(0, 1)
ax1.legend(loc="upper right", frameon=False, prop=legend_props)
ax1.set_ylabel("Mode 1\n$N$ Fraction Activated", multialignment="center")

ax2.plot(Ns, mode2_act, 'k-', label="Parameterized")
ax2.set_ylim(0, 1)
ax2.legend(loc="upper right", frameon=False, prop=legend_props)
ax2.set_xlabel("Mode 2 Aerosol Number (cm$^{-3}$)")
ax2.set_ylabel("Mode 2\n$N$ Fraction Activated", multialignment="center")

draw()
