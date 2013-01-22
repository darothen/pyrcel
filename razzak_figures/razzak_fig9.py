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
S0 = -0.005 # Supersaturation. 1-RH from wv term

z_top = 40.0 # meters
dt = 0.001 # seconds

Vs = np.logspace(np.log10(0.05), np.log10(5), 6)
#Vs = [1.0, ]

aerosol1 = AerosolSpecies('Nuclei Mode', Lognorm(mu=0.008, sigma=1.6, N=1000.), bins=200, kappa=0.6)
aerosol2 = AerosolSpecies('Accumulation Mode', Lognorm(mu=0.034, sigma=2.1, N=800.), bins=200, kappa=0.05)
aerosol3 = AerosolSpecies('Coarse Mode', Lognorm(mu=.46, sigma=2.2, N=0.72), bins=200, kappa=0.6)
initial_aerosols = [aerosol1, aerosol2, aerosol3]
aer_species = [a.species for a in initial_aerosols]
aer_dict = dict()
for aerosol in initial_aerosols:
    aer_dict[aerosol.species] = aerosol

fig, axes = subplots(3, 1, sharex=True, num=9, figsize=(8,11))

for V in Vs:
    zs = np.linspace(0, z_top, 20001)
    print zs[zs%1 == 0]
    
    ts = zs/V
    print "delta t =", np.diff(ts)[0]
    
    print V, "(%d)" % len(ts)

    print "   ... model run",
    pm = ParcelModel(initial_aerosols, V, T0, S0, P0, console=True)
    parcel, aerosols = pm.run(z_top, ts=ts)

    xs = np.arange(501)
    parcel = parcel.ix[parcel.index % 1. == 0]
    aero_subset = {}
    for key in aerosols:
        aerosol = aerosols[key]
        subset = aerosol.ix[aerosol.index % 1. == 0]
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
        ## Select correct axis
        if species == "Nuclei Mode": ax = axes[0]
        elif species == "Accumulation Mode": ax = axes[1]
        else: ax = axes[2]

        if V == Vs[0]:
            ax.plot(V, kn_frac, color="r", marker="D", linestyle="None", label="Kinetic")
            ax.plot(V, eq_frac, color="b", marker="D", linestyle="None", label="Equilibrium")
        else:
            ax.plot(V, kn_frac, color="r", marker="D", linestyle="None")
            ax.plot(V, eq_frac, color="b", marker="D", linestyle="None")
    print " done"

Vs = np.logspace(np.log10(0.01), np.log10(5.0), 100)
actfracs = []
for V in Vs:
    Smax, act = activation(V, T0, P0, initial_aerosols)
    actfracs.append(act)
actfracs = np.array(actfracs)
mode_acts = actfracs[:, 0], actfracs[:, 1], actfracs[:, 2]

legend_props = {'size':12}

ax1, ax2, ax3 = axes
for ax, act, aerosol in zip(axes, mode_acts, initial_aerosols):
    ax.plot(Vs, act, 'k-', label="Parameterized")
    ax.set_ylim(0, 1)
    ax.legend(loc="best", frameon=False, prop=legend_props)
    ax.set_ylabel("%s\n$N$umber Fraction Activated" % aerosol.species, multialignment="center")
    ax.semilogx()
ax.set_xlim(0.01, 10.0)
ax.set_xlabel("Updraft Velocity (m s$^{-1}$)")

draw()
savefig("fig9.pdf", transparent=True, bbox_inches="tight")
