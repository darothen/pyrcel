from parcel_model.parcel import AerosolSpecies
from parcel_model.lognorm import Lognorm

import numpy as np

n_type = 1
n_bins = 40
n_species = 1
n_modes = 1

densities = [1.7418e3, ] # kg/m^3
molws = [0.132, ] # kg/mol
ionns = [3., ]
Ns = [200., ] # number concentrations 1/cm^3
dms = [0.02, ] # mean diameter, microns
sigs = [2.5, ]


if __name__ == "__main__":
## SINGLE AEROSOL

    for rho, molw, ion, N, dm, sig in zip(densities, molws, ionns, Ns, dms, sigs):
        aerosol = AerosolSpecies('for_ming',
                                 Lognorm(mu=dm/2., N=N, sigma=sig), kappa=0.1, bins=n_bins)
        print aerosol
    dry_dps = 2.e6*aerosol.r_drys
    masses = (np.pi/6.)*rho*((dry_dps*1e-6)**3)*aerosol.Nis
    masses = masses*1e3

    with open("Mass.dat", "wb") as f:
        for d, m in zip(dry_dps, masses):
            f.write("%r %r \n" % (d, m))
