from lognorm import Lognorm
from parcel import ParcelModel, AerosolSpecies
from micro import kohler_crit, Rd, r_eff, activation

#from disco.core import result_iterator
#from disco.job import Job

import numpy as np

from collections import OrderedDict

import os, csv, operator

class ActivatedFractionJob(Job):
    required_modules = ['numpy',
                        ('lognorm', "/home/darothen/workspace/parcel_model/lognorm.py"),
                        ('parcel', "/home/darothen/workspace/parcel_model/parcel.py"),
                        ('parcel_aux', "/home/darothen/workspace/parcel_model/parcel_aux.so"),
                        ('micro', "/home/darothen/workspace/parcel_model/micro.py"),
                        'collections']
    partitions = 100

    @staticmethod
    def map_reader(fd, size, url, params):
        reader = csv.reader(fd, delimiter=',')
        for i, row in enumerate(reader):
            if i == 0: continue
            yield row

    def map(self, row, params):
        from razzak_disco import compute_act_fraction
        key, aerosol_rs = row[0], row[1:]
        aerosol_rs = np.array(aerosol_rs, dtype=np.float)
        S, T = params[key]['S'], params[key]['T']
        aer_meta = params['aer_meta']
        #yield key, compute_act_fraction(S, T, aer_meta, aerosol_rs)
        yield key, (S, T, aer_meta, aerosol_rs)

    def reduce(self, rows_iter, out, params):
        for key, args in rows_iter:
            out.add(float(key), compute_act_fraction(*args))

def compute_act_fraction(S, T, aer_meta, aerosol_rs):
    kappa = aer_meta.kappa
    r_drys = aer_meta.r_drys
    Nis = aer_meta.Nis

    r_crits, s_crits = zip(*[kohler_crit(T, r_dry, kappa) for r_dry in  r_drys])

    s_crits = numpy.array(s_crits)
    r_crits = numpy.array(r_crits)

    #big_s =  S_max >= s_crits
    big_s =  S >= s_crits

    rstep = aerosol_rs
    active_radii = (rstep > r_crits)

    N_tot = numpy.sum(Nis)

    Neq = numpy.sum(Nis[big_s])/N_tot
    Nkn = numpy.sum(Nis[active_radii])/N_tot
    Nunact = numpy.sum(Nis[(rstep < r_crits)])/N_tot

    return Neq, Nkn, Nunact

if __name__ == "__main__":

    P0 = 100000. # Pressure, Pa
    T0 = 294.0 # Temperature, K
    S0 = -0.00 # Supersaturation. 1-RH from wv term

    z_top = 30.0 # meters
    dt = 0.01 # seconds

    Vs = [0.5]

    aerosol1 = AerosolSpecies('Mode 1', Lognorm(mu=0.2, sigma=2.0, N=100.),
                              bins=200, kappa=0.6)
    aerosol2 = AerosolSpecies('Mode 2', Lognorm(mu=0.02, sigma=2.0, N=100.),
                              bins=200, kappa=0.05)
    initial_aerosols = [aerosol1, aerosol2, ]
    aer_species = [a.species for a in initial_aerosols]
    aer_dict = dict()
    for aerosol in initial_aerosols:
        aer_dict[aerosol.species] = aerosol

    for V in Vs:
        zs = np.linspace(0, z_top, 6001)
        print zs[zs%1 == 0]

        ts = zs/V
        print "delta t =", np.diff(ts)[0]

        print V, "(%d)" % len(ts)

        print "   ... model run",

        pm = ParcelModel(initial_aerosols, V, T0, S0, P0, console=False)
        parcel, aerosols = pm.run(z_top, ts=ts)

        parcel = parcel.ix[parcel.index % 1. == 0]
        aero_subset = {}
        for key in aerosols:
            aerosol = aerosols[key]
            subset = aerosol.ix[aerosol.index % 1. == 0]
            aero_subset[key] = subset
        aerosols = aero_subset
        pm.write_csv(parcel, aerosols, "./temp_data")

        print "done"

        print "   ... activation"
        params = OrderedDict()
        with open("./temp_data/parcel.csv", "r") as param_file:
            for i, line in enumerate(param_file.readlines()):
                if i == 0: continue
                key, P, S, T, wc, wv = line.strip().split(',')
                params[key] = {'T': float(T), 'S': float(S)}

        ## Compute Activation stats
        from razzak_disco import ActivatedFractionJob
        for species in aer_species:
            if species == "Mode 2": continue
            aer_meta = aer_dict[species]

            params['aer_meta'] = aer_meta

            Neq = []
            Nkn = []
            Nunact = []
            S_max = S0

            mode_fn = "/home/darothen/workspace/parcel_model/temp_data/%s.csv" % species
            print mode_fn
            disco_job = ActivatedFractionJob().run(input=[mode_fn],
                                                   params=params)
            for key, result in sorted(result_iterator(disco_job.wait(show=True))):
                print key, result
                eq, kn, unact = result
                Neq.append(eq)
                Nkn.append(kn)
                Nunact.append(unact)

            ###

            Neq = np.array(Neq)
            Nkn = np.array(Nkn)
            Nunact = np.array(Nunact)

            print Neq
            print Nkn
            print Nunact

        print " done"