'''Using MapReduce, explore the 'alpha' parameter space for the mixing experiment'''

from disco.job import Job
from disco.worker.classic.func import chain_reader, nop_map, nop_reduce
from disco.core import result_iterator
import os, sys

imports = ['numpy',
           ('parcel_model', "/home/darothen/workspace/parcel_model"),
           #'mixing_state',
]

class GenerateParameters(Job):
    """Prescribes the different values for alpha to generate different aerosol populations
    in order to run the model further down MapReduce line"""

    from disco_mixing import imports
    required_modules = imports

    reduce = staticmethod(nop_reduce)

    def __init__(self, partitions=1):
        super(GenerateParameters, self).__init__()
        self.partitions = partitions
        # Prescribe a single unused input so that only a single map instance is generated
        self.input = ["raw://0"]

    @staticmethod
    def map(_, params):
        """Generate the sets of parameters to use for creating the aerosol populations"""
        import numpy as np
        from itertools import product
        from random import shuffle

        if 'param_set' in params:
            parameter_sets = params['param_set']
        else:
            alphas = params['alphas']
            Vs = params['Vs']
            gammas = params['gammas']
            parameter_sets = [item for item in product(alphas, gammas, Vs)]
        shuffle(parameter_sets)

        ## discretize the parameter configurations and equitably distribute
        ## them for the next map instance to deal with.
        chunk_length = len(parameter_sets)/params['nprocs']
        leftover = len(parameter_sets) % params['nprocs']
        for n in xrange(params['nprocs']):
            if n < leftover:
                left = n*(1+chunk_length)
                to_yield = parameter_sets[left:left+1+chunk_length]
            else:
                left = leftover*(1+chunk_length) + (n-leftover)*chunk_length
                to_yield = parameter_sets[left:left+chunk_length]
            #print n, to_yield, len(to_yield)
            yield (n, to_yield)

class MakeAerosols(Job):
    """Given the different parameter sets, initialize the aerosol population to be run
    through the parcel model"""

    from disco_mixing import imports
    required_modules = imports

    ## Chain this phase to the earlier MapReduce phase
    map_reader = staticmethod(chain_reader)
    reduce = staticmethod(nop_reduce)

    def __init__(self, partitions=1):
        super(MakeAerosols, self).__init__()
        self.partitions = partitions

    @staticmethod
    def calc_mu(M, N, rho, sigma):
        """Compute mode radius (cm) from given M (kg/m^3), N (1/cm^3), rho (kg/m^3) and kappa"""
        import numpy as np
        mu_cubed = M*(3./(4.*np.pi))*(1./rho)*(1./N)*np.exp((-9./2.)*np.log(sigma)**2)
        return mu_cubed**(1./3.)

    @staticmethod
    def mixing2(M_tot, N_tot, alpha, gamma,  sigma_mixed, sigma_sulfate, sigma_carbon,
                rho_sulfate, rho_carbon, kappa_sulfate, kappa_carbon, diam_ratio=1.4):
        """
        Mass concentration in kg m**-3
        Number concentration in cm**-3

        diam_ratio = ratio of sulfur mode radius to carbon mode radius. Assert > 1!
        """
        import numpy as np
        from parcel_model.parcel import AerosolSpecies
        from parcel_model.lognorm import Lognorm

        # 1) Compute mixed mode density from prescribed values
        epsilon = 1./(gamma+1)
        rho_mixed = (1.-gamma)*rho_sulfate + gamma*rho_carbon
        kappa_mixed = (1.-gamma)*kappa_sulfate + gamma*kappa_carbon

        # 2) Separate internal/external masses
        M_ext = alpha*M_tot
        M_mixed = M_tot - M_ext

        # 3) Apportion between sulfate and carbon external modes
        M_sulfate = (epsilon/(1.+epsilon))*M_ext
        M_carbon = M_ext - M_sulfate

        # 4) Compute original (alpha = 0) mixed distribution parameters
        #mu_cubed = M_int*(3./(4.*np.pi))*(1./rho_mixed)*(1./N_int)*np.exp((-9./2.)*np.log(sigma_mixed)**2)
        #mu_mixed = mu_cubed**(1./3.) # cm
        mu_mixed = MakeAerosols.calc_mu(M_tot, N_tot, rho_mixed, sigma_mixed)

        # Compute N_mixed
        N_mixed = M_mixed/((4.*np.pi/3.)*rho_mixed*mu_mixed**3)*np.exp(-(9./2.)*np.log(sigma_mixed)**2)

        # 5) Compute number cocentration of external modes
        weighting_factor = (rho_carbon/rho_sulfate)*(diam_ratio**-3.)
        N_external = N_tot - N_mixed
        N_carbon = N_external/(1. + epsilon*weighting_factor)
        N_sulfate = N_external - N_carbon

        ## Finalize distributions
        # Mixed
        mixed = AerosolSpecies('mixed',
                               Lognorm(mu=mu_mixed*1e4, sigma=sigma_mixed, N=N_mixed),
                               kappa=kappa_mixed, bins=200)
        mixed.rho = rho_mixed

        ## Sulfate
        mu_sulfate = MakeAerosols.calc_mu(M_sulfate, N_sulfate, rho_sulfate, sigma_sulfate)
        sulfate = AerosolSpecies('sulfate',
                                 Lognorm(mu=mu_sulfate*1e4, sigma=sigma_sulfate, N=N_sulfate),
                                 kappa=kappa_sulfate, bins=200)
        sulfate.rho = rho_sulfate

        ## Carbon
        mu_carbon = MakeAerosols.calc_mu(M_carbon, N_carbon, rho_carbon, sigma_carbon)
        carbon = AerosolSpecies('carbon',
                                Lognorm(mu=mu_carbon*1e4, sigma=sigma_carbon, N=N_carbon),
                                kappa=kappa_carbon, bins=200)
        carbon.rho = rho_carbon

        return mixed, sulfate, carbon

    @staticmethod
    def mixing1(M_tot, N_tot, alpha, gamma, sigma_mixed, sigma_sulfate, sigma_carbon,
                rho_sulfate, rho_carbon, kappa_sulfate, kappa_carbon):
        """
        Mass concentration in kg m**-3
        Number concentration in cm**-3
        """
        import numpy as np
        from parcel_model.parcel import AerosolSpecies
        from parcel_model.lognorm import Lognorm

        # 1) Compute mixed mode density from prescribed values
        epsilon = 1./(gamma+1)
        rho_mixed = (1.-gamma)*rho_sulfate + gamma*rho_carbon
        kappa_mixed = (1.-gamma)*kappa_sulfate + gamma*kappa_carbon

        # 2) Separate internal/external masses
        M_ext = alpha*M_tot
        M_mixed = M_tot - M_ext

        # 3) Apportion between sulfate and carbon external modes
        M_sulfate = (epsilon/(1.+epsilon))*M_ext
        M_carbon = M_ext - M_sulfate

        # 4) Compute original (alpha = 0) mixed distribution parameters
        #mu_cubed = M_int*(3./(4.*np.pi))*(1./rho_mixed)*(1./N_int)*np.exp((-9./2.)*np.log(sigma_mixed)**2)
        #mu_mixed = mu_cubed**(1./3.) # cm
        mu_mixed = MakeAerosols.calc_mu(M_tot, N_tot, rho_mixed, sigma_mixed)

        # Compute N_mixed
        N_mixed = M_mixed/((4.*np.pi/3.)*rho_mixed*mu_mixed**3)*np.exp(-(9./2.)*np.log(sigma_mixed)**2)

        # 5) Compute number cocentration of external modes
        N_carbon = N_tot - N_mixed
        N_sulfate = N_carbon

        ## Finalize distributions
        # Mixed
        mixed = AerosolSpecies('mixed',
                               Lognorm(mu=mu_mixed*1e4, sigma=sigma_mixed, N=N_mixed),
                               kappa=kappa_mixed, bins=200)
        mixed.rho = rho_mixed

        ## Sulfate
        mu_sulfate = MakeAerosols.calc_mu(M_sulfate, N_sulfate, rho_sulfate, sigma_sulfate)
        sulfate = AerosolSpecies('sulfate',
                                 Lognorm(mu=mu_sulfate*1e4, sigma=sigma_sulfate, N=N_sulfate),
                                 kappa=kappa_sulfate, bins=200)
        sulfate.rho = rho_sulfate

        ## Carbon
        mu_carbon = MakeAerosols.calc_mu(M_carbon, N_carbon, rho_carbon, sigma_carbon)
        carbon = AerosolSpecies('carbon',
                                Lognorm(mu=mu_carbon*1e4, sigma=sigma_carbon, N=N_carbon),
                                kappa=kappa_carbon, bins=200)
        carbon.rho = rho_carbon

        return mixed, sulfate, carbon

    @staticmethod
    def map(parameter_sets, params):
        """Make the initial aerosol population given the parameter set"""
        from parcel_model.parcel import AerosolSpecies
        from parcel_model.lognorm import Lognorm
        import numpy as np
        import random

        M_tot = params['M_tot']
        N_tot = params['N_tot']
        flag = params['flag']
        rho_sulfate, rho_carbon = [params['rho'][s] for s in ['sulfate', 'carbon']]
        kappa_sulfate, kappa_carbon = [params['kappa'][s] for s in ['sulfate', 'carbon']]
        sigma_mixed = params['sigma']['mixed']
        sigma_sulfate = params['sigma']['sulfate']
        sigma_carbon = params['sigma']['carbon']


        ## Unload the parameters and iterate over them to drive the generation of aerosol
        ## populations and the model
        initial_aerosol_pops = []
        n, parameter_sets = parameter_sets
        for ps in parameter_sets:
            alpha, gamma, V = ps

            mixed, sulf, carbon = MakeAerosols.mixing2(M_tot, N_tot, alpha, gamma, sigma_mixed, sigma_sulfate, sigma_carbon, rho_sulfate, rho_carbon, kappa_sulfate, kappa_carbon)

            initial_aerosols = []
            if alpha > 0: # there is an external mode
                initial_aerosols.append(sulf)
                if flag == 'OC': ## if there is hygroscopic carbon
                    initial_aerosols.append(carbon)
            if alpha < 1: # there is leftover in the internal mode
                initial_aerosols.append(mixed)

            initial_aerosol_pops.append((initial_aerosols, ps))

        yield (n, initial_aerosol_pops)

class RunParcelModels(Job):
    """Reduce-only phase to run the parcel model given the aerosol populations
    previously generated"""

    from disco_mixing import imports
    required_modules = imports

    def __init__(self, partitions=1):
        super(RunParcelModels, self).__init__()
        self.partitions = partitions

    ## Chain this phase to the earlier MapReduce phase
    map_reader = staticmethod(chain_reader)
    reduce = staticmethod(nop_reduce)

    @staticmethod
    def simulation_pair(ps, initial_aerosols, V, T0, S0, P0,
                        z_top, dt, max_steps):
        from parcel_model.parcel import AerosolSpecies, ParcelModel, ParcelModelError
        from parcel_model.lognorm import Lognorm
        from parcel_model.micro import activate_ARG, act_fraction
        import numpy as np
        import random
        import os

        try:
            alpha, gamma, V = ps

            fact = z_top/30.
            zs = np.linspace(0, z_top, 100001*fact)
            ts = zs/V

            #print "I would've run with parameters alpha %0.2f and cbv %0.2f" % ps
            pm = ParcelModel(initial_aerosols, V, T0, S0, P0)
            parcel_out, aerosol_out = pm.run(z_top, dt=dt, max_steps=max_steps)
            #run_results = (parcel_out, aerosol_out)
            aerosols = pm.aerosols
        except ParcelModelError, e:
            #results.append((ps, None))
            return None

        activation_results = []

        ## Compute activation
        # 1) Parcel
        S_max = parcel_out['S'].max()
        S_max_idx = np.argmax(parcel_out.S)
        T_at_S_max = parcel_out['T'].ix[S_max_idx]

        N_act = 0.0
        N_total = 0.0
        all_results = {}
        for aerosol in initial_aerosols:
            species = aerosol.species
            #rs = aerosols_out[species].ix[S_max_idx].values
            rs = aerosol_out[species].ix[-1].values
            eq_frac, kn_frac = act_fraction(S_max, T_at_S_max, rs,
                                            aerosol.kappa, aerosol.r_drys, aerosol.Nis)
            N_act += eq_frac*aerosol.N
            N_total += aerosol.N

            all_results[species] = (eq_frac, aerosol.N)
        act_frac = N_act/N_total

        activation_results.append(('explicit', S_max, act_frac, N_act, all_results))

        # 2) ARG-Python
        S_max, fn = activate_ARG(V, T0, P0, initial_aerosols)
        N_total = 0.0
        N_act = 0.0
        all_results = {}
        for aerosol, f in zip(initial_aerosols, fn):
            N_act += f*aerosol.N
            N_total += aerosol.N

            all_results[species] = (f, aerosol.N)
        act_frac = N_act/N_total
        activation_results.append(('ARG-Python', S_max, act_frac, N_act, all_results))

        del pm

        return activation_results

    @staticmethod
    def map(initial_aerosol_pops, params):
        """Run the parcel model given given the aerosol population"""
        from parcel_model.parcel import ParcelModelError
        import time

        ## Pull model settings from params
        T0, S0, P0 = [params[s] for s in ('T0', 'S0', 'P0')]
        z_top, dt, max_steps = params['z_top'], params['dt'], params['max_steps']

        ## Helper method for re-submitting jobs which fail.
        def resubmit(ps, initial_aerosols, dt, max_steps):
            x = time.time()

            alpha, gamma, V = ps

            while dt >= 0.001:
                ## Try to run the model
                activation_results = RunParcelModels.simulation_pair(ps, initial_aerosols, V, T0, S0, P0, z_top, dt, max_steps)
                ## If it didn't work, report his and cut the timestep in half
                if not activation_results:
                    print "resubmitting %r with dt=%1.2e" % (ps, dt/2.,)
                    dt = dt/2.
                    max_steps = int(max_steps*3.)
                ## If it did work, we're done
                else:
                    break
            ## If we still don't have a good result after cutting dt several times,
            ## then report this.
            elapsed = time.time() - x
            if not activation_results:
                print "FAILED (%1.2e seconds) %r" % (elapsed, ps)
            else:
                print "SUCCESS (%1.2e seconds) %r" % (elapsed, ps)
            return activation_results

        results = []
        n, initial_aerosol_pops = initial_aerosol_pops
        n_runs = len(initial_aerosol_pops)
        for i, (initial_aerosols, ps) in enumerate(initial_aerosol_pops):
            print "EXECUTING RUN %d/%d" % (i+1, n_runs)
            component_results = {}

            ## INDIVIDUAL SPECIES
            param_fail = False
            for aerosol in initial_aerosols:
                species = aerosol.species

                #activation_results = RunParcelModels.simulation_pair(ps, [aerosol, ], V, T0, S0, P0, z_top, dt, max_steps)
                activation_results = resubmit(ps, [aerosol, ], dt, max_steps)

                if not activation_results:
                    results.append((ps, None))
                    param_fail = True
                    break
                else:
                    component_results[species] = activation_results

            if not param_fail:
            ## FULL MIXTURE
                #activation_results = RunParcelModels.simulation_pair(ps, initial_aerosols, V, T0, S0, P0, z_top, dt, max_steps)
                activation_results = resubmit(ps, initial_aerosols, dt, max_steps)

                if not activation_results:
                    results.append((ps, None))
                    continue

                component_results['mixture'] = activation_results
            results.append((ps, component_results))
        yield (n, results)

class ProcessRuns(Job):
    """Reduce-only phase to process the model results"""

    from disco_mixing import imports
    required_modules = imports

    def __init__(self, partitions=1):
        super(ProcessRuns, self).__init__()
        self.partitions = partitions

    ## Chain this phase to the earlier MapReduce phase
    map_reader = staticmethod(chain_reader)
    #map = staticmethod(nop_map)

    @staticmethod
    def map(model_results, params):
        n, model_results = model_results
        for result in model_results:
            yield (n, result)

    @staticmethod
    def reduce(model_results, out, params):
        """Compute the total number of activated CCN"""
        import numpy as np
        for n, (ps, activation_results) in model_results:

            if not activation_results:
                out.add(n, (ps, None))
            else:
                out.add(n, (ps, activation_results))

if __name__ == "__main__":

    from disco_mixing import GenerateParameters, MakeAerosols, RunParcelModels, ProcessRuns
    from parcel_model.parcel import AerosolSpecies
    from parcel_model.lognorm import Lognorm

    import numpy as np
    import time

    import pickle

    ## Original, totally mixed aerosol
    N_tot = 1000.0 # total number concentration, cm^-3
    M_tot = 400. # total mass concentration, ng m^-3
    M_tot = M_tot*1e-12 # convert to kg m^-3

    ## Prescribe gamma
    #gamma = 0.30
    #epsilon = 1./(gamma+1.)
    flag = "OC"
    if flag == "OC":
        kappa_carbon = 0.01
        rho_carbon = 2000. # kg/m^3
    if flag == "BC":
        kappa_carbon = 0.00
        rho_carbon = 1100.

    rho_sulfate = 1760. # kg m^-3
    kappa_sulfate = 0.54

    sigma_carbon = 2.0
    sigma_mixed = 2.0
    sigma_sulfate = 1.59

    log_sigma_carbon = np.log(sigma_carbon)
    log_sigma_mixed = np.log(sigma_mixed)
    log_sigma_sulfate = np.log(sigma_sulfate)

    rho_sulfate = 1760. # kg m^-3

    P0 = 80000. # Pressure, Pa
    T0 = 283.0 # Temperature, K
    S0 = -0.00 # Supersaturation. 1-RH from wv term

    z_top = 50.0 # meters
    dt = 0.016 # seconds
    max_steps = 100
    #V = 3. # meters/second

    nprocs = 60
    params = { 'kappa': {'sulfate': kappa_sulfate, 'carbon': kappa_carbon, },
               'rho':   {'sulfate': rho_sulfate,   'carbon': rho_carbon,   },
               'sigma': {'sulfate': sigma_sulfate, 'carbon': sigma_carbon, 'mixed': sigma_mixed},
               'M_tot': M_tot, 'N_tot': N_tot, 'flag': flag,
               'P0': P0, 'T0': T0, 'S0': S0, 'z_top': z_top, 'dt': dt,
               'nprocs': nprocs, 'max_steps': max_steps,
               'out_dir': "/net/legion/storage01/darothen/parcel/",
    }

    ##  Aerosol mixture computations
    # GENERATE PARAMETER SETS FROM SCRATCH

    alphas = np.linspace(0.0, 1.0, 11)
    gammas = np.linspace(0.1, 0.6, 6)
    Vs = np.logspace(np.log10(0.1), np.log10(1.5), 5)
    #Vs = [0.05 ,]
    params['alphas'] = alphas
    params['Vs'] = Vs
    params['gammas'] = gammas

    # SPECIFY PARAMETER SETS
    '''
    with open("run1.info", "rb") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith("failed_params"): break
        param_set_str = lines[i]
        param_set_str = param_set_str.split("params")[1].strip()
        param_set = eval(param_set_str)
    params['param_set'] = param_set
    raw_input("%d simulations - Continue?" % len(param_set))
    '''

    x = time.time()
    parameter_sets = GenerateParameters(partitions=nprocs).run(params=params)
    aerosols = MakeAerosols(partitions=nprocs).run(input=parameter_sets.wait(show=True), params=params)
    model_runs = RunParcelModels(partitions=nprocs).run(input=aerosols.wait(show=True), params=params)
    activation_results = ProcessRuns(partitions=nprocs).run(input=model_runs.wait(show=True), params=params)

    results = []
    for proc, result in result_iterator(activation_results.wait(show=True)):
        print "CPU ID - " + str(proc) + "\n" + "--"*20
        print result
        results.append(result)

    failed_params = []
    good_results = []
    for result in results:
        if not result[1]:
            failed_params.append(result[0])
        else:
            good_results.append(result)

    elapsed = time.time() - x
    failed_runs = len(failed_params)
    all_runs = len(results)
    print "Elapsed time - %r seconds" % elapsed
    print "   %d/%d runs did not finish" % (failed_runs, all_runs)

    out_f = open("run.info", "wb")
    out_f.write("params %r \n" % params)
    out_f.write("failed_params %r \n" % failed_params)
    out_f.write("elapsed %r \n" % elapsed)
    out_f.close()
    out_f = open("run.dat", "wb")
    pickle.dump(good_results, out_f)
    out_f.close()
