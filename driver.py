"""
.. module:: driver
    :synopsis: Routines for running/driving the parcel model

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""

from parcel import ParcelModel, ParcelModelError
from activation import fn2005, arg2000

def run_model(V, initial_aerosols, T, P, dt, max_steps=1000, t_end=500.,
              solver='lsoda', output='smax', solver_args={}):
    """
    Setup and run the parcel model with the given set of model and 
    integrator parameters.
    """
    if V <= 0:
        return 0.

    try: 
        model = ParcelModel(initial_aerosols, V, T, -0.0, P)
        Smax = model.run(t_end, dt, max_steps, solver=solver, 
                         output=output, solver_args=solver_args)
    except ParcelModelError:
        return None
    return Smax

def iterate_runs(V, initial_aerosols, T, P, dt=0.01, dt_iters=2, 
                 t_end=500., max_steps=500, output='smax'):
    """
    Iterate through several different strategies for integrating the parcel model.
    """
    aerosols = initial_aerosols
    if V <= 0:
        return 0., 0., 0.

    ## Check that there are actually aerosols to deal with
    aerosol_N = [a.distribution.N for a in initial_aerosols]
    if len(aerosol_N) == 1:
        if aerosol_N[0] < 0.01: return (-9999., -9999., -9999.)
    else:
        new_aerosols = []
        for i in xrange(len(aerosol_N)):
            if aerosol_N[i] > 0.01: 
                new_aerosols.append(initial_aerosols[i])
        aerosols = new_aerosols[:]

    S_max_arg, _ = arg2000(V, T, P, aerosols)
    S_max_fn, _ = fn2005(V, T, P, aerosols)

    dt_orig = dt*1.
    finished = False
    S_max = None

    ## Strategy 1: Try CVODE with modest tolerances.
    print " Trying CVODE with default tolerance"
    S_max = run_model(V, aerosols, T, P, dt, max_steps=2000, solver='cvode',
                      t_end=t_end, output=output, 
                      solver_args={'iter': 'Newton', 'time_limit': 10.0, 
                                   'linear_solver': "DENSE"})

    ## Strategy 2: Iterate over some increasingly relaxed tolerances for LSODA.
    if not S_max:
        while dt > dt_orig/(2**dt_iters):
            print " Trying LSODA, dt = %1.3e, max_steps = %d" % (dt, max_steps)
            S_max = run_model(V, aerosols, T, P, dt, max_steps, solver='lsoda',
                              t_end=t_end)
            if not S_max:
                dt = dt/2.
                print "    Retrying..."
            else:
                finished = True
                break

    ## Strategy 3: Last ditch numerical integration with LSODE. This will likely take a
    ##             a very long time.
    if not finished and not S_max:
        print " Trying LSODE"
        S_max = run_model(V, aerosols, T, P, dt_orig, max_steps=1000, solver='lsode',
                          t_end=t_end)

    ## Strategy 4: If all else fails return -9999.
    if not S_max:
        S_max = -9999.
        print " failed", V, dt

    return S_max, S_max_arg, S_max_fn
