""" Utilities for driving sets of parcel model integration strategies.

Occasionally, a pathological set of input parameters to the parcel model
will really muck up the ODE solver's ability to integrate the model.
In that case, it would be nice to quietly adjust some of the numerical
parameters for the ODE solver and re-submit the job. This module includes a
workhorse function :func:`iterate_runs` which can serve this purpose and can
serve as an example for more complex integration strategies. Alternatively,
:func:`run_model`is a useful shortcut for building/running a model and snagging
its output.

"""
from numpy import empty, nan
from pandas import DataFrame

from pyrcel.util import ParcelModelError

from .activation import arg2000, mbn2014
from .parcel import ParcelModel


def run_model(
    V,
    initial_aerosols,
    T,
    P,
    dt,
    S0=-0.0,
    max_steps=1000,
    t_end=500.0,
    solver="lsoda",
    output_fmt="smax",
    terminate=False,
    solver_kws=None,
    model_kws=None,
):
    """Setup and run the parcel model with given solver configuration.

    Parameters
    ----------
    V, T, P : float
        Updraft speed and parcel initial temperature and pressure.
    S0 : float, optional, default 0.0
        Initial supersaturation, as a percent. Defaults to 100% relative humidity.
    initial_aerosols : array_like of :class:`AerosolSpecies`
        Set of aerosol populations contained in the parcel.
    dt : float
        Solver timestep, in seconds.
    max_steps : int, optional, default 1000
        Maximum number of steps per solver iteration. Defaults to 1000; setting
        excessively high could produce extremely long computation times.
    t_end : float, optional, default 500.0
        Model time in seconds after which the integration will stop.
    solver : string, optional, default 'lsoda'
        Alias of which solver to use; see :class:`Integrator` for all options.
    output_fmt : string, optional, default 'smax'
        Alias indicating which output format to use; see :class:`ParcelModel` for
        all options.
    solver_kws, model_kws : dicts, optional
        Additional arguments/configuration to pass to the numerical integrator or model.

    Returns
    -------
    Smax : (user-defined)
        Output from parcel model simulation based on user-specified `output_fmt` argument. See
        :class:`ParcelModel` for details.

    Raises
    ------
    ParcelModelError
        If the model fails to initialize or breaks during runtime.

    """
    # Setup kw dicts
    if model_kws is None:
        model_kws = {}
    if solver_kws is None:
        solver_kws = {}

    if V <= 0:
        return 0.0

    try:
        model = ParcelModel(initial_aerosols, V, T, S0, P, **model_kws)
        Smax = model.run(
            t_end,
            dt,
            max_steps,
            solver=solver,
            output_fmt=output_fmt,
            terminate=terminate,
            **solver_kws
        )
    except ParcelModelError:
        return None
    return Smax


def iterate_runs(
    V,
    initial_aerosols,
    T,
    P,
    S0=-0.0,
    dt=0.01,
    dt_iters=2,
    t_end=500.0,
    max_steps=500,
    output_fmt="smax",
    fail_easy=True,
):
    """Iterate through several different strategies for integrating the parcel model.

    As long as `fail_easy` is set to `False`, the strategies this method implements are:

    1. **CVODE** with a 10 second time limit and 2000 step limit.
    2. **LSODA** with up to `dt_iters` iterations, where the timestep `dt` is
       halved each time.
    3. **LSODE** with coarse tolerance and the original timestep.

    If these strategies all fail, the model will print a statement indicating such
    and return either -9999 if `output_fmt` was 'smax', or an empty array or DataFrame
    accordingly.

    Parameters
    ----------
    V, T, P : float
        Updraft speed and parcel initial temperature and pressure.
    S0 : float, optional, default 0.0
        Initial supersaturation, as a percent. Defaults to 100% relative humidity.
    initial_aerosols : array_like of :class:`AerosolSpecies`
        Set of aerosol populations contained in the parcel.
    dt : float
        Solver timestep, in seconds.
    dt_iters : int, optional, default 2
        Number of times to halve `dt` when attempting **LSODA** solver.
    max_steps : int, optional, default 1000
        Maximum number of steps per solver iteration. Defaults to 1000; setting
        excessively high could produce extremely long computation times.
    t_end : float, optional, default 500.0
        Model time in seconds after which the integration will stop.
    output : string, optional, default 'smax'
        Alias indicating which output format to use; see :class:`ParcelModel` for
        all options.
    fail_easy : boolean, optional, default `True`
        If `True`, then stop after the first strategy (**CVODE**)

    Returns
    -------
    Smax : (user-defined)
        Output from parcel model simulation based on user-specified `output` argument. See
        :class:`ParcelModel` for details.

    """
    aerosols = initial_aerosols
    if V <= 0:
        return 0.0, 0.0, 0.0

    # Check that there are actually aerosols to deal with
    aerosol_N = [a.distribution.N for a in initial_aerosols]
    if len(aerosol_N) == 1:
        if aerosol_N[0] < 0.01:
            return -9999.0, -9999.0, -9999.0
    else:
        new_aerosols = []
        for i in range(len(aerosol_N)):
            if aerosol_N[i] > 0.01:
                new_aerosols.append(initial_aerosols[i])
        aerosols = new_aerosols[:]

    S_max_arg, _, _ = arg2000(V, T, P, aerosols)
    S_max_fn, _, _ = mbn2014(V, T, P, aerosols)

    dt_orig = dt * 1.0
    finished = False
    S_max = None

    # Strategy 1: Try CVODE with modest tolerances.
    print(" Trying CVODE with default tolerance")
    S_max = run_model(
        V,
        aerosols,
        T,
        P,
        dt,
        S0=S0,
        max_steps=2000,
        solver="cvode",
        t_end=t_end,
        output_fmt=output_fmt,
        solver_kws={
            "iter": "Newton",
            "time_limit": 10.0,
            "linear_solver": "DENSE",
        },
    )

    # Strategy 2: Iterate over some increasingly relaxed tolerances for LSODA.
    if (S_max is None) and not fail_easy:
        while dt > dt_orig / (2**dt_iters):
            print(" Trying LSODA, dt = %1.3e, max_steps = %d" % (dt, max_steps))
            S_max = run_model(
                V,
                aerosols,
                T,
                P,
                dt,
                S0,
                max_steps,
                solver="lsoda",
                t_end=t_end,
            )
            if not S_max:
                dt /= 2.0
                print("    Retrying...")
            else:
                finished = True
                break

    # Strategy 3: Last ditch numerical integration with LSODE. This will likely take a
    #             a very long time.
    if (not finished) and (S_max is None) and (not fail_easy):
        print(" Trying LSODE")
        S_max = run_model(
            V,
            aerosols,
            T,
            P,
            dt_orig,
            max_steps=1000,
            solver="lsode",
            t_end=t_end,
            S0=S0,
        )

    # Strategy 4: If all else fails return -9999.
    if S_max is None:
        if output_fmt == "smax":
            S_max = -9999.0
        elif output_fmt == "arrays":
            S_max = empty([0]), empty([0])
        elif output_fmt == "dataframes":
            S_max = (
                DataFrame(data={"S": [nan]}),
                DataFrame(data={"aerosol1": [nan]}),
            )
        else:
            S_max = nan
        print(" failed", V, dt)

    return S_max, S_max_arg, S_max_fn
    return S_max, S_max_arg, S_max_fn
    return S_max, S_max_arg, S_max_fn
