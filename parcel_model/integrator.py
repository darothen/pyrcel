""" Interface to numerical ODE solvers.
"""
__docformat__ = 'reStructuredText'

available_integrators = ['odeint']

try:
    from odespy.odepack import Lsode, Lsoda
    from odespy import Vode
    available_integrators.extend(['lsode', 'lsoda', 'vode'])
except ImportError:
    print "Could not import odespy package; invoking the 'lsoda' or 'lsode' options will fail!"
    pass

try: 
    from assimulo.problem import Explicit_Problem
    from assimulo.solvers.sundials import CVode, CVodeError
    from assimulo.exception import TimeLimitExceeded
    available_integrators.extend(['cvode'])
except ImportError:
    print "Could not import Assimulo; invoking the CVode solver will fail!"
    pass

from scipy.integrate import odeint
from numpy import zeros

class Integrator(object):
    """
    Container class for the various integrators to use in the parcel model.

    All defined integrators should return a tuple whose first value ``x`` is either
    ``None`` or vector containing the parcel state at all requested timestamps, and
    whose second value is a boolean indicating whether the model run was successful.

    """
    @staticmethod
    def solver(method):
        """ Maps a solver name to a function.
        """
        solvers = {
            'odeint': Integrator._solve_odeint,
            'lsoda': Integrator._solve_lsoda,
            'lsode': Integrator._solve_lsode,
            'vode': Integrator._solve_vode,
            'cvode': Integrator._solve_cvode,
        }

        if method in available_integrators:
            return solvers[method]
        else:
            ## solver name is not defined, or the module containing
            ## it is unavailable
            raise ValueError("integrator for %s is not available" % method)

    @staticmethod
    def _solve_lsode(f, t, y0, args, console=False, max_steps=1000, terminate=False, 
                     **kwargs):
        """Wrapper for odespy.odepack.Lsode
        """
        nr = args[0]
        atol = [1e-4, 1e-4, 1e-4, 1e-10, 1e-10, 1e-8] + [1e-12]*nr
        rtol = 1e-7
        kwargs = { 'atol': atol, 'rtol': rtol, 'nsteps':max_steps }
        f_w_args = lambda u, t: f(u, t, *args)
        f_terminate = lambda u, t, step_no: u[step_no][5] < u[step_no -1][5]
        solver = Lsode(f_w_args, **kwargs)
        solver.set_initial_condition(y0)
        #solver.set(f_args=args)

        try:
            if terminate:
                x, t = solver.solve(t, f_terminate)
            else:
                x, t = solver.solve(t)
        except ValueError, e:
            raise ValueError("something broke in LSODE: %r" % e)
            return None, None, False

        return x, t, True

    @staticmethod
    def _solve_lsoda(f, t, y0, args, console=False, max_steps=1000, terminate=False, 
                     **kwargs):
        """Wrapper for odespy.odepack.Lsoda
        """
        nr = args[0]
        atol = [1e-4, 1e-4, 1e-4, 1e-10, 1e-10, 1e-8] + [1e-12]*nr
        rtol = 1e-7

        kwargs = { 'atol': atol, 'rtol': rtol, 'nsteps':max_steps }
        f_w_args = lambda u, t: f(u, t, *args)
        f_terminate = lambda u, t, step_no: u[step_no][5] < u[step_no-1][5]
        solver = Lsoda(f_w_args, **kwargs)
        solver.set_initial_condition(y0)
        #solver.set(f_args=args)

        try:
            if terminate:
                x, t = solver.solve(t, f_terminate)
            else:
                x, t = solver.solve(t)
        except ValueError, e:
            raise ValueError("something broke in LSODA: %r" % e)
            return None, None, False

        return x, t, True

    @staticmethod
    def _solve_vode(f, t, y0, args, console=False, max_steps=1000, terminate=False, 
                     **kwargs):
        """Wrapper for odespy.Vode
        """
        nr = args[0]
        atol = [1e-4, 1e-4, 1e-4, 1e-10, 1e-10, 1e-8] + [1e-12]*nr
        rtol = 1e-7

        #kwargs = { 'atol':1e-10, 'rtol':1e-8, 'nsteps':max_steps,
        #           'adams_or_bdf': 'bdf',  'order':5}
        kwargs = { 'nsteps': max_steps, 'adams_or_bdf': 'bdf', 'order': 5,
                   'atol': atol, 'rtol': rtol }
        f_w_args = lambda u, t: f(u, t, *args)
        f_terminate = lambda u, t, step_no: u[step_no][5] < u[step_no-1][5]
        solver = Vode(f_w_args, **kwargs)
        solver.set_initial_condition(y0)
        #solver.set(f_args=args)

        try:
            if terminate:
                x, t = solver.solve(t, f_terminate)
            else:
                x, t = solver.solve(t)
        except ValueError, e:
            raise ValueError("something broke in LSODA: %r" % e)
            return None, None, False

        return x, t, True

    @staticmethod
    def _solve_cvode(f, t, y0, args, console=False, max_steps=1000, terminate=False, 
                     **kwargs):
        """Wrapper for Assimulo's CVODE-Sundials routine
        """
        def rhs(t, u):
            return f(u, t, *args)

        ## Setup solver
        prob = Explicit_Problem(rhs, y0)
        sim = CVode(prob)
        sim.discr = 'BDF'
        sim.maxord = 5 
        sim.maxsteps = max_steps

        if "iter" in kwargs:
            sim.iter = kwargs['iter']
        else:
            sim.iter = 'Newton'

        if "time_limit" in kwargs:
            sim.time_limit = kwargs['time_limit']
            sim.report_continuously = True
        else:
            sim.time_limit = 0.0

        if "linear_solver" in kwargs:
            sim.linear_solver = kwargs['linear_solver']

        ## Setup tolerances
        nr = args[0]
        #ny = 6 + nr # z, P, T, wv, wc, S, *droplet_sizes        
        sim.rtol = 1e-7
        sim.atol = [1e-4, 1e-4, 1e-4, 1e-10, 1e-10, 1e-8] + [1e-12]*nr

        if not console:
            sim.verbosity = 50
        #sim.verbosity = 10
        #sim.report_continuously = True

        t_end = t[-1]
        steps = len(t)

        try:
            #print t_end, steps
            t, x = sim.simulate(t_end, steps)
        except CVodeError, e:
            raise ValueError("Something broke in CVode: %r" % e)
            return None, None,False
        except TimeLimitExceeded, e:
            raise ValueError("CVode took too long to complete")
            return None, None, False

        return x, t, True

    @staticmethod
    def _solve_odeint(f, t, y0, args, console=False, max_steps=1000, terminate=False, 
                     **kwargs):
        """Wrapper for scipy.integrate.odeint
        """
        nr = args[0]
        atol = [1e-4, 1e-4, 1e-4, 1e-10, 1e-10, 1e-8] + [1e-12]*nr
        rtol = 1e-7

        x, info = odeint(f, y0, t, args=args, full_output=1, mxhnil=0,
                         mxstep=max_steps, atol=atol, rtol=rtol)

        success = info['message'] == "Integration successful."

        if not success:
            print info
            raise ValueError("something broke in odeint: %r" % info['message'])
            return None, None, False

        return x, t, success
