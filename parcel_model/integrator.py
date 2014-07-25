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
    from assimulo.solvers.odepack import LSODAR
    from assimulo.exception import TimeLimitExceeded
    available_integrators.extend(['cvode', 'lsodar'])
except ImportError:
    print "Could not import Assimulo; invoking the CVode solver will fail!"
    pass

from functools import partial
from scipy.integrate import odeint
import numpy as np

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
            #'cvode': Integrator._solve_cvode,
            #'lsodar': Integrator._solve_lsodar,
            'cvode': partial(Integrator._solve_with_assimulo, method='cvode'),
            'lsodar': partial(Integrator._solve_with_assimulo, method='lsodar'),
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
    def _solve_with_assimulo(f, t, y0, args, console=False, max_steps=1000, terminate=False, 
                             method='cvode', **kwargs):
        """ Wrapper for Assimulo's solver routines
        """
        def user_rhs(t, y):
            dode_dt = f(y, t, *args)
            return dode_dt

        class Extended_Problem(Explicit_Problem):

            name  = 'Parcel model ODEs'
            sw0   = [True,  # Normal integration switch
                     False, # Past cut-off switch
            ] 
            t_cutoff = 1e5
            dS_dt = 1.0

            def rhs(self, t, y, sw):
                if not sw[1]: # Normal integration before cutoff
                    dode_dt = f(y, t, *args)
                    self.dS_dt = dode_dt[5]
                else:
                    dode_dt = np.zeros(6 + args[0])
                return dode_dt

            # The event function
            def state_events(self, t, y, sw):
                """ Check whether an 'event' has occurred. We want to see if the
                supersaturation is decreasing or not. """
                if sw[0]: 
                    #dode_dt = f(y, t, *args)
                    #smax_event = dode_dt[5]
                    smax_event = self.dS_dt
                else:
                    smax_event = -1.0

                t_cutoff_event = t - self.t_cutoff

                return np.array([smax_event, t_cutoff_event])

            # Event handling function
            def handle_event(self, solver, event_info):
                """ Event handling. This function is called when Assimulo finds
                an event as specified by the event function. """
                event_info = event_info[0] # Only state events, event_info[1] is time events
                if event_info[0] != 0:
                    solver.sw[0] = False
                    self.t_cutoff = solver.t + 5.0

            def handle_result(self, solver, t, y):
                if t < self.t_cutoff:
                    Explicit_Problem.handle_result(self, solver, t, y)

        ## Setup solver
        if terminate:
            prob = Extended_Problem(y0=y0)
        else:
            prob = Explicit_Problem(user_rhs, y0)

        ## Choose simulator
        if method == "cvode":
            sim = CVode(prob)
            sim.discr = 'BDF'
            sim.maxord = 5 

            if "iter" in kwargs:
                sim.iter = kwargs['iter']
            else:
                sim.iter = 'Newton'

            if "linear_solver" in kwargs:
                sim.linear_solver = kwargs['linear_solver']

        elif method == "lsodar":
            sim = LSODAR(prob)
            sim.maxords = 5

        else:
            raise ValueError("Passed method (%r) must be 'cvode' or 'lsodar'" % method)

        sim.maxsteps = max_steps

        if "time_limit" in kwargs:
            sim.time_limit = kwargs['time_limit']
            sim.report_continuously = True
        else:
            sim.time_limit = 0.0

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
            return None, None, False
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
