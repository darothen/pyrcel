# -*- coding: utf-8 -*-
""" Interface to numerical ODE solvers.
"""
import sys

# Compatibility - timer functions
# In Python 3, the more accurate `time.process_time()` method is available. But
# for legacy support, can default instead to `time.clock()`
import time
import warnings
from abc import ABCMeta, abstractmethod

import numpy as np

if sys.version_info[0] < 3:
    timer = time.clock
else:
    timer = time.process_time

from . import constants as c

available_integrators = ["odeint"]

try:
    from odespy import Vode
    from odespy.odepack import Lsoda, Lsode

    available_integrators.extend(["lsode", "lsoda", "vode"])
except ImportError:
    warnings.warn(
        "Could not import odespy package; "
        "invoking the 'lsoda' or 'lsode' options will fail!"
    )


try:
    # from assimulo.solvers.odepack import LSODAR
    from assimulo.exception import TimeLimitExceeded
    from assimulo.problem import Explicit_Problem
    from assimulo.solvers.sundials import CVode, CVodeError

    available_integrators.extend(["cvode", "lsodar"])
except ImportError:
    warnings.warn("Could not import Assimulo; " "invoking the CVode solver will fail!")


__all__ = ["Integrator"]

state_atol = [1e-4, 1e-4, 1e-4, 1e-10, 1e-10, 1e-4, 1e-8]
state_rtol = 1e-7


class Integrator(metaclass=ABCMeta):
    """
    Container class for the various integrators to use in the parcel model.

    All defined integrators should return a tuple whose first value ``x`` is either
    ``None`` or vector containing the parcel state at all requested timestamps, and
    whose second value is a boolean indicating whether the model run was successful.

    """

    def __init__(self, rhs, output_dt, solver_dt, y0, args, t0=0.0, console=False):
        self.output_dt = output_dt
        self.solver_dt = solver_dt
        self.y0 = y0
        self.t0 = t0
        self.console = console

        self.args = args

        def _user_rhs(t, y):
            dode_dt = rhs(y, t, *self.args)
            return dode_dt

        self.rhs = _user_rhs

    @abstractmethod
    def integrate(self, t_end, **kwargs):
        pass

    @abstractmethod
    def __repr__(self):
        pass

    @staticmethod
    def solver(method):
        """Maps a solver name to a function."""
        solvers = {
            # # SciPy interfaces
            # 'odeint': Integrator._solve_odeint,
            # # ODESPY interfaces
            # 'lsoda': partial(Integrator._solve_with_odespy, method='lsoda'),
            # 'lsode': partial(Integrator._solve_with_odespy, method='lsode'),
            # 'vode': partial(Integrator._solve_with_odespy, method='vode'),
            # # Assimulo interfaces
            # 'cvode': partial(Integrator._solve_with_assimulo, method='cvode'),
            # 'lsodar': partial(Integrator._solve_with_assimulo, method='lsodar'),
            "cvode": CVODEIntegrator
        }

        if method in available_integrators:
            return solvers[method]
        else:
            # solver name is not defined, or the module containing
            # it is unavailable
            raise ValueError("integrator for %s is not available" % method)


class ExtendedProblem(Explicit_Problem):
    """This extension of the Assimulo 'Explicit_Problem' class
    encodes some of the logic particular to the parcel model simulation,
    specifically rules for terminating the simulation and detecting
    events such as the maximum supersaturation occurring"""

    name = "Parcel model ODEs"
    sw0 = [True, False]  # Normal integration switch  # Past cut-off switch
    t_cutoff = 1e5
    dS_dt = 1.0

    def __init__(self, rhs_fcn, rhs_args, terminate_depth, *args, **kwargs):
        self.rhs_fcn = rhs_fcn
        self.rhs_args = rhs_args
        self.V = rhs_args[3]
        self.terminate_time = terminate_depth / self.V
        super(Explicit_Problem, self).__init__(*args, **kwargs)

    def rhs(self, t, y, sw):
        if not sw[1]:  # Normal integration before cutoff
            dode_dt = self.rhs_fcn(t, y)  # FROM THE CVODEINTEGRATOR
            self.dS_dt = dode_dt[c.N_STATE_VARS - 1]
        else:
            # There may be a bug here. I can't recall when this branch is ever run; it
            # seems to zero out the state derivative, but to construct that array it
            # should be looking at self.rhs_args, not self.args (which isn't saved).
            # I'm going to comment out this line which I think is broken and replace it
            # with the correct one for now, but leave a record of this change
            # Daniel Rothenberg <daniel@danielrothenberg.com> - 2/15/2016
            # dode_dt = np.zeros(c.N_STATE_VARS + self.args[0])  # FROM INIT ARGS
            dode_dt = np.zeros(c.N_STATE_VARS + self.rhs_args[0])
        return dode_dt

    # The event function
    def state_events(self, t, y, sw):
        """Check whether an 'event' has occurred. We want to see if the
        supersaturation is decreasing or not."""
        if sw[0]:
            smax_event = self.dS_dt
        else:
            smax_event = -1.0

        t_cutoff_event = t - self.t_cutoff

        return np.array([smax_event > 0, t_cutoff_event < 0])

    # Event handling function
    def handle_event(self, solver, event_info):
        """Event handling. This function is called when Assimulo finds
        an event as specified by the event function."""
        event_info = event_info[0]  # Only state events, event_info[1] is time events
        if event_info[0] != 0:
            solver.sw[0] = False
            self.t_cutoff = solver.t + self.terminate_time

    def handle_result(self, solver, t, y):
        if t < self.t_cutoff:
            Explicit_Problem.handle_result(self, solver, t, y)


class CVODEIntegrator(Integrator):
    kwargs = None  # Save the kwargs used for setting up the interface to CVODE!

    def __init__(
        self,
        rhs,
        output_dt,
        solver_dt,
        y0,
        args,
        t0=0.0,
        console=False,
        terminate=False,
        terminate_depth=100.0,
        **kwargs
    ):
        self.terminate = terminate
        super(CVODEIntegrator, self).__init__(
            rhs, output_dt, solver_dt, y0, args, t0, console
        )

        # Setup solver
        if terminate:
            self.prob = ExtendedProblem(
                self.rhs, self.args, terminate_depth, y0=self.y0
            )
        else:
            self.prob = Explicit_Problem(self.rhs, self.y0)

        self.sim = self._setup_sim(**kwargs)

        self.kwargs = kwargs

    def _setup_sim(self, **kwargs):
        """Create a simulation interface to Assimulo using CVODE, given
        a problem definition

        """

        # Create Assimulo interface
        sim = CVode(self.prob)
        sim.discr = "BDF"
        sim.maxord = 5

        # Setup some default arguments for the ODE solver, or override
        # if available. This is very hackish, but it's fine for now while
        # the number of anticipated tuning knobs is small.
        if "maxh" in kwargs:
            sim.maxh = kwargs["maxh"]
        else:
            sim.maxh = np.min([0.1, self.output_dt])

        if "minh" in kwargs:
            sim.minh = kwargs["minh"]
        # else: sim.minh = 0.001

        if "iter" in kwargs:
            sim.iter = kwargs["iter"]
        else:
            sim.iter = "Newton"

        if "linear_solver" in kwargs:
            sim.linear_solver = kwargs["linear_solver"]

        if "max_steps" in kwargs:  # DIFFERENT NAME!!!!
            sim.maxsteps = kwargs["max_steps"]
        else:
            sim.maxsteps = 1000

        if "time_limit" in kwargs:
            sim.time_limit = kwargs["time_limit"]
            sim.report_continuously = True
        else:
            sim.time_limit = 0.0

        # Don't save the [t_-, t_+] around events
        sim.store_event_points = False

        # Setup tolerances
        nr = self.args[0]
        sim.rtol = state_rtol
        sim.atol = state_atol + [1e-12] * nr

        if not self.console:
            sim.verbosity = 50
        else:
            sim.verbosity = 40
        # sim.report_continuously = False

        # Save the Assimulo interface
        return sim

    def integrate(self, t_end, **kwargs):
        # Compute integration logic. We need to know:
        # 1) How are we iterating the solver loop?
        t_increment = self.solver_dt
        # 2) How many points do we want to interpolate for output?
        n_out = int(self.solver_dt / self.output_dt)
        t_current = self.t0

        if self.console:
            print()
            print("Integration Loop")
            print()
            print("  step     time  walltime  Δwalltime |     z       T       S")
            print(" " "------------------------------------|----------------------")
            step_fmt = (
                " {:5d} {:7.2f}s  {:7.2f}s  {:8.2f}s |" " {:5.1f} {:7.2f} {:6.2f}%"
            )

        txs, xxs = [], []
        n_steps = 1
        total_walltime = 0.0
        now = timer()
        while t_current < t_end:
            if self.console:
                # Update timing estimates
                delta_walltime = timer() - now
                total_walltime += delta_walltime

                # Grab state vars
                state = self.y0 if n_steps == 1 else xxs[-1][-1]
                _z = state[c.STATE_VAR_MAP["z"]]
                _T = state[c.STATE_VAR_MAP["T"]]
                _S = state[c.STATE_VAR_MAP["S"]] * 100
                print(
                    step_fmt.format(
                        n_steps,
                        t_current,
                        total_walltime,
                        delta_walltime,
                        _z,
                        _T,
                        _S,
                    )
                )
            try:
                now = timer()
                out_list = np.linspace(t_current, t_current + t_increment, n_out + 1)
                tx, xx = self.sim.simulate(t_current + t_increment, 0, out_list)
            except CVodeError as e:
                raise ValueError("Something broke in CVode: %r" % e)
            except TimeLimitExceeded:
                raise ValueError("CVode took too long to complete")

            if n_out == 1:
                txs.append(tx[-1])
                xxs.append(xx[-1])
            else:
                txs.extend(tx[:-1])
                xxs.append(xx[:-1])
            t_current = tx[-1]

            # Has the max been found and can we terminate?
            if self.terminate:
                if not self.sim.sw[0]:
                    if self.console:
                        print("---- termination condition reached ----")
                    break

            n_steps += 1
        if self.console:
            print("---- end of integration loop ----")

        # Determine output information
        t = np.array(txs)
        if n_out == 1:  # can just merge the outputs
            x = np.array(xxs)
        else:  # Need to concatenate lists of outputs
            x = np.concatenate(xxs)

        return x, t, True

    def __repr__(self):
        return "CVODE integrator - direct Assimulo interface"
