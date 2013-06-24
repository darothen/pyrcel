"""
.. module:: parcel
    :synopsis: Numerical integrators for solving the parcel model.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""
__docformat__ = 'reStructuredText'

try:
    from odespy.odepack import Lsode, Lsoda
except ImportError:
    pass
from scipy.integrate import odeint

class Integrator(object):
    """
    Container class for the various integrators to use in the parcel model.

    All defined integrators should return a tuple whose first value ``x`` is either
    ``None`` or vector containing the parcel state at all requested timestamps, and
    whose second value is a boolean indicating whether the model run was successful.

    """
    @staticmethod
    def solver(method):
        """Maps a solver name to a function.
        """
        solvers = {
            'odeint': ParcelIntegrator._solve_odeint,
            'lsoda': ParcelIntegrator._solve_lsoda,
            'lsode': ParcelIntegrator._solve_lsode,
        }

        return solvers[method]

    @staticmethod
    def _solve_lsode(f, t, y0, args, console=False, max_steps=1000):
        """Wrapper for odespy.odepack.Lsode
        """
        solver = Lsode(f, fargs=args, atol=1e-15, rtol=1e-12, nsteps=max_steps)
        solver.set_initial_condition(y0)
        solver.set(f_args=args)

        try:
            x, t = solver.solve(t)
        except ValueError, e:
            raise ValueError("something broke in LSODE")
            return None, False

        return x, True

    @staticmethod
    def _solve_lsoda(f, t, y0, args, console=False, max_steps=1000):
        """Wrapper for odespy.odepack.Lsoda
        """
        solver = Lsoda(f, fargs=args, atol=1e-15, rtol=1e-12, nsteps=max_steps)
        solver.set_initial_condition(y0)
        solver.set(f_args=args)

        try:
            x, t = solver.solve(t)
        except ValueError, e:
            raise ValueError("something broke in LSODA")
            return None, False

        return x, True

    @staticmethod
    def _solve_odeint(f, t, y0, args, console=False, max_steps=1000):
        """Wrapper for scipy.integrate.odeint
        """
        if console:
            x, info = odeint(der, y0, t, args=args,
                             full_output=1, printmessg=1, ixpr=1, mxstep=max_steps,
                             mxhnil=0, atol=1e-15, rtol=1e-12)

        else:
            x, info = odeint(der, y0, t, args=args,
                             full_output=1, mxstep=max_steps,
                             mxhnil=0, atol=1e-15, rtol=1e-12)

        success = info['message'] == "Integration successful."

        return x, success
