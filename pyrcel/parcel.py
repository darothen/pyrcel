""" Main implementation of parcel model.
"""
import os

import numpy as np
from scipy.optimize import bisect

# Parcel model imports
from . import constants as c
from . import output

# Special import of derivative ODE rhs to main namespace
from ._parcel_aux_numba import parcel_ode_sys
from .aerosol import AerosolSpecies
from .thermo import Seq, es, kohler_crit, rho_air
from .util import ParcelModelError

__all__ = ["ParcelModel"]


class ParcelModel(object):
    """ Wrapper class for instantiating and running the parcel model.

    The parcel model has been implemented in an object-oriented format to facilitate
    easy extensibility to different aerosol and meteorological conditions. A
    typical use case would involve specifying the initial conditions such as:

    >>> import pyrcel as pm
    >>> P0 = 80000.
    >>> T0 = 283.15
    >>> S0 = 0.0
    >>> V = 1.0
    >>> aerosol1 = pm.AerosolSpecies('sulfate',
    ...                              Lognorm(mu=0.025, sigma=1.3, N=2000.),
    ...                              bins=200, kappa=0.54)
    >>> initial_aerosols = [aerosol1, ]
    >>> z_top = 50.
    >>> dt = 0.01

    which initializes the model with typical conditions at the top of the boundary
    layer (800 hPa, 283.15 K, 100% Relative Humidity, 1 m/s updraft), and a simple
    sulfate aerosol distribution which will be discretized into 200 size bins to
    track. Furthermore the model was specified to simulate the updraft for 50
    meters (`z_top`) and use a time-discretization of 0.01 seconds. This
    timestep is used in the model output -- the actual ODE solver will generally
    calculate the trace of the model at many more times.

    Running the model and saving the output can be accomplished by invoking::

    >>> model = pm.ParcelModel(initial_aerosols, V, T0, S0, P0)
    >>> par_out, aer_out = pm.run(z_top, dt)

    This will yield `par_out`, a  :class:`pandas.DataFrame` containing the meteorological
    conditions in the parcel, and `aerosols`, a dictionary of :class:`DataFrame` objects
    for each species in `initial_aerosols` with the appropriately tracked size
    bins and their evolution over time.

    Attributes
    ----------
    V, T0, S0, P0, aerosols : floats
        Initial parcel settings (see **Parameters**).
    _r0s : array_like of floats
        Initial equilibrium droplet sizes.
    _r_drys : array_like of floats
        Dry radii of aerosol population.
    _kappas : array_like of floats
        Hygroscopicity of each aerosol size.
    _Nis : array_like of floats
        Number concentration of each aerosol size.
    _nr : int
        Number of aerosol sizes tracked in model.
    _model_set : boolean
        Flag indicating whether or not at any given time the model
        initialization/equilibration routine has been run with the current
        model settings.
    _y0 : array_like
        Initial state vector.

    Methods
    -------
    run(t_end, dt, max_steps=1000, solver="odeint", output_fmt="dataframes",\
        terminate=False, solver_args={})
        Execute model simulation.
    set_initial_conditions(V=None, T0=None, S0=None, P0=None, aerosols=None)
        Re-initialize a model simulation in order to run it.

    See Also
    --------
    _setup_run : companion routine which computes equilibrium droplet sizes
                 and sets the model's state vectors.

    """

    def __init__(
        self,
        aerosols,
        V,
        T0,
        S0,
        P0,
        console=False,
        accom=c.ac,
        truncate_aerosols=False,
    ):
        """Initialize the parcel model.

        Parameters
        ----------
        aerosols : array_like sequence of :class:`AerosolSpecies`
            The aerosols contained in the parcel.
        V, T0, S0, P0 : floats
            The updraft speed and initial temperature (K), pressure (Pa),
            supersaturation (percent, with 0.0 = 100% RH).
        console : boolean, optional
            Enable some basic debugging output to print to the terminal.
        accom : float, optional (default=:const:`constants.ac`)
            Condensation coefficient
        truncate_aerosols : boolean, optional (default=**False**)
            Eliminate extremely small aerosol which will cause numerical problems

        """
        self._model_set = False
        self.console = console
        self.trunc = truncate_aerosols

        self.V = V
        self.T0 = T0
        self.S0 = S0
        self.P0 = P0
        if self.trunc:
            self.aerosols = self._replace_aerosol(aerosols)
        else:
            self.aerosols = aerosols
        self.accom = accom

        # To be set by call to "self._setup_run()"
        self._r0s = None
        self._r_drys = None
        self._kappas = None
        self._Nis = None
        self._nr = None

        self._setup_run()

    def set_initial_conditions(self, V=None, T0=None, S0=None, P0=None, aerosols=None):
        """Set the initial conditions and parameters for a new parcel
        model run without having to create a new :class:`ParcelModel` instance.

        Based on the aerosol population which has been stored in the model, this
        method will finish initializing the model. This has three major parts:

        1. concatenate the aerosol population information (their dry radii,
           hygroscopicities, etc) into single arrays which can be placed into the
           state vector for forward integration.
        2. Given the initial ambient water vapor concentration (computed from the
           temperature, pressure, and supersaturation), determine how much water
           must already be coated on the aerosol particles in order for their
           size to be in equilibrium.
        3. Set-up the state vector with these initial conditions.

        Once the state vector has been set up, the setup routine will record
        attributes in the parent instance of the :class:`ParcelModel`.

        Parameters
        ----------
        V, T0, S0, P0 : floats
            The updraft speed and initial temperature (K), pressure (Pa),
            supersaturation (percent, with 0.0 = 100% RH).
        aerosols : array_like sequence of :class:`AerosolSpecies`
            The aerosols contained in the parcel.

        Raises
        ------
        ParcelModelError
            If an equilibrium droplet size distribution could not be calculated.

        Notes
        -----
        The actual setup occurs in the private method `_setup_run()`; this
        method is simply an interface that can be used to modify an existing
        :class:`ParcelModel`.

        """

        if V:
            self.V = V
        if T0:
            self.T0 = T0
        if P0:
            self.P0 = P0
        if S0:
            self.S0 = S0
        if aerosols:
            if self.trunc:
                self.aerosols = self._replace_aerosol(aerosols)
            else:
                self.aerosols = aerosols

        if T0 or P0 or S0 or aerosols:
            self._setup_run()

    def _replace_aerosol(self, aerosols, smallest_rdry=1e-10):
        """Truncates the bins with the smallest particles so that none
        have a dry radius of less than 0.1 nm.

        Parameters
        ----------
        aerosols : array_like sequence of :class:`AerosolSpecies`
            The aerosols contained in the parcel.
        smallest_rdry: float
            Smallest allowable dry radius bin.

        Returns
        -------
        fixed_aerosols : array_like sequence of :class:`AerosolSpecies`
            The modified sequence of aerosols.

        """
        fixed_aerosols = []

        for aer in aerosols:
            N_old = np.sum(aer.Nis)
            bin_count = np.count_nonzero(aer.r_drys[aer.r_drys < smallest_rdry])
            if bin_count > 0:
                aer_new = AerosolSpecies(
                    aer.species,
                    aer.distribution,
                    kappa=aer.kappa,
                    bins=aer.bins,
                    r_min=smallest_rdry * 1e6,
                )
                N_new = np.sum(aer_new.Nis)
                if self.console:
                    print(
                        "%s: Removed %03d bins, N change: %5.1f -> %5.1f (%2.1f%%)"
                        % (
                            aer.species,
                            bin_count,
                            N_old,
                            N_new,
                            (N_new - N_old) / N_new,
                        )
                    )
                fixed_aerosols.append(aer_new)
            else:
                fixed_aerosols.append(aer)

        return fixed_aerosols

    def _setup_run(self):
        """Perform equilibration and initialization calculations for the
        parcel model simulation.

        .. note:: See `set_initial_conditions` for full details.

        """
        T0, S0, P0 = self.T0, self.S0, self.P0
        z0 = 0.0
        out = dict()

        if self.console:
            print("Setting up parcel model initial conditions...")

        # 1) Setup aerosols
        # a) grab all the initial aerosol size/concentrations
        species = []
        r_drys, Nis, kappas = [], [], []
        for aerosol in self.aerosols:
            r_drys.extend(aerosol.r_drys)
            kappas.extend([aerosol.kappa] * aerosol.nr)
            Nis.extend(aerosol.Nis)
            species.extend([aerosol.species] * aerosol.nr)

        r_drys = np.array(r_drys)
        kappas = np.array(kappas)
        Nis = np.array(Nis)

        out["r_drys"] = r_drys
        out["kappas"] = kappas
        out["Nis"] = Nis

        if self.console:
            # TODO: Fix print table formatting
            print("AEROSOL DISTRIBUTION")
            print("%8s %6s" % ("r", "N"))
            for sp, r, N in zip(species, r_drys, Nis):
                print("%10s %2.2e %4.1f" % (sp, r, N))
            print("\n" + "-" * 44)

        # 2) Setup parcel initial conditions
        # a) water vapor
        #        RH    * (Rd/Rv = epsilon) * ( es / P - es )
        wv0 = (S0 + 1.0) * (
            c.epsilon * es(T0 - 273.15) / (P0 - es(T0 - 273.15))
        )  # Water Vapor mixing ratio, kg/kg

        # b) find equilibrium wet particle radius
        # wrapper function for quickly computing deviation from chosen
        # equilibrium supersaturation given current size and dry size
        f = lambda r, r_dry, kappa: Seq(r, r_dry, T0, kappa) - S0
        # Compute the equilibrium wet particle radii
        r0s = []
        # r_b, _ = kohler_crit(T0, r_drys[-1], kappas[-1])
        # for r_dry , kappa in zip(r_drys, kappas)[::-1]:
        for r_dry, kappa in zip(reversed(r_drys), reversed(kappas)):
            # Locate the critical value (f(r_crit) > 0), and use bisection from it and
            # r_dry (f(r_dry) < 0) to find equilibrium wet particle radius for a given S0
            r_b, _ = kohler_crit(T0, r_dry, kappa)
            r_a = r_dry

            r0 = bisect(f, r_a, r_b, args=(r_dry, kappa), xtol=1e-30, maxiter=500)
            r0s.append(r0)

        r0s = np.array(r0s[::-1])

        # Console logging output, if requested, of the equilibrium calcuations. Useful for
        # checking if the computations worked
        raised = False
        for r, r_dry, sp, kappa in zip(r0s, r_drys, species, kappas):
            ss = Seq(r, r_dry, T0, kappa)
            # rc, _ = kohler_crit(T0, r_dry, kappa)
            if r < 0:
                if self.console:
                    print("Found bad r", r, r_dry, sp)
                raised = True
            if np.abs(ss - S0) > 1e-4:
                if self.console:
                    print("Found S discrepancy", ss, S0, r_dry)
                raised = True
        if raised:
            raise ParcelModelError(
                "Couldn't calculate initial aerosol population wet sizes."
            )
        out["r0s"] = r0s

        # c) compute equilibrium droplet water content
        water_vol = (
            lambda r0, r_dry, Ni: (4.0 * np.pi / 3.0)
            * c.rho_w
            * Ni
            * (r0**3 - r_dry**3)
        )
        wc0 = np.sum(
            [water_vol(r0, r_dry, Ni) for r0, r_dry, Ni in zip(r0s, r_drys, Nis)]
        )
        wc0 /= rho_air(T0, P0, 0.0)

        # d) compute initial ice water content
        wi0 = 0.0

        # e) concatenate into initial conditions arrays
        y0 = [z0, P0, T0, wv0, wc0, wi0, S0]
        if self.console:
            print("PARCEL INITIAL CONDITIONS")
            print(
                ("    " + "{:>9} " * 6).format(
                    "P (hPa)",
                    "T (K)",
                    "wv (g/kg)",
                    "wc (g/kg)",
                    "wi (g/kg)",
                    "S",
                )
            )
            print(
                "    "
                + "{:9.1f} {:9.2f} {:9.1e} {:9.1e} {:9.1e} {:9.3f}".format(
                    P0 / 100.0, T0, wv0 * 1e3, wc0 * 1e3, wi0 * 1e3, S0
                )
            )
        y0.extend(r0s)
        y0 = np.array(y0)
        out["y0"] = y0
        self.y0 = y0

        # Store the model configuration
        self._r0s = r0s
        self._r_drys = r_drys
        self._kappas = kappas
        self._Nis = Nis
        self._nr = len(r_drys)

        self._model_set = True

        if self.console:
            print("Initial conditions set successfully.")

    def run(
        self,
        t_end,
        output_dt=1.0,
        solver_dt=None,
        max_steps=1000,
        solver="odeint",
        output_fmt="dataframes",
        terminate=False,
        terminate_depth=100.0,
        **solver_args
    ):
        """Run the parcel model simulation.

        Once the model has been instantiated, a simulation can immediately be
        performed by invoking this method. The numerical details underlying the
        simulation and the times over which to integrate can be flexibly set
        here.

        **Time** -- The user must specify two timesteps: `output_dt`, which is the
        timestep between output snapshots of the state of the parcel model, and
        `solver_dt`, which is the the interval of time before the ODE integrator
        is paused and re-started. It's usually okay to use a very large `solver_dt`,
        as `output_dt` can be interpolated from the simulation. In some cases though
        a small `solver_dt` could be useful to force the solver to use smaller
        internal timesteps.

        **Numerical Solver** -- By default, the model will use the `odeint` wrapper
        of LSODA shipped by default with scipy. Some fine-tuning of the solver tolerances
        is afforded here through the `max_steps`. For other solvers, a set of optional
        arguments `solver_args` can be passed.

        **Solution Output** -- Several different output formats are available by default.
        Additionally, the output arrays are saved with the `ParcelModel` instance so they
        can be used later.

        Parameters
        ----------
        t_end : float
            Total time over interval over which the model should be integrated
        output_dt : float
            Timestep intervals to report model output.
        solver_dt : float
            Timestep interval for calling solver integration routine.
        max_steps : int
            Maximum number of steps allowed by solver to satisfy error tolerances
            per timestep.
        solver : {'odeint', 'lsoda', 'lsode', 'vode', cvode'}
            Choose which numerical solver to use:
            * `'odeint'`: LSODA implementation from ODEPACK via
                SciPy's integrate module
            * `'lsoda'`: LSODA implementation from ODEPACK via odespy
            * `'lsode'`: LSODE implementation from ODEPACK via odespy
            * `'vode'` : VODE implementation from ODEPACK via odespy
            * `'cvode'` : CVODE implementation from Sundials via Assimulo
            * `'lsodar'` : LSODAR implementation from Sundials via Assimulo
        output_fmt : str, one of {'dataframes', 'arrays', 'smax'}
            Choose format of solution output.
        terminate : boolean
            End simulation at or shortly after a maximum supersaturation has been achieved
        terminate_depth : float, optional (default=100.)
            Additional depth (in meters) to integrate after termination criterion
            eached.

        Returns
        -------
        DataFrames, array, or float
            Depending on what was passed to the *output* argument, different
            types of data might be returned:

            - `dataframes': (default) will process the output into
               two pandas DataFrames - the first one containing profiles
               of the meteorological quantities tracked in the model,
               and the second a dictionary of DataFrames with one for
               each AerosolSpecies, tracking the growth in each bin
               for those species.
            - 'arrays': will return the raw output from the solver
               used internally by the parcel model - the state vector
               `y` and the evaluated timesteps converted into height
               coordinates.
            - 'smax': will only return the maximum supersaturation
               value achieved in the simulation.

        Raises
        ------
        ParcelModelError
            The parcel model failed to complete successfully or failed to initialize.

        See Also
        --------
        der : right-hand side derivative evaluated during model integration.

        """
        from .integrator import Integrator

        if output_fmt not in ["dataframes", "arrays", "smax"]:
            raise ParcelModelError(
                "Invalid value ('%s') specified for output format." % output
            )

        if solver_dt is None:
            solver_dt = 10.0 * output_dt

        if not self._model_set:
            self._setup_run()

        y0 = self.y0
        r_drys = self._r_drys
        kappas = self._kappas
        Nis = self._Nis
        nr = self._nr

        # Setup/run integrator
        try:
            # Cython
            # from .parcel_aux import der as rhs_fcn
            # Numba - JIT
            from ._parcel_aux_numba import parcel_ode_sys as rhs_fcn

            # Numba - AOT
            # from .parcel_aux_numba import parcel_ode_sys as rhs_fcn
        except ImportError:
            print("Could not load Cython derivative; using Python version.")
            from .parcel import parcel_ode_sys as rhs_fcn
        # Hack in Python derivative function
        # rhs_fcn = der

        # Is the updraft speed a function of time?
        v_is_func = hasattr(self.V, "__call__")
        if v_is_func:  # Re-wrap the function to correctly figure out V
            orig_rhs_fcn = rhs_fcn

            def rhs_fcn(y, t, *args):
                V_t = self.V(t)
                args[3] = V_t
                return orig_rhs_fcn(y, t, *args)

        # Will the simulation terminate early?
        if not terminate:
            terminate_depth = 0.0
        else:
            if terminate_depth <= 0.0:
                raise ParcelModelError("`terminate_depth` must be greater than 0!")

        if self.console:
            print()
            print("Integration control")
            print("----------------------------")
            print("        output dt: ", output_dt)
            print("    max solver dt: ", solver_dt)
            print(" solver int steps: ", int(solver_dt / output_dt))
            print("      termination: %r (%5dm)" % (terminate, terminate_depth))

        args = [nr, r_drys, Nis, self.V, kappas, self.accom]
        integrator_type = Integrator.solver(solver)
        integrator = integrator_type(
            rhs_fcn,
            output_dt,
            solver_dt,
            y0,
            args,
            terminate=terminate,
            terminate_depth=terminate_depth,
            console=self.console,
            **solver_args
        )
        success = False
        try:
            # Pack args as tuple for solvers
            args = tuple(args)

            if self.console:
                print("\nBEGIN INTEGRATION ->\n")
            x, t, success = integrator.integrate(t_end)
        except ValueError as e:
            raise ParcelModelError("Something failed during model integration: %r" % e)
        finally:
            if not success:
                raise ParcelModelError("Something failed during model integration.")

            # Success if reached this point!
            if self.console:
                print("\nEND INTEGRATION <-\n")

            self.x = x
            self.heights = self.x[:, c.STATE_VAR_MAP["z"]]
            self.time = t

            if output_fmt == "dataframes":
                return output.parcel_to_dataframes(self)
            elif output_fmt == "arrays":
                return self.x, self.heights
            elif output_fmt == "smax":
                S = self.x[:, c.STATE_VAR_MAP["S"]]
                return S.max()

    def write_summary(self, parcel_data, aerosol_data, out_filename):
        """Write a quick and dirty summary of given parcel model output to the
        terminal.
        """
        from .activation import lognormal_activation

        # Check if parent dir of out_filename exists, and if not,
        # create it
        out_dir = os.path.dirname(out_filename)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Open a file to write to
        with open(out_filename, "w") as out_file:
            # Write header
            out_file.write("PARCEL MODEL\n")
            out_file.write("--------------------------------------\n")

            # Write initial conditions
            out_file.write("V = %f\n" % self.V)
            out_file.write("P0 = %f\n" % self.P0)
            out_file.write("T0 = %f\n" % self.T0)
            out_file.write("S0 = %f\n" % self.S0)
            out_file.write("--------------------------------------\n")

            # Write aerosol details
            for aerosol in self.aerosols:
                out_file.write(aerosol.species + " = " + aerosol.summary_str() + "\n")
            out_file.write("--------------------------------------\n")

            # Write simulation summary results
            # 1) Maximum supersaturation in parcel
            S_max = parcel_data["S"].max()
            S_max_idx = np.argmax(parcel_data.S)
            out_file.write("S_max = %f\n" % S_max)

            # 2) Activated fraction of each species
            T_at_S_max = parcel_data["T"].iloc[S_max_idx]
            total_number = 0.0
            total_activated = 0.0
            for aerosol in self.aerosols:
                act_frac = lognormal_activation(
                    S_max,
                    aerosol.mu * 1e-6,
                    aerosol.sigma,
                    aerosol.N,
                    aerosol.kappa,
                    T=T_at_S_max,
                )
                act_num = act_frac * aerosol.N
                out_file.write(
                    "%s - eq_act_frac = %f (%3.2f/%3.2f)\n"
                    % (aerosol.species, act_frac, act_num, aerosol.N)
                )

                total_number += aerosol.N
                total_activated += act_num
            total_act_frac = total_activated / total_number
            out_file.write(
                "Total activated fraction = %f (%3.2f/%3.2f)\n"
                % (total_act_frac, total_activated, total_number)
            )

    def save(self, filename=None, format="nc", other_dfs=None):
        output.write_parcel_output(
            filename=filename, format=format, parcel=self, other_dfs=other_dfs
        )

    @staticmethod
    def write_csv(parcel_data, aerosol_data, output_dir=None):
        """Write output to CSV files.

        Utilize pandas fast output procedures to write the model run output to a
        set of CSV files. Written as a static method so that any prior saved
        parcel model output can be saved to disk in batch, even after its
        associated model has been destroyed.

        **Args:
            * *parcel_data* -- Pandas DataFrame of the parcel thermodynamic profile
            * *aerosol_data* -- dictionary of pandas DataFrames with the aerosol radii \
                at each model step
            * *output_dir* -- String to location where output should be saved; if not \
                provided then the model will save to the current path.

        """
        if not output_dir:
            output_dir = os.getcwd()

        # Write parcel data
        parcel_data.to_csv(os.path.join(output_dir, "parcel.csv"))

        # Write aerosol data
        for species, data in list(aerosol_data.items()):
            data.to_csv(os.path.join(output_dir, "%s.csv" % species))


#
# def parcel_ode_sys(y, t, nr, r_drys, Nis, V, kappas, accom=c.ac):
#     """ Calculates the instantaneous time-derivate of the parcel model system.
#
#     Given a current state vector `y` of the parcel model, computes the tendency
#     of each term including thermodynamic (pressure, temperature, etc) and aerosol
#     terms. The basic aerosol properties used in the model must be passed along
#     with the state vector (i.e. if being used as the callback function in an ODE
#     solver).
#
#     This function is implemented in NumPy and Python, and is likely *very* slow
#     compared to the available Cython version.
#
#     Parameters
#     ----------
#     y : array_like
#         Current state of the parcel model system,
#             * y[0] = altitude, m
#             * y[1] = Pressure, Pa
#             * y[2] = temperature, K
#             * y[3] = water vapor mass mixing ratio, kg/kg
#             * y[4] = cloud liquid water mass mixing ratio, kg/kg
#             * y[5] = cloud ice water mass mixing ratio, kg/kg
#             * y[6] = parcel supersaturation
#             * y[7:] = aerosol bin sizes (radii), m
#     t : float
#         Current simulation time, in seconds.
#     nr : Integer
#         Number of aerosol radii being tracked.
#     r_drys : array_like
#         Array recording original aerosol dry radii, m.
#     Nis : array_like
#         Array recording aerosol number concentrations, 1/(m**3).
#     V : float
#         Updraft velocity, m/s.
#     kappas : array_like
#         Array recording aerosol hygroscopicities.
#     accom : float, optional (default=:const:`constants.ac`)
#         Condensation coefficient.
#
#     Returns
#     -------
#     x : array_like
#         Array of shape (``nr``+7, ) containing the evaluated parcel model
#         instaneous derivative.
#
#     Notes
#     -----
#     This Python sketch of the derivative function shouldn't really be used for
#     any computational purposes. Instead, see the cythonized version in the auxiliary
#     file, **parcel_aux.pyx**. In the default configuration, once the code has been
#     built, you can set the environmental variable **OMP_NUM_THREADS** to control
#     the parallel for loop which calculates the condensational growth rate for each
#     bin.
#
#     """
#     from . import thermo
#
#     z = y[0]
#     P = y[1]
#     T = y[2]
#     wv = y[3]
#     wc = y[4]
#     wi = y[5]
#     S = y[6]
#     rs = np.asarray(y[c.N_STATE_VARS :])
#
#     T_c = T - 273.15  # convert temperature to Celsius
#     pv_sat = thermo.es(T - T_c)  # saturation vapor pressure
#     wv_sat = wv / (S + 1.0)  # saturation mixing ratio
#     Tv = (1.0 + 0.61 * wv) * T  # virtual temperature given parcel humidity
#     e = (1.0 + S) * pv_sat  # water vapor pressure
#     rho_air = P / (c.Rd * Tv)  # current air density accounting for humidity
#     rho_air_dry = (P - e) / c.Rd / T  # dry air density
#
#     # 1) Pressure
#     dP_dt = -1.0 * rho_air * c.g * V
#
#     # 2/3) Wet particle growth rates and droplet liquid water
#     drs_dt = np.zeros(shape=(nr,))
#     dwc_dt = 0.0
#     for i in range(nr):
#         r = rs[i]
#         r_dry = r_drys[i]
#         kappa = kappas[i]
#
#         dv_r = thermo.dv(T, r, P, accom)
#         ka_r = thermo.ka(T, rho_air, r)
#
#         G_a = (c.rho_w * c.R * T) / (pv_sat * dv_r * c.Mw)
#         G_b = (c.L * c.rho_w * ((c.L * c.Mw / (c.R * T)) - 1.0)) / (ka_r * T)
#         G = 1.0 / (G_a + G_b)
#
#         delta_S = S - thermo.Seq(r, r_dry, T, kappa)
#         dr_dt = (G / r) * delta_S
#         Ni = Nis[i]
#         dwc_dt += Ni * (r * r) * dr_dt
#         drs_dt[i] = dr_dt
#
#         # if i == nr-1:
#         #    print 1, r, r_dry, Ni
#         #    print 2, dv(T, r, P, accom), ka(T, rho_air, r)
#         #    print 3, G_a, G_b
#         #    print 4, Seq(r, r_dry, T, kappa, 1.0)
#         #    print 5, dr_dt
#
#     dwc_dt *= 4.0 * np.pi * c.rho_w / rho_air_dry
#
#     # 4) ice water content
#     dwi_dt = 0.0
#
#     # 5) Water vapor content
#     dwv_dt = -1.0 * (dwc_dt + dwi_dt)
#
#     # 6) Temperature
#     dT_dt = -c.g * V / c.Cp - c.L * dwv_dt / c.Cp
#
#     # 7) Supersaturation
#     alpha = (c.g * c.Mw * c.L) / (c.Cp * c.R * (T ** 2))
#     alpha -= (c.g * c.Ma) / (c.R * T)
#
#     gamma = (P * c.Ma) / (c.Mw * thermo.es(T - 273.15))
#     gamma += (c.Mw * c.L * c.L) / (c.Cp * c.R * T * T)
#     dS_dt = alpha * V - gamma * dwc_dt
#
#     dz_dt = V
#
#     # Repackage tendencies for feedback to numerical solver
#     x = np.zeros(shape=(nr + c.N_STATE_VARS,))
#     x[0] = dz_dt
#     x[1] = dP_dt
#     x[2] = dT_dt
#     x[3] = dwv_dt
#     x[4] = dwc_dt
#     x[5] = dwi_dt
#     x[6] = dS_dt
#     x[c.N_STATE_VARS :] = drs_dt[:]
#
#     return x
