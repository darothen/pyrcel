""" Main implementation of parcel model.
"""

import os
import time

from scipy.optimize import bisect
import numpy as np
import pandas as pd

## Parcel model imports
import constants as c
from integrator import Integrator
from thermo import *


class ParcelModelError(Exception):
    """ Custom exception to throw during parcel model execution.
    """
    def __init__(self, error_str):
        self.error_str = error_str

    def __str__(self):
        return repr(self.error_str)


class ParcelModel(object):
    """ Wrapper class for instantiating and running the parcel model.

    The parcel model has been implemented in an object-oriented format to facilitate
    easy extensibility to different aerosol and meteorological conditions. A
    typical use case would involve specifying the initial conditions such as:

    >>> import parcel_model as pm
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
    run(t_end, dt, max_steps=1000, solver="odeint", output="dataframes",\
        terminate=False, solver_args={})
        Execute model simulation.     
    set_initial_conditions(V=None, T0=None, S0=None, P0=None, aerosols=None)   
        Re-initialize a model simulation in order to run it.

    """

    def __init__(self, aerosols, V, T0, S0, P0, console=False, accom=c.ac):
        """ Initialize the parcel model.

        During the construction of the parcel model instance, an equilibrium
        droplet calculation and other basic setup routines will run. 

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
            condensation coefficient

        Returns
        -------
        self : ParcelModel
            Returns self for running the actual parcel model simulation

        See Also
        --------
        _setup_run : companion routine which computes equilibrium droplet sizes
                     and sets the model's state vectors.

        """
        self._model_set = False
        self.console = console

        self.V  = V
        self.T0 = T0
        self.S0 = S0
        self.P0 = P0
        self.aerosols = aerosols
        self.accom = accom

        ## To be set by call to "self._setup_run()"
        self._r0s    = None
        self._r_drys = None
        self._kappas = None
        self._Nis    = None
        self._nr     = None

        self._setup_run()

    def set_initial_conditions(self, V=None, T0=None, S0=None, P0=None, aerosols=None):
        """ Set the initial conditions and parameters for a new parcel
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
            self.aerosols = aerosols

        if T0 or P0 or S0 or aerosols:
            self._setup_run()

    def _setup_run(self):
        """ Perform equilibration and initialization calculations for the
        parcel model simulation.

        .. note:: See `set_initial_conditions` for full details.

        """
        T0, S0, P0 = self.T0, self.S0, self.P0
        z0 = 0.0
        out = dict()

        if self.console:
            print "Setting up parcel model initial conditions..."

        ## 1) Setup aerosols
        # a) grab all the initial aerosol size/concentrations
        species = []
        r_drys, Nis, kappas = [], [], []
        for aerosol in self.aerosols:
            r_drys.extend(aerosol.r_drys)
            kappas.extend([aerosol.kappa]*aerosol.nr)
            Nis.extend(aerosol.Nis)
            species.extend([aerosol.species]*aerosol.nr)

        r_drys = np.array(r_drys)
        kappas = np.array(kappas)
        Nis = np.array(Nis)

        out['r_drys'] = r_drys
        out['kappas'] = kappas
        out['Nis'] = Nis

        if self.console:
            print "AEROSOL DISTRIBUTION"
            print "%8s %6s" % ("r", "N")
            for sp, r, N in zip(species, r_drys, Nis):
                print "%10s %2.2e %4.1f" % (sp, r, N)
            print "\n"+"-"*44

        ## 2) Setup parcel initial conditions
        # a) water vapor
        wv0 = (1.-S0)*0.622*es(T0-273.15)/(P0-es(T0-273.15)) # Water Vapor mixing ratio, kg/kg

        # b) find equilibrium wet particle radius
        # wrapper function for quickly computing deviation from chosen
        # equilibrium supersaturation given current size and dry size
        f = lambda r, r_dry, kappa: (Seq(r, r_dry, T0, kappa) - S0)
        ## Compute the equilibrium wet particle radii
        r0s = []
        #r_b, _ = kohler_crit(T0, r_drys[-1], kappas[-1])
        for r_dry , kappa in zip(r_drys, kappas)[::-1]:
        # Locate the critical value (f(r_crit) > 0), and use bisection from it and
        # r_dry (f(r_dry) < 0) to find equilibrium wet particle radius for a given S0
            r_b, _ = kohler_crit(T0, r_dry, kappa)
            r_a = r_dry

            r0 = bisect(f, r_a, r_b, args=(r_dry, kappa), xtol=1e-30)
            r0s.append(r0)

            r_b = r0
        r0s = np.array(r0s[::-1])

        ## Console logging output, if requested, of the equilibrium calcuations. Useful for
        ## checking if the computations worked
        raised = False
        for (r,  r_dry, sp, kappa) in zip(r0s, r_drys, species, kappas):
            ss = Seq(r, r_dry, T0, kappa)
            rc, _ = kohler_crit(T0, r_dry, kappa)
            if r < 0:
                if self.console: print "Found bad r", r, r_dry, sp
                raised = True
            if np.abs(ss-S0) > 1e-8:
                if self.console: print "Found S discrepancy", ss, S0, r_dry
                raised = True
        if raised:
            raise ParcelModelError("Couldn't calculate initial aerosol population wet sizes.")
        out['r0s'] = r0s

        # c) compute equilibrium droplet water content
        wc0 = np.sum([(4.*np.pi/3.)*rho_w*Ni*(r0**3 - r_dry**3) for r0, r_dry, Ni in zip(r0s, r_drys, Nis)])

        # d) concatenate into initial conditions arrays
        y0 = [z0, P0, T0, wv0, wc0, S0]
        if self.console:
            print "PARCEL INITIAL CONDITIONS"
            print "    {:>8} {:>8} {:>8} {:>8} {:>8}".format("P (hPa)", "T (K)", "wv", "wc", "S")
            print "      %4.1f   %3.2f  %3.1e   %3.1e   %01.3f" % (P0/100., T0, wv0, wc0, S0)
        y0.extend(r0s)
        y0 = np.array(y0)
        out['y0'] = y0
        self.y0 = y0

        ## Store the model configuration
        self._r0s = r0s
        self._r_drys = r_drys
        self._kappas = kappas
        self._Nis = Nis
        self._nr = len(r_drys)

        self._model_set = True

        if self.console:
            "Initial conditions set successfully."

    def _integrate_increment(self, der_fcn, y0, dt, args, integrator, console=False, max_steps=1000, solver='odeint', **kwargs):
        """ Integrate the parcel model simulation over discrete height increments.

        .. note:: Deprecated in v1.0

        """

        V = args[3]
        t_meter = 1./V # time to travel 1 meter
        z = 0.0
        t0 = 0.0

        xs, ts = [], []
        S_old = y0[4]

        if console:
            print           "SIMULATION STATUS"
            print           "    z (m) |    T (K) |        S |  time elapsed (s)"
            status_format = " %8.3f | %8.2f | %8.5f | %8.2f/%-9.2f"

        # Loop over 1 meter increments until S begins to decrease
        S_increasing = True
        initial_time = time.time()
        while S_increasing:

            if console and z == 0:
                print status_format % (z, y0[0], S_old, 0.0, 0.0)

            t = np.arange(t0, t0+t_meter+dt, dt)

            begin = time.time()
            try:
                x, success = integrator(der_fcn, t, y0, args, console, max_steps)
            except ValueError, e:
                raise ParcelModelError("Solver '%s' failed; check console output" % solver)
            end = time.time()
            elapsed = end-begin
            total = end-initial_time

            xs.append(x[:-1])
            ts.append(t[:len(x)-1]) # If we terminate early, don't append
                                    # unused timesteps!
            t0 = t[-1]

            # Check if S is still increasing
            S_new = x[-1, 4]
            T_new = x[-1, 0]

            #print S_new, S_old
            S_increasing = S_new > S_old
            S_old = S_new

            y0 = x[-1, :].copy()

            z += 1.0
            if console:
                print status_format % (z, T_new, S_new, elapsed, total)

        x = np.concatenate(xs)
        t = np.concatenate(ts)
        return t, x

    def run(self, t_end, dt, max_steps=1000, solver="odeint", output="dataframes",
            terminate=False, solver_args={}):
        """ Run the parcel model simulation.

        Once the model has been instantiated, a simulation can immediately be 
        performed by invoking this method. The numerical details underlying the 
        simulation and the times over which to integrate can be flexibly set 
        here.

        **Time** -- The user must specify the timesteps used to evaluate the trace 
        of the parcel by passng the `dt` arg. Note that this will not necessarily dictate the 
        numerical timesteps used in the model integration, just the interpolation of
        the output (depending on which solver was requested).

        **Numerical Solver** -- By default, the model will use the `odeint` wrapper
        of LSODA shipped by default with scipy. Some fine-tuning of the solver tolerances
        is afforded here through the `max_steps`. For other solvers, a set of optional
        arguments `solver_args` can be passed as a dictionary. 

        **Solution Output** -- Several different output formats are available by default. Also,
        for the LSODA and LSODE solvers accessed via odespy, the flag `terminate` can bet
        set to *True* so that the solver will stop once it has reached a supersaturation
        maximum.

        Parameters
        ----------
        z_top : float
            Vertical distance over which the model should be integrated.
        dt : float 
            Timestep intervals to report model output. Typical value is 1/100 to 1/20
            seconds.
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
            * `'cvode'` : CVODE implementation from Assimulo
        output : {'dataframes', 'arrays', 'smax'}
            Choose format of solution output.

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
        if not output in ["dataframes", "arrays", "smax"]:
            raise ParcelModelError("Invalid value ('%s') specified for output format." % output)

        if not self._model_set:
            _ = self._setup_run()

        y0 = self.y0
        r_drys = self._r_drys
        kappas = self._kappas
        Nis = self._Nis
        nr = self._nr

        ## Setup run time conditions
        t = np.arange(0., t_end+dt, dt)
        if self.console:
            print "\n"+"n_steps = %d" % (len(t))+"\n"

        ## Setup/run integrator
        try:
            from parcel_aux import der as der_fcn
        except ImportError:
            print "Could not load Cython derivative; using Python version."
            from parcel import der as der_fcn

        args = [nr, r_drys, Nis, self.V, kappas, self.accom]
        integrator = Integrator.solver(solver)

        ## Is the updraft speed a function?
        v_is_func = hasattr(self.V, '__call__')
        if v_is_func: # Re-wrap the function to correctly figure out V            
            orig_der_fcn = der_fcn
            def der_fcn(y, t, *args):
                V_t = self.V(t)
                args[3] = V_t
                return orig_der_fcn(y, t, *args)

        try:
            ## Pack args as tuple for solvers
            args = tuple(args)

            ## Call integration routine
            x, t, success = integrator(der_fcn, t, y0, args, self.console, max_steps, terminate,
                                       **solver_args)
        except ValueError, e:
            raise ParcelModelError("Something failed during model integration: %r" % e)

        if not success:
            raise ParcelModelError("Something failed during model integration.")

        self.x = x
        self.heights = self.x[:,0]
        self.time = t

        if output == "dataframes":
            return self._convert_to_dataframes()
        elif output == "arrays":
            return self.x, self.heights
        elif output == "smax":
            S = self.x[:,5]
            return S.max()
        else: # Shouldn't ever get here; invalid output specified
            raise ParcelModelError("Invalid value (%s) specified for \
                                    output format." % output)

    def _convert_to_dataframes(self):
        """ Process the output arrays into DataFrames for returning.
        """
        x = self.x
        heights = self.heights
        time = self.time

        parcel = pd.DataFrame( {'P':x[:,1], 'T':x[:,2], 'wv':x[:,3],
                                'wc':x[:,4], 'S':x[:,5], 'z':heights},
                                 index=time)
                                 #index=heights[offset:])

        aerosol_dfs = {}
        species_shift = 0 # increment by nr to select the next aerosol's radii
        for aerosol in self.aerosols:
            nr = aerosol.nr
            species = aerosol.species

            labels = ["r%03d" % i for i in xrange(nr)]
            radii_dict = dict()
            for i, label in enumerate(labels):
                radii_dict[label] = x[:,6+species_shift+i]

            aerosol_dfs[species] = pd.DataFrame( radii_dict,
                                                 index=time  )
                                                     #index=heights[offset:])
            species_shift += nr

        return parcel, aerosol_dfs

    def write_summary(self, parcel_data, aerosol_data, out_filename):
        """ Write a quick and dirty summary of given parcel model output to the 
        terminal.
        """
        from activation import activate_lognormal_mode
        ## Check if parent dir of out_filename exists, and if not,
        ## create it
        out_dir = os.path.dirname(out_filename)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        ## Open a file to write to
        with open(out_filename, 'w') as out_file:
            ## Write header
            out_file.write("PARCEL MODEL\n")
            out_file.write("--------------------------------------\n")

            ## Write initial conditions
            out_file.write("V = %f\n" % self.V)
            out_file.write("P0 = %f\n" % self.P0)
            out_file.write("T0 = %f\n" % self.T0)
            out_file.write("S0 = %f\n" % self.S0)
            out_file.write("--------------------------------------\n")

            ## Write aerosol details
            for aerosol in self.aerosols:
                out_file.write(aerosol.species+" = "+aerosol.summary_str()+"\n")
            out_file.write("--------------------------------------\n")

            ## Write simulation summary results
            # 1) Maximum supersaturation in parcel
            S_max = parcel_data['S'].max()
            S_max_idx = np.argmax(parcel_data.S)
            out_file.write("S_max = %f\n" % S_max)

            # 2) Activated fraction of each species
            T_at_S_max = parcel_data['T'].ix[S_max_idx]
            total_number = 0.0
            total_activated = 0.0
            for aerosol in self.aerosols:
                act_frac = activate_lognormal_mode(S_max, aerosol.mu*1e-6,
                                                   aerosol.N, aerosol.kappa,
                                                   T=T_at_S_max)
                act_num = act_frac*aerosol.N
                out_file.write("%s - eq_act_frac = %f (%3.2f/%3.2f)\n" % (aerosol.species, act_frac, act_num, aerosol.N))

                total_number += aerosol.N
                total_activated += act_num
            total_act_frac = total_activated/total_number
            out_file.write("Total activated fraction = %f (%3.2f/%3.2f)\n" % (total_act_frac, total_activated, total_number))

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
        for species, data in aerosol_data.iteritems():
            data.to_csv(os.path.join(output_dir, "%s.csv" % species))

    @staticmethod
    def der(y, t, nr, r_drys, Nis, V, kappas, accom=c.ac):
        """ Calculates the instantaneous time-derivate of the parcel model system.

        Given a current state vector `y` of the parcel model, computes the tendency
        of each term including thermodynamic (pressure, temperature, etc) and aerosol
        terms. The basic aerosol properties used in the model must be passed along
        with the state vector (i.e. if being used as the callback function in an ODE
        solver).

        This function is implemented in NumPy and Python, and is likely *very* slow
        compared to the available Cython version.

        Parameters
        ----------
        y : array_like
            Current state of the parcel model system,
                * y[0] = altitude, m
                * y[1] = pressure, Pa
                * y[2] = temperature, K
                * y[3] = water vapor mass mixing ratio, kg/kg
                * y[4] = droplet liquid water mass mixing ratio, kg/kg
                * y[5] = parcel supersaturation
                * y[`nr`:] = aerosol bin sizes (radii), m
        t : float
            Current simulation time, in seconds.
        nr : Integer
            Number of aerosol radii being tracked.
        r_drys : array_like
            Array recording original aerosol dry radii, m.
        Nis : array_like
            Array recording aerosol number concentrations, 1/(m**3).
        V : float
            Updraft velocity, m/s.
        kappas : array_like
            Array recording aerosol hygroscopicities.
        accom : float, optional (default=:const:`constants.ac`)
            Condensation coefficient.

        Returns
        -------
        x : array_like
            Array of shape (`nr`+6, ) containing the evaluated parcel model
            instaneous derivative.

        Notes
        -----
        This Python sketch of the derivative function shouldn't really be used for
        any computational purposes. Instead, see the cythonized version in the auxiliary
        file, **parcel_aux.pyx**. In the default configuration, once the code has been
        built, you can set the environmental variable **OMP_NUM_THREADS** to control
        the parallel for loop which calculates the condensational growth rate for each
        bin.

        """
        z  = y[0]
        P  = y[1]
        T  = y[2]
        wv = y[3]
        wc = y[4]
        S  = y[5]
        rs = np.asarray(y[6:])

        T_c = T - 273.15   # convert temperature to Celsius
        pv_sat = es(T-T_c) # saturation vapor pressure
        wv_sat = wv/(S+1.) # saturation mixing ratio
        Tv = (1.+0.61*wv)*T # virtual temperature given parcel humidity
        rho_air = P/(Rd*Tv) # current air density accounting for humidity

        ## Calculate tendency terms
        # 1) Pressure
        dP_dt = (-g*P*V)/(Rd*Tv)

        # 2/3) Wet particle growth rates and droplet liquid water
        drs_dt = np.zeros(shape=(nr, ))
        dwc_dt = 0.
        for i in range(nr):
            r = rs[i]
            r_dry = r_drys[i]
            kappa = kappas[i]

            G_a = (rho_w*R*T)/(pv_sat*dv(T, r, P, accom)*Mw)
            G_b = (L*rho_w*((L*Mw/(R*T))-1.))/(ka(T, rho_air, r)*T)
            G = 1./(G_a + G_b)

            ## Remove size-dependence from G
            #G_a = (rho_w*R*T)/(pv_sat*Dv*Mw)
            #G_b = (L*rho_w*((L*Mw/(R*T))-1.))/(Ka*T)

            delta_S = S - Seq(r, r_dry, T, kappa, 1.0)

            dr_dt = (G/r)*delta_S

            ## ---

            Ni = Nis[i]
            dwc_dt = dwc_dt + Ni*(r**2)*dr_dt
            drs_dt[i] = dr_dt
        dwc_dt = (4.*np.pi*rho_w/rho_air)*dwc_dt

        # 4) Water vapor content
        dwv_dt = -dwc_dt

        # 5) Temperature
        dT_dt = -g*V/Cp - L*dwv_dt/Cp

        # 6) Supersaturation

        ''' Alternative methods for calculation supersaturation tendency
        # Used eq 12.28 from Pruppacher and Klett in stead of (9) from Nenes et al, 2001
        #S_a = (S+1.0)

        ## NENES (2001)
        #S_b_old = dT_dt*wv_sat*(17.67*243.5)/((243.5+(Tv-273.15))**2.)
        #S_c_old = (rho_air*g*V)*(wv_sat/P)*((0.622*L)/(Cp*Tv) - 1.0)
        #dS_dt_old = (1./wv_sat)*(dwv_dt - S_a*(S_b_old-S_c_old))

        ## PRUPPACHER (PK 1997)
        #S_b = dT_dt*0.622*L/(Rd*T**2.)
        #S_c = g*V/(Rd*T)
        #dS_dt = P*dwv_dt/(0.622*es(T-273.15)) - S_a*(S_b + S_c)

        ## SEINFELD (SP 1998)
        #S_b = L*Mw*dT_dt/(R*T**2.)
        #S_c = V*g*Ma/(R*T)
        #dS_dt = dwv_dt*(Ma*P)/(Mw*es(T-273.15)) - S_a*(S_b + S_c)
        '''

        ## GHAN (2011) - prefer to use this!
        alpha = (g*Mw*L)/(Cp*R*(T**2)) - (g*Ma)/(R*T)
        gamma = (P*Ma)/(Mw*es(T-273.15)) + (Mw*L*L)/(Cp*R*T*T)
        dS_dt = alpha*V - gamma*dwc_dt

        dz_dt = V

        ## Repackage tendencies for feedback to numerical solver
        x = np.zeros(shape=(nr+6, ))
        x[0] = dz_dt
        x[1] = dP_dt
        x[2] = dT_dt
        x[3] = dwv_dt
        x[4] = dwc_dt
        x[5] = dS_dt
        x[6:] = drs_dt[:]

        return x
