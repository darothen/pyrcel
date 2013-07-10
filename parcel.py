"""Adiabatic Cloud Parcel Model based on *Nenes et al, 2001* and implemented by
Daniel Rothenberg (darothen@mit.edu). The core logic to initialize and run the
model is contained in this script.

.. module:: parcel
    :synopsis: Parcel model main script.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""
__docformat__ = 'reStructuredText'

import os
import time

import numpy as np
from scipy.optimize import bisect
import pandas

## Parcel model imports
from aerosol import AerosolSpecies
from integrator import Integrator
from micro import *
from parcel_aux import der

class ParcelModelError(Exception):
    """Custom exception to throw during parcel model execution.
    """
    def __init__(self, error_str):
        self.error_str = error_str
    def __str__(self):
        return repr(self.error_str)

class ParcelModel(object):
    """Implementation of the logic for setting up and running the parcel model.

    The parcel model has been implemented in an object-oriented format to facilitate
    easy extensibility to different aerosol and meteorological conditions. A
    typical use case would involve specifying the initial conditions such as::

    >>> P0 = 80000.
    >>> T0 = 283.15
    >>> S0 = 0.0
    >>> V = 1.0
    >>> aerosol1 = AerosolSpecies('sulfate', Lognorm(mu=0.025, sigma=1.3, N=2000.),
                                  bins=200, kappa=0.54)
    >>> initial_aerosols = [aerosol1, ]
    >>> z_top = 50.
    >>> dt = 0.01

    which initializes the model with typical conditions at the top of the boundary
    layer (800 hPa, 283.15 K, 100% Relative Humidity, 1 m/s updraft), and a simple
    sulfate aerosol distribution which will be discretized into 200 size bins to
    track. Furthermore the model was specified to simulate the updraft for 50
    meters (``z_top``) and use a time-discretization of 0.01 seconds. This
    timestep is used in the model output -- the actual ODE solver will generally
    calculate the trace of the model at many more times.

    Running the model and saving the output can be accomplished by invoking::

    >>> pm = ParcelModel(initial_aerosols, V, T0, S0, P0)
    >>> parcel, aerosols = pm.run(z_top, dt)

    This will yield ``parcel``, a Pandas DataFrame containing the meteorological
    conditions in the parcel, and ``aerosols``, a dictionary of DataFrames
    for each species in ``initial_aerosols`` with the appropriately tracked size
    bins and their evolution over time.

    **Attributes**:
        * *aerosols* -- A list of AerosolSpecies objects representing the particle distributions to be simulated in the parcel model run.
        * *V* -- Parcel updraft speed, m/s
        * *T0* -- Parcel initial temperature, K
        * *S0* -- Parcel initial supersaturation
        * *P0* -- Parcel initial pressure, Pa
        * *console*: Flag indicating whether the model should print diagnostic output to the console while it runs.
    """

    def __init__(self, aerosols, V, T0, S0, P0, console=False):
        """
        Initialize the model with information about updraft speed and the aerosol
        distribution.
        """
        self.aerosols = aerosols
        self.V = V
        self.T0 = T0
        self.S0 = S0
        self.P0 = P0

        self.console = console

    def _setup_run(self):
        """Setup the initial parcel state for the run, given the initial
        temperature (K), pressure (Pa), and supersaturation.

        Based on the aerosol population which has been stored in the model, this
        method will finish initializing the model. This has three major parts:

        1. Concatenate the aerosol population information (their dry radii, hygroscopicities, etc) into single arrays which can be placed into the
            state vector for forward integration.
        2. Given the initial ambient water vapor concentration (computed from the
            temperature, pressure, and supersaturation), determine how much water
            must already be coated on the aerosol particles in order for their
            size to be in equilibrium.
        3. Set-up the state vector with these initial conditions.

        Once the state vector has been set up, it will return that vector as well
            as the aerosol information arrays so that they can be used in the model
            run. The ``ParcelModel`` instance itself does not save that information
            for future runs; this method is invoked on every single run of the
            parcel model.

        **Returns**:
            A dictionary containing the named arrays:

                * *r_drys* -- an array with all the aerosol dry radii concatenated
                    together
                * *kappas* -- an array with all the aerosol hygroscopicities
                    concatenated together
                * *Nis* -- an array with all the aerosol number concentrations
                    for each dry radii concatenated together
                * *y0* -- the initial state vector used in the parcel model run

        **Raises**:
            ``ParcelModelError``: An equilibrium droplet size distribution could not be calculated.
        """
        T0, S0, P0 = self.T0, self.S0, self.P0
        out = dict()

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
        for r_dry , kappa in zip(r_drys, kappas):
        # Locate the critical value (f(r_crit) > 0), and use bisection from it and
        # r_dry (f(r_dry) < 0) to find equilibrium wet particle radius for a given S0
            r_b, _ = kohler_crit(T0, r_dry, kappa)
            r_a = r_dry

            r0 = bisect(f, r_a, r_b, args=(r_dry, kappa), xtol=1e-30)
            r0s.append(r0)
        r0s = np.array(r0s)

        ## Console logging output, if requested, of the equilibrium calcuations. Useful for
        ## checking if the computations worked
        raised = False
        for (r,  r_dry, sp, kappa) in zip(r0s, r_drys, species, kappas):
            ss = Seq(r, r_dry, T0, kappa)
            rc, _ = kohler_crit(T0, r_dry, kappa)
            if r < 0:
                if self.console: print "Found bad r", r, r_dry, sp
                raised = True
            if np.abs(ss-S0) > 1e-10:
                if self.console: print "Found S discrepancy", ss, S0, r_dry
                raised = True
        if raised:
            raise ParcelModelError("Couldn't calculate initial aerosol population wet sizes.")
        out['r0s'] = r0s

        # c) compute equilibrium droplet water content
        wc0 = np.sum([(4.*np.pi/3.)*rho_w*Ni*(r0**3 - r_dry**3) for r0, r_dry, Ni in zip(r0s, r_drys, Nis)])

        # d) concatenate into initial conditions arrays
        y0 = [P0, T0, wv0, wc0, S0]
        if self.console:
            print "PARCEL INITIAL CONDITIONS"
            print "    {:>8} {:>8} {:>8} {:>8} {:>8}".format("P (hPa)", "T (K)", "wv", "wc", "S")
            print "      %4.1f   %3.2f  %3.1e   %3.1e   %01.3f" % (P0/100., T0, wv0, wc0, S0)
        y0.extend(r0s)
        y0 = np.array(y0)
        out['y0'] = y0
        self.y0 = y0

        return out

    def _integrate_increment(self, der_fcn, y0, dt, args, integrator, console=False, max_steps=1000, **kwargs):

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
            x, success = integrator(der_fcn, t, y0, args, console, max_steps)
            if not success:
                raise ParcelModelError("Solver '%s' failed; check console output" % solver)
            end = time.time()
            elapsed = end-begin
            total = end-initial_time

            xs.append(x[:-1])
            ts.append(t[:-1])
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

    def run(self, z_top, incremental=False, dt=None, ts=None, max_steps=1000, solver="odeint"):
        """Run the parcel model.

        After initializing the parcel model, it can be immediately run by
        calling this function. Before the model is integrated, a routine
        :func:`_setup_run()` is performed to equilibrate the initial aerosol
        population to the ambient meteorology. Then, the initial conditions are
        passed to a user-specified solver which integrates the system forward
        in time. By default, the integrator wraps ODEPACK's LSODA routine through
        SciPy's :func:`odeint` method, but extensions to use other solvers can be
        written easily (for instance, two methods from :mod:`odespy` are given
        below).

        The user can specify the timesteps to evaluate the trace of the parcel
        in one of two ways:

            #. Setting ``dt`` will automatically specify a timestep to use and \
                the model will use it to automatically populate the array to \
                pass to the solver.
            #. Otherwise, the user can specify ``ts`` as a list or array of the \
                timesteps where the model should be evaluated.

        **Args**:
            * *z_top* -- Vertical extent to which the model should be integrated.
            * *dt* -- Timestep intervals to report model output.
            * *ts* -- Pre-computed array of timestamps where output is requested.
            * *max_steps* -- Maximum number of steps allowed by solver to achieve \
                tolerance during integration.
            * *solver* -- A string indicating which solver to use:
                * ``"odeint"``: LSODA implementation from ODEPACK via \
                    SciPy's ``integrate`` module
                * ``"lsoda"``: LSODA implementation from ODEPACK via odespy
                * ``"lsode"``: LSODE implementation from ODEPACK via odespy

        **Returns**:
            A tuple whose first element is a DataFrame containing the profiles
            of the meteorological quantities tracked in the parcel model, and
            whose second element is a dictionary of DataFrames corresponding to
            each aerosol species and their droplet sizes.

        **Raises**:
            ``ParcelModelError``: The parcel model failed to complete successfully.

        """
        setup_results = self._setup_run()

        y0 = setup_results['y0']
        r_drys = setup_results['r_drys']
        kappas = setup_results['kappas']
        Nis = setup_results['Nis']
        nr = len(r_drys)

        ## Setup run time conditions
        if dt:
            t0 = 0.
            if self.V:
                t_end = z_top/self.V
            else:
                t_end = dt*1000
            t = np.arange(t0, t_end+dt, dt)
            if self.console:
                print "\n"+"n_steps = %d" % (len(t))+"\n"
                #raw_input("Continue run?")
        else:
            t = ts[:]

        ## Setup/run integrator
        try:
            from parcel_aux import der as der_fcn
        except ImportError:
            print "Could not load Cython derivative; using Python version."
            from parcel import der as der_fcn

        args = (nr, r_drys, Nis, self.V, kappas)
        integrator = Integrator.solver(solver)

        if incremental:
            t, x = self._integrate_increment(der_fcn, y0, dt, args, integrator, self.console, max_steps)

        else:

            x, success = integrator(der_fcn, t, y0, args, self.console, max_steps)
            if not success:
                raise ParcelModelError("Solver '%s' failed; check console output" % solver)

        heights = t*self.V
        offset = 0
        if len(heights) > x.shape[0]:
            offset = 1

        parcel = pandas.DataFrame( {'P':x[:,0], 'T':x[:,1], 'wv':x[:,2],
                                'wc':x[:,3], 'S':x[:,4]} , index=heights[offset:])

        aerosol_dfs = {}
        species_shift = 0 # increment by nr to select the next aerosol's radii
        for aerosol in self.aerosols:
            nr = aerosol.nr
            species = aerosol.species

            labels = ["r%03d" % i for i in xrange(nr)]
            radii_dict = dict()
            for i, label in enumerate(labels):
                radii_dict[label] = x[:,5+species_shift+i]

            aerosol_dfs[species] = pandas.DataFrame( radii_dict, index=heights[offset:])
            species_shift += nr

        return parcel, aerosol_dfs

    def write_summary(self, parcel_data, aerosol_data, out_filename):
        """Write a quick and dirty summary of given parcel model output to the \
        terminal.
        """
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
                act_frac = eq_act_fraction(S_max, T_at_S_max, aerosol.kappa, aerosol.r_drys, aerosol.Nis)
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

    '''
    TODO: Implement this functionality
    @staticmethod
    def write_to_hdf(name, parcel_data, aerosol_data, store_loc, meta=None):
        pass

    @staticmethod
    def retrieve_from_hdf(store_loc, run_name):
        pass
    '''

def der(y, t, nr, r_drys, Nis, V, kappas):
    """Calculates the instantaneous time-derivate of the parcel model system.

    Given a current state vector ``y`` of the parcel model, computes the tendency
    of each term including thermodynamic (pressure, temperature, etc) and aerosol
    terms. The basic aerosol properties used in the model must be passed along
    with the state vector (i.e. if being used as the callback function in an ODE
    solver).

    This function is implemented in NumPy and Python, and is likely *very* slow
    compared to the available Cython version.

    **Args**:
        * *y* -- NumPy array containing the current state of the parcel model system,
            * y[0] = pressure, Pa
            * y[1] = temperature, K
            * y[2] = water vapor mass mixing ratio, kg/kg
            * y[3] = droplet liquid water mass mixing ratio, kg/kg
            * y[3] = parcel supersaturation
            * y[nr:] = aerosol bin sizes (radii), m
        * *t* -- Current decimal model time
        * *nr* -- Integer number of aerosol radii being tracked
        * *r_drys* -- NumPy array with original aerosol dry radii, m
        * *Nis* -- NumPy array with aerosol number concentrations, m**-3
        * *V* -- Updraft velocity, m/s
        * *kappas* -- NumPy array containing all aerosol hygroscopicities

    **Returns**:
        A NumPy array with the same shape and term order as y, but containing
            all the computed tendencies at this time-step.

    """
    P = y[0]
    T = y[1]
    wv = y[2]
    wc = y[3]
    S = y[4]
    #P, T, wv, wc, S = y[:5]
    rs = np.array(y[5:])

    pv_sat = es(T-273.15) # saturation vapor pressure
    #wv_sat = wv/(S+1.) # saturation mixing ratio
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

        G_a = (rho_w*R*T)/(pv_sat*dv(T, r, P)*Mw)
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
    dT_dt = 0.
    dT_dt = -g*V/Cp - L*dwv_dt/Cp

    # 6) Supersaturation

    ## GHAN (2011) - prefer to use this!
    alpha = (g*Mw*L)/(Cp*R*(T**2)) - (g*Ma)/(R*T)
    gamma = (P*Ma)/(Mw*es(T-273.15)) + (Mw*L*L)/(Cp*R*T*T)
    dS_dt = alpha*V - gamma*dwc_dt

    ## Repackage tendencies for feedback to numerical solver
    x = np.zeros(shape=(nr+5, ))
    x[0] = dP_dt
    x[1] = dT_dt
    x[2] = dwv_dt
    x[3] = dwc_dt
    x[4] = dS_dt
    x[5:] = drs_dt[:]

    ## Kill off unused variables to get rid of numba warnings
    extra = 0.*t*wc
    if extra > 1e6:
        print "used"

    return x
