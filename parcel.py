"""Cloud Parcel Model based on Nenes et al, 2001; implemented by
Daniel Rothenberg (darothen@mit.edu) as part of research undertaken for the
General examination in the Program in Atmospheres, Oceans, and Climate

.. module:: parcel
    :synopsis: Parcel model setup and driver.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""

__docformat__ = 'reStructuredText'
import pandas
from lognorm import Lognorm, MultiModeLognorm

from pylab import *
ion()

from scipy.optimize import fsolve
from scipy.integrate import odeint

from parcel_aux import der, guesses
from micro import Seq, es, rho_w, Mw, sigma_w, R, kohler_crit

import os

class AerosolSpecies(object):
    """Container class for organizing and passing around important details about
    aerosol species present in the parcel model.

    To allow flexibility with how aerosols are defined in the model, this class is
    meant to act as a wrapper to contain metadata about aerosols (their species name, etc),
    their chemical composition (particle mass, hygroscopicity, etc), and the particular
    size distribution chosen for the initial dry aerosol. Because the latter could be very
    diverse - for instance, it might be desired to have a monodisperse aerosol population,
    or a bin representation of a canonical size distribution - the core of this class
    is designed to take those representations and homogenize them for use in the model.

    An :class:`AerosolSpecies` instance has the following attributes:

    .. attribute :: species

        A string representing a name for the particular aerosol species. This is purely
        metadata and doesn't serve any function in the parcel model except for tagging
        aerosols, so this can be set to anything.

    .. attribute :: kappa

        The hygroscopicity parameter :math:`\kappa` of the aerosol species used in
        :math:`\kappa`-Kohler theory. This should be a `float`, and is non-dimensional.

    .. attribute :: nr

        The number of bins in the size distribution. Can be 1, for a monodisperse aerosol.

    .. attribute :: r_drys

        A :mod:`numpy` array instance containing the representative dry radii for each bin
        in the aerosol size distribution. Has length equal to `nr`, and units (m).

    .. attribute :: rs

        A :mod:`numpy` array instance containing the edge of each bin in the aerosol size
        distribution. Has length equal to `nr + 1` and units (cm).

    .. attribute :: Nis

        A :mod:`numpy` array instance of length `nr` with the number concentration in
        (m**-3) of each aerosol size bin.

    .. attribute :: N

        The total aerosol species number concentration in (m**-3)

    To construct an :class:`AerosolSpecies`, only the metadata (`species` and `kappa`) and
    the size distribution needs to be specified. The size distribution (`distribution`) can be
    an instance of :class:`parcel_model.lognorm.Lognorm`, as long as an extra parameter `bins`,
    which is an `int` representing how many bins into which the distribution should be divided,
    is also passed to the constructor. In this case, the constructor will figure out how to slice
    the size distribution to calculate all the aerosol dry radii and their number concentrations.

    Alternatively, a :class:`dict` can be passed as `distribution` where that slicing has
    already occurred. In this case, `distribution` must have 2 keys: *r_drys* and *Nis*. Each
    of the values stored to those keys should fit the attribute descriptors above (although
    they don't need to be :mod:`numpy` arrays - they can be any iterable.)

    :Example:

    Constructing sulfate aerosol with a specified lognormal distribution -

    >>> aerosol1 = AerosolSpecies('(NH4)2SO4', Lognorm(mu=0.05, sigma=2.0, N=300.),
    >>>                           bins=200, kappa=0.6)

    Constructing a monodisperse sodium chloride distribution -

    >>> aerosol2 = AerosolSpecies('NaCl', {'r_drys': [0.25, ], 'Nis': [1000.0, ]}, kappa=0.2)

    .. todo:: Expand functionality for MultiModeLognorm distributions

    .. warning ::

        Throws a :class:`ValueError` if an unknown type of `distribution` is passed to the
        constructor, or if `bins` isn't present when `distribution` is an instance of
        :class:`parcel_model.lognorm.Lognorm`

    """

    def __init__(self, species, distribution, kappa, bins=None):

        self.species = species # Species molecular formula
        self.kappa = kappa # Kappa hygroscopicity parameter

        ## Handle the size distribution passed to the constructor
        self.distribution = distribution
        if isinstance(distribution, dict):
            self.r_drys = np.array(distribution['r_drys'])*1e-6
            self.rs = np.array([self.r_drys[0]*0.9, self.r_drys[0]*1.1, ])*1e6
            self.Nis = np.array(distribution['Nis'])
            self.N = np.sum(self.Nis)

        elif isinstance(distribution, Lognorm):
            # Check for missing keyword argument
            if bins is None:
                raise ValueError("Need to specify `bins` argument if passing a Lognorm distribution")

            self.mu = distribution.mu
            self.sigma = distribution.sigma
            self.N = distribution.N
            lr, rr = np.log10(self.mu/(10.*self.sigma)), np.log10(self.mu*10.*self.sigma)

            self.rs = np.logspace(lr, rr, num=bins+1)[:]
            mids = np.array([np.sqrt(a*b) for a, b in zip(self.rs[:-1], self.rs[1:])])[0:bins]
            self.Nis = np.array([0.5*(b-a)*(distribution.pdf(a) + distribution.pdf(b)) for a, b in zip(self.rs[:-1], self.rs[1:])])[0:bins]
            self.r_drys = mids*1e-6

        else:
            raise ValueError("Could not work with size distribution of type %r" % type(distribution))

        ## Correct to SI units
        # Nis: cm**-3 -> m**-3
        self.Nis *= 1e6
        self.nr = len(self.r_drys)

    def __repr__(self):
        return "%s - %r" % (self.species, self.distribution)

class ParcelModel(object):
    """Container class to set-up and run the parcel model.

    (full description of parcel model)

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

    def _setup_run(self, P0, T0, S0):
        """Setup the initial parcel state for the run, given the initial
        temperature (K), pressure (Pa), and supersaturation.

        Based on the aerosol population which has been stored in the model, this
        method will finish initializing the model. This has three major parts:

        1. Concate the aerosol population information (their dry radii, hygroscopicites,
            etc) into single arrays which can be placed into the state vector for forward
            integration.
        2. Given the initial ambient water vapor concentration (computed from the temperature,
            pressure, and supersaturation), determine how much water must already be coated
            on the aerosol particles in order for their size to be in equilibrium.
        3. Set-up the state vector with these initial conditions.

        Once the state vector has been set-up, it will return that vector as well as the
        aerosol information arrays so that they can be used in the model run. The `ParcelModel`
        instance itself does not save that information for future runs; this method is invoked
        on every single run of the parcel model.

        :param P0:
            Parcel initial pressure (Pa)
        :type P0: float

        :param T0:
            Parcel initial temperature (K)
        :type T0: float

        :param S0:
            Parcel initial supersaturation (decimal with respect to 1.0)
        :type S0: float

        :returns:
            A dictionary containinng the following arrays:

            * **r_drys** - an array with all the aerosol dry radii concatenated together
            * **kappas** - an array with all the aerosol hygroscopicities concatenated together
            * **Nis** - an array with all the aerosol number concentrations for each dry radii
              concatenated together
            * **y0** - the initial state vector used in the parcel model run

        """
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
        if self.console: print "calculating seeds for equilibrium solver..."
        r_guesses = np.array(guesses(T0, S0, r_drys, kappas))
        #r_guesses = []
        #for aerosol in self.aerosols:
        #    r_guesses.extend(guesses(T0, S0, r_drys, kappas))
        #r_guesses = np.array(r_guesses)
        if self.console: print " done"

        # wrapper function for quickly computing deviation from chosen equilibrium supersaturation given
        # current size and dry size
        f = lambda r, r_dry, kappa: (Seq(r, r_dry, T0, kappa) - S0)
        ## Compute the equilibrium wet particle radii
        r0s = np.array([fsolve(f, (rd+guess)/2., args=(rd, kappa), xtol=1e-10)[0] for guess, rd, kappa in zip(r_guesses, r_drys, kappas)])
        ## Console logging output, if requested, of the equilibrium calcuations. Useful for
        ## checking if the computations worked
        if self.console:
            for (r,  r_dry, sp, kappa) in zip(r0s, r_drys, species, kappas):
                ss = Seq(r, r_dry, T0, kappa)
                rc, _ = kohler_crit(T0, r_dry, kappa)
                if r < 0 or r > 1e-3: print "Found bad r", r, r_dry, sp
                if np.abs(ss-S0)/S0 > 0.02: print "Found S discrepancy", ss, r_dry
        out['r0s'] = r0s

        # c) compute equilibrium droplet water content
        wc0 = np.sum([(4.*np.pi/3.)*rho_w*Ni*(r0**3 - r_dry**3) for r0, r_dry, Ni in zip(r0s, r_drys, Nis)])

        # d) concat into initial conditions arrays
        y0 = [P0, T0, wv0, wc0, S0]
        if self.console:
            print "PARCEL INITIAL CONDITIONS"
            print "    {:>8} {:>8} {:>8} {:>8} {:>8}".format("P (hPa)", "T (K)", "wv", "wc", "S")
            print "      %4.1f   %3.2f  %3.1e   %3.1e   %01.2f" % (P0/100., T0, wv0, wc0, S0)
        y0.extend(r0s)
        y0 = np.array(y0)
        out['y0'] = y0
        self.y0 = y0

        return out

    def run(self, z_top, dt=None, ts=None, max_steps=1000):

        P0, T0, S0 = self.P0, self.T0, self.S0

        setup_results = self._setup_run(P0, T0, S0)
        y0 = setup_results['y0']
        r_drys = setup_results['r_drys']
        kappas = setup_results['kappas']
        Nis = setup_results['Nis']
        nr = len(r_drys)

        aerosol = self.aerosols[0]

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

        ## Setup integrator
        if self.console:
            x, info = odeint(der, y0, t, args=(nr, r_drys, Nis, self.V, kappas),
                             full_output=1, printmessg=1, ixpr=1, mxstep=max_steps,
                             mxhnil=0, atol=1e-15, rtol=1e-12)
        else:
            x = odeint(der, y0, t, args=(nr, r_drys, Nis, self.V, kappas),
                       mxstep=max_steps, mxhnil=0, atol=1e-15, rtol=1e-12)

        heights = t*self.V
        offset = 0
        if len(heights) > x.shape[0]:
            offset = 1

        df1 = pandas.DataFrame( {'P':x[:,0], 'T':x[:,1], 'wv':x[:,2],
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


        return df1, aerosol_dfs

    def write_csv(self, parcel_data, aerosol_data, output_dir=None):
        """Write output to CSV files.

        Utilize pandas fast output procedures to write the model run output to a set of CSV files.

        :param parcel_data:
            pandas DataFrame of the parcel thermodynamic properties
        :type parcel_data: pandas.DataFrame

        :param aerosol_data:
            dictionary of pandas DataFrames with the aerosol radii at each model step
        :type aerosol_data: dictionary
        """

        if not output_dir:
            output_dir = os.getcwd()

        # Write parcel data
        parcel_data.to_csv(os.path.join(output_dir, "parcel.csv"))

        # Write aerosol data
        for species, data in aerosol_data.iteritems():
            data.to_csv(os.path.join(output_dir, "%s.csv" % species))

if __name__ == "__main__":

    ## Initial conditions
    P0 = 100000. # Pressure, Pa
    T0 = 294.0 # Temperature, K
    S0 = -0.00 # Supersaturation. 1-RH from wv term
    V = 0.5 # m/s

    aerosol1 = AerosolSpecies('(NH4)2SO4', Lognorm(mu=0.05, sigma=2.0, N=100.),
                          bins=300, kappa=0.6)
    #aerosol2 = AerosolSpecies('NaCl', {'r_drys': [0.02, ], 'Nis': [1000.0, ]}, kappa=0.1)

    #initial_aerosols = [aerosol1, aerosol2]
    initial_aerosols = [aerosol1, ]
    print initial_aerosols

    aer_species = [a.species for a in initial_aerosols]
    aer_dict = dict()
    for aerosol in initial_aerosols:
        aer_dict[aerosol.species] = aerosol

    for aerosol in initial_aerosols: print np.sum(aerosol.Nis) *1e-6

    figure(2, figsize=(12,10))
    #clf()
    subplot(3,2,1)
    colors = 'bgrcmyk'
    for i, aerosol in enumerate(initial_aerosols):
        rs, Nis = aerosol.rs, aerosol.Nis
        bar(rs[:-1], Nis/np.diff(rs)*1e-6, diff(rs), color=colors[i], alpha=0.5)
    #semilogy()
    #vlines(mids, 0, ylim()[1], color='red', linestyle='dotted')
    semilogx()


    pm = ParcelModel(initial_aerosols, V, T0, S0, P0, console=True)

    ## Run model
    #rdt = np.max([V/100., 0.01])
    dt = 0.01
    parcel, aerosols = pm.run(z_top=100.0, dt=dt, max_steps=500)

    xs = np.arange(501)
    parcel = parcel.ix[parcel.index % 1 == 0]
    aero_subset = {}
    for key in aerosols:
        print key
        aerosol = aerosols[key]
        subset = aerosol.ix[aerosol.index % 1 == 0]
        aero_subset[key] = subset
    aerosols = pandas.Panel(aero_subset)

    subplot(3,2,4)
    p = parcel.S.plot(logx=False)
    print parcel.S.max()
    max_idx = np.argmax(parcel.S)
    max_z = parcel.index[max_idx]
    vlines([max_z], ylim()[0], ylim()[1], color='k', linestyle='dashed')
    xlabel("Height"); ylabel("Supersaturation")

    #subplot(3,2,2)
    #step = 1 if len(aerosols) < 20000 else 50
    #for r in aerosols[::step]:
    #    (aerosols[r]*1e6).plot(logy=True)
    #vlines([max_z], ylim()[0], ylim()[1], color='k', linestyle='dashed')
    #xlabel("Height"); ylabel("Wet Radius, micron")


    subplot(3,2,3)
    plot(parcel['T'], parcel['P']*1e-2)
    ylim(ylim()[::-1])
    hlines([parcel.P[max_idx]*1e-2], xlim()[0], xlim()[1], color='k', linestyle='dashed')
    xlabel("Temperature"); ylabel("Pressure (hPa)")

    ## PLOT AEROSOLS!!
    show()
    n_species = len(aer_species)
    fig = figure(1, figsize=(9, 5*n_species))
    clf()
    for n, key in enumerate(aerosols):
        subplot(2, n_species, n+1)
        aerosol = aerosols[key]
        for r in aerosol:
            (aerosol[r]*1e6).plot(logy=True)
        vlines([max_z], ylim()[0], ylim()[1], color='k', linestyle='dashed')
        xlabel("Height"); ylabel("Wet Radius, micron")
        title(key)

        subplot(2, n_species, n+1 + n_species)
        initial_aerosol = initial_aerosols[n]
        rs, Nis = initial_aerosol.rs, initial_aerosol.Nis
        bar(rs[:-1], Nis, diff(rs), color=colors[n], alpha=0.2)

        rs, Nis = aerosol.ix[-1]*1e6, initial_aerosol.Nis
        if not rs.shape == Nis.shape:
            new_Nis = np.zeros_like(rs)
            new_Nis[:len(Nis)] = Nis[:]
            new_Nis[len(Nis):] = 0.0
            Nis = new_Nis
        plot(rs, Nis, color='k', alpha=0.5)
        semilogx(); #semilogy()

        rs, Nis = aerosol.ix[0]*1e6, initial_aerosol.Nis
        if not rs.shape == Nis.shape:
            new_Nis = np.zeros_like(rs)
            new_Nis[:len(Nis)] = Nis[:]
            new_Nis[len(Nis):] = 0.0
            Nis = new_Nis
        plot(rs, Nis, color='r', alpha=0.5)

        #rs, Nis = pm.y0[-initial_aerosol.nr:]*1e6, initial_aerosol.Nis
        #plot(rs, Nis, color='b', alpha=0.5)


    raw_input("N analysis...")

    for species in aer_species:
        print species
        aerosol = aerosols[species]
        aer_meta = aer_dict[species]
        Nis = aer_meta.Nis

        Neq = []
        Nkn = []
        Nunact = []
        S_max = S0

        for S, T, i in zip(parcel.S, parcel['T'], xrange(len(parcel.S))):

            r_crits, s_crits = zip(*[kohler_crit(T, r_dry, aer_meta.kappa) for r_dry in aer_meta.r_drys])
            s_crits = np.array(s_crits)
            r_crits = np.array(r_crits)
            if S > S_max: S_max = S

            big_s =  S_max >= s_crits
            Neq.append(np.sum(Nis[big_s]))

            rstep = np.array(aerosol.ix[i])
            #active_radii = (S > s_crits) | (rstep > r_crits)
            active_radii = (rstep > r_crits)
            #sar = np.min(active_radii) if len(active_radii) > 0 else 1e99
            if len(active_radii) > 0:
                Nkn.append(np.sum(Nis[active_radii]))
                Nunact.append(np.sum(Nis[(rstep < r_crits)]))
            else:
                Nkn.append(0.0)
                Nunact.append(np.sum(Nis))

            print parcel.index[i], Neq[i], Nkn[i], Nunact[i], S_max, S

        Neq = np.array(Neq)
        Nkn = np.array(Nkn)
        Nunact = np.array(Nunact)

        parcel[species+'_Neq'] = Neq
        parcel[species+'_Nkn'] = Nkn
        parcel[species+'_Nunact'] = Nunact

        alphaz = Nkn/Neq
        alphaz[isnan(alphaz)] = 0.
        phiz = Nunact/Nkn
        phiz[phiz == inf] = 1.

        parcel[species+'_alpha'] = alphaz
        parcel[species+'_phi'] = phiz

        figure(2)
        ax = subplot(3,2,5)
        parcel[[species+'_Neq', species+'_Nkn']].plot(ax=ax, grid=True)
        xlabel("Height")

        subplot(3,2,6)
        parcel[species+'_alpha'].plot()
        parcel[species+'_phi'].plot()
        ylim(0, 1)
        xlabel("Height"); ylabel(r'$\alpha(z),\quad\phi(z)$')
        print alphaz[-1]

        print "=="*35
        print species + " Summary - "
        print "Max activated fraction"
        print "   Eq: ", Neq.max()/np.sum(aer_meta.Nis)
        print "  Kin: ", Nkn.max()/np.sum(aer_meta.Nis)
        print ""
        print "Alpha maximum: %2.2f" % alphaz.max()
        print "  Phi maximum: %2.2f" % phiz.max()
        print "=="*35
