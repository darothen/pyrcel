""" Collection of activation parameterizations.

"""

import numpy as np
from numpy import min as nmin
from scipy.special import erfc, erf, erfinv

from thermo import es, rho_air, ka, ka_cont, dv, dv_cont, sigma_w, kohler_crit
import constants as c

def _unpack_aerosols(aerosols):
    """ Convert a list of :class:`AerosolSpecies` into lists of aerosol properties.

    Parameters
    ----------
    aerosols : list of :class:`AerosolSpecies`

    Returns
    -------
    dictionary of lists of aerosol properties

    """

    species, mus, sigmas, kappas, Ns = [], [], [], [], []
    for a in aerosols:
        species.append(a.species)
        mus.append(a.distribution.mu)
        sigmas.append(a.distribution.sigma)
        Ns.append(a.distribution.N)
        kappas.append(a.kappa)

    species = np.asarray(species)
    mus     = np.asarray(mus)
    sigmas  = np.asarray(sigmas)
    kappas  = np.asarray(kappas)
    Ns      = np.asarray(Ns)

    return dict(species=species, mus=mus, sigmas=sigmas, Ns=Ns, kappas=kappas)

def activate_lognormal_mode(smax, mu, sigma, N, kappa, sgi=None, T=None, approx=True):
    """ Compute the activated number/fraction from a lognormal mode

    Parameters
    ----------
    smax : float
        Maximum parcel supersaturation
    mu, sigma, N : floats
        Lognormal mode parameters; ``mu`` should be in meters
    kappa : float
        Hygroscopicity of material in aerosol mode
    sgi :float, optional
        Modal critical supersaturation; if not provided, this method will
        go ahead and compute them, but a temperature ``T`` must also be passed
    T : float, optional
        Parcel temperature; only necessary if no ``sgi`` was passed
    approx : boolean, optional (default=False)
        If computing modal critical supersaturations, use the approximated
        Kohler theory

    Returns
    -------
    N_act, act_frac : floats
        Activated number concentration and fraction for the given mode

    """

    if not sgi:
        assert T
        _, sgi = kohler_crit(T, mu, kappa, approx)

    ui = 2.*np.log(sgi/smax)/(3.*np.sqrt(2.)*np.log(sigma))
    N_act = 0.5*N*erfc(ui)
    act_frac = N_act/N

    return N_act, act_frac

def multi_mode_activation(Smax, T, aerosols, rss):
    """ Compute the activation statistics of a multi-mode aerosol population.

    Parameters
    ----------
    Smax : float
        Environmental maximum supersaturation.
    T : float
        Environmental temperature.
    aerosol : array of :class:`AerosolSpecies`
        The characterizations of the dry aerosols.
    rss : array of arrays of floats
        Wet radii corresponding to each aerosol/droplet population.

    Returns
    -------
    act_fracs : floats
        The activated fraction of each aerosol population.

    """
    act_fracs = []
    for rs, aerosol in zip(rss, aerosols):
        eq, _ = act_fraction(Smax, T, rs, aerosol)
        act_fracs.append(eq)
    return act_fracs

def act_fraction(Smax, T, rs, aerosol):
    """ Calculates the equilibrium activation statistics of a given aerosol and
    maximum supersaturation.

    This calculation is only implemented for a single aerosol size. See
    :func:`multimode_act_frac` for the equivalent calculation carried out
    over successive modes.

    Parameters
    ----------
    Smax : float
        Environmental maximum supersaturation.
    T : float
        Environmental temperature.
    rs : array of floats
        Wet radii of aerosol/droplet population.
    aerosol : :class:`AerosolSpecies`
        The characterization of the dry aerosol.

    Returns
    -------
    eq_frac, kn_frac : floats
        Equilibrium and kinetic activated aerosol fraction, in units corresponding
        to those of ``aerosol``'s ``total_N`` attribute.

    """

    kappa = aerosol.kappa
    r_drys = aerosol.r_drys
    Nis = aerosol.Nis
    N_tot = np.sum(Nis)

    r_crits, s_crits = zip(*[kohler_crit(T, r_dry, kappa) for r_dry in r_drys])
    s_crits = np.array(s_crits)
    r_crits = np.array(r_crits)

    ## Equilibrium calculation - all aerosol whose critical supersaturation is
    ##                           less than the environmental supersaturation
    activated_eq = (Smax >= s_crits)
    eq_frac = np.sum(Nis[activated_eq])/N_tot

    ## Kinetic calculation - find the aerosol with the smallest dry radius which has
    ##                       grown past its critical radius, and count this aerosol and
    ##                       all larger ones as activated. This will include some large
    ##                       particles which haven't grown to critical size yet.
    is_kn_large = rs >= r_crits
    smallest_ind = np.where(is_kn_large)[0][0]
    kn_frac = np.sum(Nis[smallest_ind:])/N_tot

    return eq_frac, kn_frac

######################################################################
## Code implementing the Nenes and Seinfeld (2003) parameterization,
## with improvements from Fountoukis and Nenes (2005), Barahona
## et al (2010), and Morales Betancourt and Nenes (2014)
######################################################################

def _vpres(T):
    """ Polynomial approximation of saturated water vapour pressure as
    a function of temperature.

    Parameters
    ----------
    T : float
        Ambient temperature, in Kelvin

    Returns
    -------
    float
        Saturated water vapor pressure expressed in mb

    See Also
    --------
    es

    """
    #: Coefficients for vapor pressure approximation
    A = [ 6.107799610e+0, 4.436518521e-1, 1.428945805e-2,
          2.650648471e-4, 3.031240396e-6, 2.034080948e-8,
          6.136820929e-11 ]
    T = T-273
    vp = A[-1]*T
    for ai in reversed(A[1:-1]):
        vp = (vp + ai)*T
    vp = vp + A[0]
    return vp

def _erfp(x):
    """ Polynomial approximation to error function
    """
    AA = [0.278393, 0.230389, 0.000972, 0.078108]
    y = np.abs(1.0 * x)
    axx = 1.0 + y*(AA[0] + y*(AA[1] + y*(AA[2] + y*AA[3])))
    axx = axx*axx
    axx = axx*axx
    axx = 1.0 - (1.0/axx)

    if x <= 0.:
        axx = -1.0*axx

    return axx

def mbn2014(V, T, P, aerosols=[], accom=c.ac,
            mus=[], sigmas=[], Ns=[], kappas=[],
            xmin=1e-5, xmax=0.1, tol=1e-6, max_iters=100):
    """ Computes droplet activation using an iterative scheme.

    This method implements the iterative activation scheme under development by
    the Nenes' group at Georgia Tech. It encompasses modifications made over a
    sequence of several papers in the literature, culminating in [MBN2014]. The
    implementation here overrides some of the default physical constants and
    thermodynamic calculations to ensure consistency with a reference implementation.

    Parameters
    ----------
    V, T, P : floats
        Updraft speed (m/s), parcel temperature (K) and pressure (Pa)
    aerosols : list of :class:`AerosolSpecies`
        List of the aerosol population in the parcel; can be omitted if ``mus``,
        ``sigmas``, ``Ns``, and ``kappas`` are present. If both supplied, will
        use ``aerosols``.
    accom : float, optional (default=:const:`constants.ac`)
        Condensation/uptake accomodation coefficient
    mus, sigmas, Ns, kappas : lists of floats
        Lists of aerosol population parameters; must be present if ``aerosols``
        is not passed, but ``aerosols`` overrides if both are present
    xmin, xmax : floats, opional
        Minimum and maximum supersaturation for bisection
    tol : float, optional
        Convergence tolerance threshold for supersaturation, in decimal units
    max_iters : int, optional
        Maximum number of bisections before exiting convergence

    Returns
    -------
    smax, N_acts, act_fracs : lists of floats
        Maximum parcel supersaturation and the number concentration/activated
        fractions for each mode

    .. [MBN2014] Morales Betancourt, R. and Nenes, A.: Droplet activation
       parameterization: the population splitting concept revisited, Geosci.
       Model Dev. Discuss., 7, 2903-2932, doi:10.5194/gmdd-7-2903-2014, 2014.

    """

    if aerosols:
        d = _unpack_aerosols(aerosols)
        mus    = d['mus']
        sigmas = d['sigmas']
        kappas = d['kappas']
        Ns     = d['Ns']
    else:
        ## Assert that the aerosol was already decomposed into component vars
        assert mus
        assert sigmas
        assert Ns
        assert kappas

    # Convert sizes/number concentrations to diameters + SI units
    mus = np.asarray(mus)
    Ns = np.asarray(Ns)

    dpgs = 2*(mus*1e-6)
    Ns   = Ns*1e6

    nmodes = len(Ns)

    # Overriding using Nenes' physical constants, for numerical accuracy comparisons
    # TODO: Revert to saved 'constants' module
    c.Ma = 29e-3
    c.g  = 9.81
    c.Mw = 18e-3
    c.R  = 8.31
    c.rho_w = 1e3
    c.L  = 2.25e6
    c.Cp = 1.0061e3

    # Thermodynamic environmental values to be set
    # TODO: Revert to functions in 'thermo' module
    #rho_a    = rho_air(T, P, RH=0.0) # MBN2014: could set RH=1.0, account for moisture
    rho_a    = P*c.Ma/c.R/T
    aka      = ka_cont(T)             # MBN2014: could use thermo.ka(), include air density
    #surt     = sigma_w(T)
    surt     = (0.0761 - 1.55e-4*(T-273.))

    # Compute modal critical supersaturation (Sg) for each mode, corresponding
    # to the critical supersaturation at the median diameter of each mode
    A = 4.*c.Mw*surt/c.R/T/c.rho_w
    # There are three different ways to do this:
    #    1) original formula from MBN2014
    f = lambda T, dpg, kappa: np.sqrt((A**3.)*4./27./kappa/(dpg**3.))
    #    2) detailed kohler calculation
    #f = lambda T, dpg, kappa: kohler_crit(T, dpg/2., kappa)
    #    3) approximate kohler calculation
    #f = lambda T, dpg, kappa: kohler_crit(T, dpg/2., kappa, approx=True)
    # and the possibility of a correction factor:
    f2 = lambda T, dpg, kappa: np.exp(f(T, dpg, kappa)) - 1.0
    sgis = [f2(T, dpg, kappa) for dpg, kappa in zip(dpgs, kappas)]

    # Calculate the correction factor for the water vapor diffusivity.
    # Note that the Nenes' form is the exact same as the continuum form in this package
    #dv_orig = dv_cont(T, P)
    dv_orig = 1e-4*(0.211/(P/1.013e5)*(T/273.)**1.94)
    dv_big = 5.0e-6
    dv_low = 1e-6*0.207683*(accom**(-0.33048))

    coef     = (2.*np.pi*c.Mw/(c.R*T))**0.5
    dv_ave   = (dv_orig/(dv_big-dv_low))*((dv_big-dv_low)-(2*dv_orig/accom)*coef*      \
               (np.log((dv_big+(2*dv_orig/accom)*coef)/(dv_low+(2*dv_orig/accom)* \
                coef))))

    ## Setup constants used in supersaturation equation
    wv_pres_sat = _vpres(T)*(1e5/1e3) # MBN2014: could also use thermo.es()
    alpha = c.g*c.Mw*c.L/c.Cp/c.R/T/T - c.g*c.Ma/c.R/T
    beta1 = P*c.Ma/wv_pres_sat/c.Mw + c.Mw*c.L*c.L/c.Cp/c.R/T/T
    beta2 = c.R*T*c.rho_w/wv_pres_sat/dv_ave/c.Mw/4.0 + \
            c.L*c.rho_w/4.0/aka/T*(c.L*c.Mw/c.R/T - 1.0)   # this is 1/G
    beta  = 0.5*np.pi*beta1*c.rho_w/beta2/alpha/V/rho_a

    cf1   = 0.5*np.sqrt((1/beta2)/(alpha*V))
    cf2   = A/3.0

    def _sintegral(smax):
        """ Integrate the activation equation, using ``spar`` as the population
        splitting threshold.

        Inherits the workspace thermodynamic/constant variables from one level
        of scope higher.
        """

        zeta_c   = ((16./9.)*alpha*V*beta2*(A**2))**0.25
        delta    = 1.0 - (zeta_c/smax)**4. # spar -> smax
        critical = delta <= 0.

        if critical:
            ratio = (2e7/3.0)*A*(smax**(-0.3824) - zeta_c**(-0.3824)) # Computing sp1 and sp2 (sp1 = sp2)
            ratio = 1./np.sqrt(2.) + ratio

            if ratio > 1. : ratio = 1. # cap maximum value
            ssplt2 = smax*ratio

        else:
            ssplt1 = 0.5*(1. - np.sqrt(delta)) # min root --> sp1
            ssplt2 = 0.5*(1. + np.sqrt(delta)) # max root --> sp2
            ssplt1 = np.sqrt(ssplt1)*smax
            ssplt2 = np.sqrt(ssplt2)*smax

        ssplt = ssplt2 # secondary partitioning supersaturation

        ## Computing the condensation integrals I1 and I2
        sum_integ1 = 0.0
        sum_integ2 = 0.0

        integ1 = np.empty(nmodes)
        integ2 = np.empty(nmodes)

        sqrtwo = np.sqrt(2.)
        for i in xrange(nmodes):
            log_sigma    = np.log(sigmas[i])      # ln(sigma_i)
            log_sgi_smax = np.log(sgis[i]/smax)   # ln(sg_i/smax)
            log_sgi_sp2  = np.log(sgis[i]/ssplt2) # ln(sg_i/sp2

            u_sp2  = 2.*log_sgi_sp2/(3.*sqrtwo*log_sigma)
            u_smax = 2.*log_sgi_smax/(3.*sqrtwo*log_sigma)
            # Subtract off the integrating factor
            log_factor = 3.*log_sigma/(2.*sqrtwo)

            d_eq = A*2./sgis[i]/3./np.sqrt(3.) # Dpc/sqrt(3) - equilibrium diameter

            erf_u_sp2  = _erfp(u_sp2 - log_factor)  # ERF2
            erf_u_smax = _erfp(u_smax - log_factor) # ERF3

            integ2[i] = (np.exp(9./8.*log_sigma*log_sigma)*Ns[i]/sgis[i])*(erf_u_sp2 - erf_u_smax)

            if critical:
                u_sp_plus = sqrtwo*log_sgi_sp2/3./log_sigma
                erf_u_sp_plus = _erfp(u_sp_plus - log_factor)

                integ1[i] = 0.
                I_extra_term = Ns[i]*d_eq*np.exp((9./8.)*log_sigma*log_sigma)* \
                               (1.0 - erf_u_sp_plus)*((beta2*alpha*V)**0.5)     # 'inertially limited' particles

            else:
                g_i = np.exp((9./2.)*log_sigma*log_sigma)
                log_sgi_sp1 = np.log(sgis[i]/ssplt1)                     # ln(sg_i/sp1)

                int1_partial2  = Ns[i]*smax                              # building I1(0, sp2), eq (B4)
                int1_partial2 *= ((1. - _erfp(u_sp2)) - 0.5*((sgis[i]/smax)**2.)*g_i* \
                                  (1. - _erfp(u_sp2 + 3.*log_sigma/sqrtwo)))

                u_sp1 = 2.*log_sgi_sp1/(3.*sqrtwo*log_sigma)
                int1_partial1  = Ns[i]*smax                              # building I1(0, sp1), eq (B4)
                int1_partial1 *= ((1. - _erfp(u_sp1)) - 0.5*((sgis[i]/smax)**2.)*g_i*\
                                  (1. - _erfp(u_sp1 + 3.*log_sigma/sqrtwo)))

                integ1[i] = int1_partial2 - int1_partial1                # I1(sp1, sp2)

                u_sp1_inertial = sqrtwo*log_sgi_sp1/3./log_sigma
                erf_u_sp1 = _erfp(u_sp1_inertial - log_factor)
                I_extra_term = Ns[i]*d_eq*np.exp((9./8.)*log_sigma*log_sigma)* \
                               (1.0 - erf_u_sp1)*((beta2*alpha*V)**0.5) # 'inertially limited' particles

            ## Compute total integral values
            sum_integ1 += integ1[i] + I_extra_term
            sum_integ2 += integ2[i]

        return sum_integ1, sum_integ2

    ## Bisection routine
    # Initial calculation
    x1 = xmin # min cloud super sat -> 0.
    integ1, integ2 = _sintegral(x1)
    # Note that we ignore the contribution from the FHH integral term in this
    # implementation
    y1 = (integ1*cf1 + integ2*cf2)*beta*x1 - 1.0

    x2 = xmax # max cloud super sat
    integ1, integ2 = _sintegral(x2)
    y2 = (integ1*cf1 + integ2*cf2)*beta*x2 - 1.0

    # Iteration of bisection routine to convergence
    iter_count = 0
    for i in xrange(max_iters):
        iter_count += 1

        x3 = 0.5*(x1 + x2)
        integ1, integ2 = _sintegral(x3)
        y3 = (integ1*cf1 + integ2*cf2)*beta*x3 - 1.0

        if (y1*y3) <= 0. : # different signs
            y2 = y3
            x2 = x3
        else:
            y1 = y3
            x1 = x3
        if np.abs(x2 - x1 <= tol*x1): break
        #raw_input("continue")

    ## Finalize bisection with one more intersection
    x3 = 0.5*(x1 + x2)
    integ1, integ2 = _sintegral(x3)
    y3 = (integ1*cf1 + integ2*cf2)*beta*x3 - 1.0

    smax = x3

    n_acts, act_fracs = [], []
    for mu, sigma, N, kappa, sgi in zip(mus, sigmas, Ns, kappas, sgis):
        N_act, act_frac = activate_lognormal_mode(smax, mu*1e-6, sigma, N, kappa, sgi)
        n_acts.append(N_act)
        act_fracs.append(act_frac)

    return smax, n_acts, act_fracs

def arg2000(V, T, P, aerosols=[], accom=c.ac,
            mus=[], sigmas=[], Ns=[], kappas=[]):
    """ Computes droplet activation using a psuedo-analytical scheme.

    This method implements the psuedo-analytical scheme of [ARG2000] to
    calculate droplet activation an an adiabatically ascending parcel. It
    includes the extension to multiple lognormal modes, and the correction
    for non-unity condensation coefficient [GHAN2011].

    Parameters
    ----------
    V, T, P : floats
        Updraft speed (m/s), parcel temperature (K) and pressure (Pa)
    aerosols : list of :class:`AerosolSpecies`
        List of the aerosol population in the parcel; can be omitted if ``mus``,
        ``sigmas``, ``Ns``, and ``kappas`` are present. If both supplied, will
        use ``aerosols``.
    accom : float, optional (default=:const:`constants.ac`)
        Condensation/uptake accomodation coefficient
    mus, sigmas, Ns, kappas : lists of floats
        Lists of aerosol population parameters; must be present if ``aerosols``
        is not passed, but ``aerosols`` overrides if both are present.

    Returns
    -------
    smax, N_acts, act_fracs : lists of floats
        Maximum parcel supersaturation and the number concentration/activated
        fractions for each mode

    .. [ARG2000] Abdul-Razzak, H., and S. J. Ghan (2000), A parameterization of
       aerosol activation: 2. Multiple aerosol types, J. Geophys. Res., 105(D5),
       6837-6844, doi:10.1029/1999JD901161.

    .. [GHAN2011] Ghan, S. J. et al (2011) Droplet Nucleation: Physically-based
       Parameterization and Comparative Evaluation, J. Adv. Model. Earth Syst.,
       3, doi:10.1029/2011MS000074

    """

    if aerosols:
        d = _unpack_aerosols(aerosols)
        mus    = d['mus']
        sigmas = d['sigmas']
        kappas = d['kappas']
        Ns     = d['Ns']
    else:
        ## Assert that the aerosol was already decomposed into component vars
        assert mus
        assert sigmas
        assert Ns
        assert kappas

    ## Originally from Abdul-Razzak 1998 w/ Ma. Need kappa formulation
    wv_sat = es(T-273.15)
    alpha = (c.g*c.Mw*c.L)/(c.Cp*c.R*(T**2)) - (c.g*c.Ma)/(c.R*T)
    gamma = (c.R*T)/(wv_sat*c.Mw) + (c.Mw*(c.L**2))/(c.Cp*c.Ma*T*P)

    ## Condensation effects - base calculation
    G_a = (c.rho_w*c.R*T)/(wv_sat*dv_cont(T, P)*c.Mw)
    G_b = (c.L*c.rho_w*((c.L*c.Mw/(c.R*T))-1))/(ka_cont(T)*T)
    G_0 = 1./(G_a + G_b) # reference, no kinetic effects


    Smis = []
    Sparts = []
    for (mu, sigma, N, kappa) in zip(mus, sigmas, Ns, kappas):

        am = mu*1e-6
        N = N*1e6

        fi = 0.5*np.exp(2.5*(np.log(sigma)**2))
        gi = 1.0 + 0.25*np.log(sigma)

        A = (2.*sigma_w(T)*c.Mw)/(c.rho_w*c.R*T)
        rc_mode, Smi2 = kohler_crit(T, am, kappa, approx=True)

        ## Scale ``G`` to account for differences in condensation coefficient
        if accom == 1.0:
            G = G_0
        else:
            ## Scale using the formula from [GHAN2011]
            # G_ac - estimate using critical radius of number mode radius,
            #        and new value for condensation coefficient
            G_a = (c.rho_w*c.R*T)/(wv_sat*dv(T, rc_mode, P, accom)*c.Mw)
            G_b = (c.L*c.rho_w*((c.L*c.Mw/(c.R*T))-1))/(ka_cont(T)*T)
            G_ac = 1./(G_a + G_b)

            # G_ac1 - estimate using critical radius of number mode radius,
            #         unity condensation coefficient; re_use G_b (no change)
            G_a = (c.rho_w*c.R*T)/(wv_sat*dv(T, rc_mode, P, accom=1.0)*c.Mw)
            G_ac1 = 1./(G_a + G_b)

            # Combine using scaling formula (40) from [GHAN2011]
            G = G_0*G_ac/G_ac1

        ## Parameterization integral solutions
        zeta = (2./3.)*A*(np.sqrt(alpha*V/G))
        etai = ((alpha*V/G)**(3./2.))/(N*gamma*c.rho_w*2.*np.pi)

        ## Contributions to maximum supersaturation
        Spa = fi*((zeta/etai)**(1.5))
        Spb = gi*(((Smi2**2)/(etai + 3.*zeta))**(0.75))
        S_part = (1./(Smi2**2))*(Spa + Spb)

        Smis.append(Smi2)
        Sparts.append(S_part)

    smax = 1./np.sqrt(np.sum(Sparts))

    n_acts, act_fracs = [], []
    for mu, sigma, N, kappa, sgi in zip(mus, sigmas, Ns, kappas, Smis):
        N_act, act_frac = activate_lognormal_mode(smax, mu*1e-6, sigma, N, kappa, sgi)
        n_acts.append(N_act)
        act_fracs.append(act_frac)

    return smax, n_acts, act_fracs

def shipwayabel2010(V, T, P, aerosol):
    """ Activation scheme following Shipway and Abel, 2010 
    (doi:10.1016/j.atmosres.2009.10.005).

    """
    rho_a = rho_air(T, P)

    # The following calculation for Dv_mean is identical to the Fountoukis and Nenes (2005)
    # implementation, as referenced in Shipway and Abel, 2010
    Dp_big = 5e-6
    Dp_low = np.min([0.207683*(ac**-0.33048), 5.0])*1e-5
    Dp_B = 2.*dv_cont(T, P)*np.sqrt(2*np.pi*Mw/R/T)/ac
    Dp_diff = Dp_big - Dp_low
    Dv_mean = (dv_cont(T, P)/Dp_diff)*(Dp_diff - Dp_B*np.log((Dp_big + Dp_B)/(Dp_low+Dp_B)))

    G = 1./rho_w/(Rv*T/es(T-273.15)/Dv_mean + (L/Ka/T)*(L/Rv/T - 1))

    ### FROM APPENDIX B
    psi1 = (g/T/Rd)*(Lv/Cp/T - 1.)
    psi2 = (2.*np.pi*rho_w/rho_a) \
         * ((2.*G)**(3./2.))      \
         * (P/epsilon/es(T-273.15) + epsilon*(L**2)/Rd/(T**2)/Cp)

    Smax = 0

    act_fracs = []
    #for Smi, aerosol in zip(Smis, aerosols):
    #    ui = 2.*np.log(Smi/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.distribution.sigma))
    #    N_act = 0.5*aerosol.distribution.N*erfc(ui)
    #    act_fracs.append(N_act/aerosol.distribution.N)

    return Smax, act_fracs

def ming2006(V, T, P, aerosol):
    """Ming activation scheme.

    NOTE - right now, the variable names correspond to the FORTRAN implementation
    of the routine. Will change in the future.

    TODO: rename variables
    TODO: docstring
    TODO: extend for multiple modes.
    """

    Num = aerosol.Nis*1e-6

    RpDry = aerosol.distribution.mu*1e-6
    kappa = aerosol.kappa

    ## pre-algorithm
    ## subroutine Kohler()... calculate things from Kohler theory, particularly critical
    ## radii and supersaturations for each bin
    r_crits, s_crits = zip(*[kohler_crit(T, r_dry, kappa) for r_dry in aerosol.r_drys])

    ## subroutine CalcAlphaGamma
    alpha = (g*Mw*L)/(Cp*R*(T**2)) - (g*Ma)/(R*T)
    gamma = (R*T)/(es(T-273.15)*Mw) + (Mw*(L**2))/(Cp*Ma*T*P)

    ## re-name variables as in Ming scheme
    Dpc = 2.*np.array(r_crits)*1e6
    Dp0 = r_crits/np.sqrt(3.)
    Sc = np.array(s_crits)+1.0
    DryDp = aerosol.r_drys*2.

    ## Begin algorithm
    Smax1 = 1.0
    Smax2 = 1.1

    iter_count = 1
    while (Smax2 - Smax1) > 1e-7:
        #print "\t", iter_count, Smax1, Smax2
        Smax = 0.5*(Smax2 + Smax1)
        #print "---", Smax-1.0

        ## subroutine Grow()

        ## subroutine CalcG()
        # TODO: implement size-dependent effects on Dv, ka, using Dpc
        #G_a = (rho_w*R*T)/(es(T-273.15)*Dv_T(T)*Mw)
        G_a = (rho_w*R*T)/(es(T-273.15)*dv(T, (Dpc*1e-6)/2.)*Mw)
        #G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka_T(T)*T)
        G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka(T, 1.007e3, (Dpc*1e-6)/2.)*T)
        G = 1./(G_a + G_b) # multiply by four since we're doing diameter this time

        Smax_large = (Smax > Sc) # if(Smax>Sc(count1,count2))
        WetDp = np.zeros_like(Dpc)
        #WetDp[Smax_large] = np.sqrt(Dpc[Smax_large]**2. + 1e12*(G[Smax_large]/(alpha*V))*((Smax-.0)**2.4 - (Sc[Smax_large]-.0)**2.4))
        WetDp[Smax_large] = 1e6*np.sqrt((Dpc[Smax_large]*1e-6)**2. + (G[Smax_large]/(alpha*V))*((Smax-.0)**2.4 - (Sc[Smax_large]-.0)**2.4))

        #print Dpc
        #print WetDp/DryDp
        #print WetDp

        ## subroutine Activity()
        def Activity(dry, wet, dens, molar_weight):
            temp1 = (dry**3)*dens/molar_weight
            temp2 = ((wet**3) - (dry**3))*1e3/0.018
            act = temp2/(temp1+temp2)*np.exp(0.66/T/wet)
            #print dry[0], wet[0], dens, molar_weight, act[0]
            return act
        # Is this just the Kohler curve?
        Act = np.ones_like(WetDp)
        WetDp_large = (WetDp > 1e-5) # if(WetDp(i,j)>1e-5)
        Act[WetDp_large] = Seq(WetDp[WetDp_large]*1e-6, DryDp[WetDp_large], T, kappa) + 1.0
        #Act[WetDp_large] = Activity(DryDp[WetDp_large]*1e6, WetDp[WetDp_large], 1.7418e3, 0.132)

        #print Act

        ## subroutine Conden()

        ## subroutine CalcG()
        # TODO: implement size-dependent effects on Dv, ka, using WetDp
        #G_a = (rho_w*R*T)/(es(T-273.15)*Dv_T(T)*Mw)
        G_a = (rho_w*R*T)/(es(T-273.15)*dv(T, (WetDp*1e-6)/2.)*Mw)
        #G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka_T(T)*T)
        G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka(T, 1.3e3, (WetDp*1e-6)/2.)*T)
        G = 1./(G_a + G_b) # multiply by four since we're doing diameter this time

        WetDp_large = (WetDp > Dpc) # (WetDp(count1,count2)>Dpc(count1,count2))
        #WetDp_large = (WetDp > 0)
        f_stre = lambda x: "%12.12e" % x
        f_strf = lambda x: "%1.12f" % x
        #for i, a in enumerate(Act):
        #    if WetDp[i] > Dpc[i]:
        #        print "      ",i+1,  Act[i], f_stre(Smax-Act[i])
        CondenRate = np.sum((np.pi/2.)*1e3*G[WetDp_large]*(WetDp[WetDp_large]*1e-6)*Num[WetDp_large]*1e6*
                             (Smax-Act[WetDp_large]))

        #print iter_count, "%r %r %r" % (Smax, CondenRate, alpha*V/gamma)
        DropletNum = np.sum(Num[WetDp_large])
        ActDp = 0.0
        for i in xrange(1, len(WetDp)):
            if (WetDp[i] > Dpc[i]) and (WetDp[i-1] < Dpc[i]):
                ActDp = DryDp[i]

        ## Iteration logic
        if CondenRate < (alpha*V/gamma):
            Smax1 = Smax*1.0
        else:
            Smax2 = Smax*1.0

        iter_count += 1

    Smax = Smax-1.0

    return Smax, None


### PCE Parameterization

def lognorm_to_norm(x, mu, sigma):
    """
    Map a value from the lognormal distribution with given mu and sigma to the
    standard normal distribution with mean 0 and std 1
    """
    return (np.log(x)-mu)/sigma
    
def uni_to_norm(x, a, b):
    """                                                                                               
    Map a value from the uniform distribution [a, b] to the normal distribution                       
    """
    return np.sqrt(2.)*erfinv(2.*(x-a)/(b-a)  - 1.0)

def pce_agu_param(V, T, P, aerosols):

    Smaxes = []
    for aerosol in aerosols:
        N = aerosol.distribution.N
        mu = aerosol.distribution.mu
        sigma = aerosol.distribution.sigma
        kappa = aerosol.kappa

        Smax = _pce_fit(N, mu, sigma, kappa, V, T, P)
        Smaxes.append(Smax)

        #print "PCE with", N, mu, sigma, kappa, V, T, P, Smax

    min_smax = nmin(Smaxes)
    if 0. <= min_smax <= 0.5: 
        Smax = min_smax
    else:
        return 0., [0.]*len(aerosols)

    ## Compute scrit of each mode
    scrits = []
    for aerosol in aerosols:
        _, scrit = kohler_crit(T, aerosol.distribution.mu*1e-6, aerosol.kappa)
        scrits.append(scrit)

    act_fracs = []
    for aerosol, scrit in zip(aerosols, scrits):
        ui = 2.*np.log(scrit/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.distribution.sigma))
        N_act = 0.5*aerosol.distribution.N*erfc(ui)
        act_fracs.append(N_act/aerosol.distribution.N)

    return Smax, act_fracs

def _pce_fit(N, mu, sigma, kappa, V, T, P):
    ## P in Pa
    dist_bounds = {
        'mu': [0.01, 0.25],
        #'N': [100., 10000.],
        'kappa': [0.1, 1.2],
        'sigma': [1.2, 3.0],
        'V': [0., 4.0],
        'T': [235., 310.],
        'P': [50000., 105000.],
    }
    dist_params = {
        'N': [ 7.5, 1.95 ], 
    }

    N = lognorm_to_norm(N, *dist_params['N'])

    mu = uni_to_norm(mu, *dist_bounds['mu'])
    kappa = uni_to_norm(kappa, *dist_bounds['kappa'])
    sigma = uni_to_norm(sigma, *dist_bounds['sigma'])
    V = uni_to_norm(V, *dist_bounds['V'])
    T = uni_to_norm(T, *dist_bounds['T'])
    P = uni_to_norm(P, *dist_bounds['P'])

    Smax = 6.4584111537e-6*N**6 - \
         2.5994976288e-5*N**5 - \
         1.7065251097e-7*N**4*P**2 + \
         1.3741352226e-5*N**4*P + \
         2.8567989557e-5*N**4*T**2 - \
         7.4876643038e-5*N**4*T - \
         2.0388391982e-6*N**4*V**2 + \
         4.3054466907e-5*N**4*V + \
         3.6504788687e-6*N**4*kappa**2 + \
         8.7165631487e-7*N**4*kappa + \
         1.6542743001e-5*N**4*mu**2 + \
         4.8195946039e-6*N**4*mu + \
         3.9282682647e-6*N**4*sigma**2 + \
         1.1137326431e-5*N**4*sigma + \
         2.795758112227e-5*N**4 + \
         1.5947545697e-6*N**3*P**2 - \
         6.9358311166e-5*N**3*P - \
         0.00014252420422*N**3*T**2 + \
         0.00039466661884*N**3*T + \
         2.15368184e-5*N**3*V**2 - \
         0.00025279065671*N**3*V + \
         4.6142483833e-6*N**3*kappa**2 - \
         2.5055687574e-5*N**3*kappa - \
         3.0424806654e-6*N**3*mu**2 - \
         4.5156027497e-5*N**3*mu - \
         1.780917608e-6*N**3*sigma**2 - \
         2.516400813e-5*N**3*sigma - \
         0.0003567127574296*N**3 + \
         5.9696014699e-7*N**2*P**4 - \
         1.3472490172e-5*N**2*P**3 - \
         1.0610551852e-6*N**2*P**2*T**2 + \
         2.0181530448e-6*N**2*P**2*T + \
         2.5327194907e-7*N**2*P**2*V**2 - \
         1.4006527233e-6*N**2*P**2*V + \
         5.4851851852e-7*N**2*P**2*kappa**2 - \
         1.320380981e-6*N**2*P**2*kappa + \
         1.7644666667e-7*N**2*P**2*mu**2 - \
         2.7894950894e-7*N**2*P**2*mu + \
         1.8201189815e-7*N**2*P**2*sigma**2 - \
         5.0510811394e-7*N**2*P**2*sigma - \
         6.88818634103e-6*N**2*P**2 + \
         5.0207581099e-5*N**2*P*T**2 - \
         0.00013814911722*N**2*P*T - \
         6.2792651121e-6*N**2*P*V**2 + \
         7.2980075931e-5*N**2*P*V - \
         3.7856114614e-6*N**2*P*kappa**2 + \
         1.2860228333e-5*N**2*P*kappa - \
         1.5691902399e-6*N**2*P*mu**2 + \
         8.2376491667e-6*N**2*P*mu - \
         1.3435745045e-6*N**2*P*sigma**2 + \
         6.0282465278e-6*N**2*P*sigma + \
         0.0001877522259389*N**2*P - \
         4.0442507595e-5*N**2*T**4 + \
         5.6586533058e-5*N**2*T**3 - \
         8.9548419306e-6*N**2*T**2*V**2 + \
         0.00014183762216*N**2*T**2*V + \
         1.7477041667e-7*N**2*T**2*kappa**2 - \
         2.2336680774e-5*N**2*T**2*kappa - \
         3.9516949861e-5*N**2*T**2*mu**2 - \
         1.428384236e-5*N**2*T**2*mu - \
         8.1085041667e-6*N**2*T**2*sigma**2 - \
         4.4004842538e-5*N**2*T**2*sigma + \
         0.00038258884934483*N**2*T**2 + \
         2.9970384599e-5*N**2*T*V**2 - \
         0.00041049796829*N**2*T*V + \
         6.5092115599e-6*N**2*T*kappa**2 + \
         3.0809800694e-5*N**2*T*kappa + \
         9.9551207477e-5*N**2*T*mu**2 + \
         1.0952167639e-5*N**2*T*mu + \
         2.1329980047e-5*N**2*T*sigma**2 + \
         8.81912525e-5*N**2*T*sigma - \
         0.0008911162845737*N**2*T + \
         6.9026802931e-6*N**2*V**4 - \
         4.7531336217e-5*N**2*V**3 + \
         2.5832318241e-6*N**2*V**2*kappa**2 - \
         1.2472907784e-6*N**2*V**2*kappa + \
         1.1149875079e-5*N**2*V**2*mu**2 - \
         2.9708960501e-6*N**2*V**2*mu + \
         3.2880035648e-7*N**2*V**2*sigma**2 + \
         7.9685785603e-6*N**2*V**2*sigma - \
         8.857197689645e-5*N**2*V**2 - \
         1.3905780926e-5*N**2*V*kappa**2 + \
         2.6425726833e-5*N**2*V*kappa - \
         4.4290453362e-5*N**2*V*mu**2 + \
         3.4602470958e-5*N**2*V*mu - \
         9.497372933e-6*N**2*V*sigma**2 - \
         8.4509070972e-6*N**2*V*sigma + \
         0.0007493009795633*N**2*V + \
         2.8884698866e-6*N**2*kappa**4 + \
         1.349739092e-6*N**2*kappa**3 + \
         1.7550156389e-5*N**2*kappa**2*mu**2 - \
         1.9786638902e-6*N**2*kappa**2*mu + \
         5.529520787e-6*N**2*kappa**2*sigma**2 + \
         1.2209966835e-5*N**2*kappa**2*sigma - \
         5.448370112109e-5*N**2*kappa**2 - \
         4.359847391e-5*N**2*kappa*mu**2 - \
         2.2737228056e-5*N**2*kappa*mu - \
         9.990113266e-6*N**2*kappa*sigma**2 - \
         5.9185131528e-5*N**2*kappa*sigma + \
         9.22206763018e-6*N**2*kappa + \
         3.0424263183e-5*N**2*mu**4 - \
         1.9098455668e-5*N**2*mu**3 + \
         8.54937625e-6*N**2*mu**2*sigma**2 - \
         4.684071842e-6*N**2*mu**2*sigma - \
         0.00035110649314667*N**2*mu**2 - \
         3.6261121147e-6*N**2*mu*sigma**2 - \
         9.2769369028e-5*N**2*mu*sigma + \
         0.00011212992202954*N**2*mu + \
         5.0891009441e-6*N**2*sigma**4 + \
         3.5893477645e-6*N**2*sigma**3 - \
         7.197212424173e-5*N**2*sigma**2 - \
         0.00011060069230486*N**2*sigma + \
         0.00151313669719111*N**2 - \
         1.1284287469e-6*N*P**4 + \
         3.0704412322e-5*N*P**3 + \
         1.6278653855e-6*N*P**2*T**2 - \
         3.3672619444e-6*N*P**2*T - \
         6.2110532065e-8*N*P**2*V**2 + \
         3.2172427639e-6*N*P**2*V - \
         2.8443321387e-7*N*P**2*kappa**2 + \
         7.4341916667e-7*N*P**2*kappa - \
         7.2252756038e-7*N*P**2*mu**2 + \
         8.4614527778e-7*N*P**2*mu - \
         1.2720654237e-7*N*P**2*sigma**2 + \
         2.7419097222e-7*N*P**2*sigma + \
         8.764064001385e-6*N*P**2 - \
         9.5804124167e-5*N*P*T**2 + \
         0.00027775604478*N*P*T + \
         1.2225588236e-5*N*P*V**2 - \
         0.0001716045343*N*P*V + \
         1.8806313889e-6*N*P*kappa**2 - \
         1.2701263287e-6*N*P*kappa + \
         1.3923449444e-5*N*P*mu**2 - \
         1.4186857698e-6*N*P*mu + \
         3.9634190278e-6*N*P*sigma**2 + \
         9.947051115e-6*N*P*sigma - \
         0.0003223815596977*N*P + \
         7.5931315992e-5*N*T**4 - \
         0.00011373913927*N*T**3 + \
         1.5275617964e-5*N*T**2*V**2 - \
         0.00027033311532*N*T**2*V - \
         9.5487976257e-6*N*T**2*kappa**2 + \
         6.1326942361e-5*N*T**2*kappa + \
         1.1264628419e-5*N*T**2*mu**2 + \
         0.00013788716208*N*T**2*mu + \
         8.4656605352e-6*N*T**2*sigma**2 + \
         9.7041865833e-5*N*T**2*sigma - \
         0.00046604057151*N*T**2 - \
         5.2633321153e-5*N*T*V**2 + \
         0.00082461082486*N*T*V + \
         1.9930343472e-5*N*T*kappa**2 - \
         0.00014021885993*N*T*kappa - \
         3.2740759028e-5*N*T*mu**2 - \
         0.00033633604715*N*T*mu - \
         2.7467490833e-5*N*T*sigma**2 - \
         0.0002267701814*N*T*sigma + \
         0.0010623078872764*N*T - \
         1.158262868e-5*N*V**4 + \
         0.00010282581987*N*V**3 + \
         2.0236346649e-7*N*V**2*kappa**2 - \
         7.0861126667e-6*N*V**2*kappa - \
         1.8571179464e-7*N*V**2*mu**2 - \
         2.9025647069e-5*N*V**2*mu - \
         2.8250510694e-6*N*V**2*sigma**2 - \
         1.0873061236e-5*N*V**2*sigma + \
         9.5035330545615e-5*N*V**2 - \
         4.7114408333e-7*N*V*kappa**2 + \
         3.5583995899e-5*N*V*kappa + \
         4.3931099764e-5*N*V*mu**2 + \
         9.4893199047e-5*N*V*mu + \
         1.9266056153e-5*N*V*sigma**2 + \
         8.172457216e-5*N*V*sigma - \
         0.00114623544365757*N*V + \
         2.5465455757e-6*N*kappa**4 - \
         1.1844938245e-5*N*kappa**3 - \
         1.2548361851e-5*N*kappa**2*mu**2 - \
         6.6498102778e-6*N*kappa**2*mu - \
         2.6413318548e-6*N*kappa**2*sigma**2 - \
         1.9743217083e-5*N*kappa**2*sigma - \
         5.237116508322e-5*N*kappa**2 + \
         2.2628180833e-5*N*kappa*mu**2 + \
         6.4409460169e-5*N*kappa*mu + \
         2.1659551389e-6*N*kappa*sigma**2 + \
         8.2294682962e-5*N*kappa*sigma + \
         0.00031141975514413*N*kappa - \
         1.1565474947e-5*N*mu**4 - \
         2.1450508636e-5*N*mu**3 + \
         6.1585477702e-6*N*mu**2*sigma**2 - \
         8.815663375e-5*N*mu**2*sigma + \
         8.443778184842e-5*N*mu**2 - \
         3.4406999306e-5*N*mu*sigma**2 + \
         0.00027943018423*N*mu*sigma + \
         0.00073132224303402*N*mu - \
         8.4378328798e-6*N*sigma**4 - \
         2.0928942447e-5*N*sigma**3 + \
         9.361717372097e-5*N*sigma**2 + \
         0.00042563590086478*N*sigma - \
         0.00207631579223133*N - \
         1.9577562243e-8*P**6 + \
         9.8981784049e-7*P**5 - \
         5.8363352597e-8*P**4*T**2 + \
         2.6457614122e-8*P**4*T + \
         9.0459993866e-8*P**4*V**2 + \
         2.1439092975e-7*P**4*V + \
         1.7814328446e-7*P**4*kappa**2 - \
         4.1686622901e-7*P**4*kappa + \
         5.3644855238e-8*P**4*mu**2 - \
         2.0156224591e-7*P**4*mu + \
         5.8210558734e-8*P**4*sigma**2 - \
         1.8962978248e-7*P**4*sigma + \
         4.46780827354e-7*P**4 - \
         1.1972281072e-6*P**3*T**2 + \
         5.3768532472e-6*P**3*T + \
         2.2417961995e-7*P**3*V**2 - \
         6.4936735747e-6*P**3*V - \
         7.4042040112e-7*P**3*kappa**2 + \
         3.1988643743e-6*P**3*kappa - \
         1.1174867493e-7*P**3*mu**2 + \
         4.9097886778e-6*P**3*mu + \
         2.5563905537e-7*P**3*sigma**2 + \
         3.0112545186e-6*P**3*sigma - \
         2.528028422697e-5*P**3 - \
         6.2461608542e-8*P**2*T**4 + \
         7.2160816962e-9*P**2*T**3 - \
         3.9231712963e-8*P**2*T**2*V**2 + \
         1.2045330835e-7*P**2*T**2*V - \
         1.2065046296e-7*P**2*T**2*kappa**2 + \
         5.5655764074e-8*P**2*T**2*kappa - \
         7.1769768519e-7*P**2*T**2*mu**2 + \
         9.214390015e-7*P**2*T**2*mu - \
         1.1352175926e-7*P**2*T**2*sigma**2 - \
         7.2727690784e-8*P**2*T**2*sigma + \
         9.20245283507e-7*P**2*T**2 + \
         8.0179117696e-8*P**2*T*V**2 + \
         4.235625e-8*P**2*T*V + \
         2.258923022e-7*P**2*T*kappa**2 - \
         4.446875e-7*P**2*T*kappa + \
         1.319233337e-6*P**2*T*mu**2 - \
         2.0462986111e-6*P**2*T*mu + \
         8.8850197051e-8*P**2*T*sigma**2 + \
         3.2751388889e-8*P**2*T*sigma - \
         3.083990086676e-7*P**2*T + \
         1.3373206908e-7*P**2*V**4 + \
         9.0363593389e-8*P**2*V**3 + \
         2.4463467593e-7*P**2*V**2*kappa**2 - \
         3.2486785978e-7*P**2*V**2*kappa + \
         2.7125416667e-7*P**2*V**2*mu**2 - \
         7.8996752457e-8*P**2*V**2*mu + \
         2.3869884259e-7*P**2*V**2*sigma**2 - \
         3.9030882886e-8*P**2*V**2*sigma - \
         1.676552162853e-6*P**2*V**2 - \
         3.2961961287e-7*P**2*V*kappa**2 + \
         1.0572459722e-6*P**2*V*kappa + \
         5.8846025249e-7*P**2*V*mu**2 - \
         1.2420694444e-7*P**2*V*mu + \
         1.1084042637e-7*P**2*V*sigma**2 + \
         5.9708680556e-7*P**2*V*sigma - \
         2.731802676607e-6*P**2*V + \
         1.6946466336e-7*P**2*kappa**4 - \
         1.697455568e-7*P**2*kappa**3 + \
         4.3087824074e-7*P**2*kappa**2*mu**2 - \
         2.6231749106e-8*P**2*kappa**2*mu + \
         4.1886990741e-7*P**2*kappa**2*sigma**2 - \
         2.1180014438e-7*P**2*kappa**2*sigma - \
         2.68356980142e-6*P**2*kappa**2 - \
         1.1445592205e-6*P**2*kappa*mu**2 + \
         2.838125e-7*P**2*kappa*mu - \
         7.6035506889e-7*P**2*kappa*sigma**2 + \
         3.2736111111e-9*P**2*kappa*sigma + \
         5.198700314056e-6*P**2*kappa + \
         2.570161515e-7*P**2*mu**4 - \
         4.2080255264e-7*P**2*mu**3 + \
         1.7831666667e-7*P**2*mu**2*sigma**2 - \
         1.3384887703e-6*P**2*mu**2*sigma - \
         2.364220102728e-6*P**2*mu**2 + \
         8.8640907579e-8*P**2*mu*sigma**2 + \
         1.2368444444e-6*P**2*mu*sigma + \
         4.855707265804e-6*P**2*mu + \
         9.059626753e-8*P**2*sigma**4 - \
         1.5254956604e-7*P**2*sigma**3 - \
         1.126043894964e-6*P**2*sigma**2 + \
         2.8315063203e-6*P**2*sigma - \
         2.778055205468e-6*P**2 - \
         2.1400747816e-6*P*T**4 + \
         4.9258105417e-6*P*T**3 + \
         5.6327936107e-7*P*T**2*V**2 + \
         9.3250965278e-6*P*T**2*V + \
         3.2710316757e-6*P*T**2*kappa**2 - \
         8.2942513889e-6*P*T**2*kappa + \
         7.9343747988e-6*P*T**2*mu**2 - \
         1.9181326389e-5*P*T**2*mu + \
         1.5642303199e-6*P*T**2*sigma**2 - \
         7.5503763889e-6*P*T**2*sigma + \
         2.871260752173e-5*P*T**2 + \
         1.4432569444e-7*P*T*V**2 - \
         3.6214883811e-5*P*T*V - \
         7.5203763889e-6*P*T*kappa**2 + \
         2.3241170134e-5*P*T*kappa - \
         1.7993693056e-5*P*T*mu**2 + \
         5.5177810267e-5*P*T*mu - \
         1.9631597222e-6*P*T*sigma**2 + \
         2.2474684728e-5*P*T*sigma - \
         0.00010697897078404*P*T + \
         4.3918577044e-7*P*V**4 - \
         6.5398867255e-6*P*V**3 - \
         1.6642529664e-6*P*V**2*kappa**2 + \
         4.5294820833e-6*P*V**2*kappa - \
         8.3413225416e-7*P*V**2*mu**2 + \
         3.3331961111e-6*P*V**2*mu - \
         5.4879252019e-7*P*V**2*sigma**2 + \
         2.5988245833e-6*P*V**2*sigma - \
         6.06167138571e-6*P*V**2 + \
         1.76840125e-6*P*V*kappa**2 - \
         1.2335091147e-5*P*V*kappa - \
         7.172135e-6*P*V*mu**2 - \
         1.4280271047e-5*P*V*mu - \
         1.7709023611e-6*P*V*sigma**2 - \
         1.5558439144e-5*P*V*sigma + \
         0.0001341572007329*P*V - \
         1.3847755834e-6*P*kappa**4 + \
         2.9663988051e-6*P*kappa**3 - \
         6.701873906e-7*P*kappa**2*mu**2 + \
         1.2602319444e-6*P*kappa**2*mu - \
         1.6684083648e-6*P*kappa**2*sigma**2 + \
         4.8915402778e-6*P*kappa**2*sigma + \
         1.830572029266e-5*P*kappa**2 + \
         7.2758208333e-6*P*kappa*mu**2 - \
         4.4294673496e-6*P*kappa*mu + \
         5.7967319444e-6*P*kappa*sigma**2 - \
         7.5561654674e-6*P*kappa*sigma - \
         6.41037679493e-5*P*kappa + \
         1.1303131108e-7*P*mu**4 + \
         4.477378242e-6*P*mu**3 - \
         8.1756005619e-8*P*mu**2*sigma**2 + \
         1.8019847222e-5*P*mu**2*sigma + \
         3.591958513889e-6*P*mu**2 + \
         3.5776194444e-6*P*mu*sigma**2 - \
         2.1832827595e-5*P*mu*sigma - \
         0.000104836086236*P*mu + \
         4.4303746065e-7*P*sigma**4 + \
         3.5105365953e-6*P*sigma**3 - \
         5.903007579301e-6*P*sigma**2 - \
         6.46457674367e-5*P*sigma + \
         0.000248152385516871*P + \
         9.3486231047e-7*T**6 - \
         1.8168892496e-6*T**5 - \
         1.2017117971e-7*T**4*V**2 - \
         6.0623117829e-6*T**4*V - \
         2.2138037142e-6*T**4*kappa**2 + \
         5.3907485972e-6*T**4*kappa - \
         1.5012731379e-5*T**4*mu**2 + \
         2.9320254519e-5*T**4*mu - \
         9.9591672455e-7*T**4*sigma**2 + \
         4.5407528726e-6*T**4*sigma - \
         2.1701625406048e-5*T**4 - \
         8.86080478e-8*T**3*V**2 + \
         1.4786332237e-5*T**3*V + \
         3.9370511477e-6*T**3*kappa**2 - \
         1.1123904654e-5*T**3*kappa + \
         1.9629299957e-5*T**3*mu**2 - \
         4.3941049344e-5*T**3*mu + \
         1.1089594981e-6*T**3*sigma**2 - \
         9.8896573442e-6*T**3*sigma + \
         4.92291899313038e-5*T**3 - \
         1.8208011851e-7*T**2*V**4 - \
         3.3830247075e-6*T**2*V**3 + \
         1.5875634259e-7*T**2*V**2*kappa**2 - \
         8.2594310191e-7*T**2*V**2*kappa - \
         7.0372356481e-6*T**2*V**2*mu**2 + \
         1.1314601854e-5*T**2*V**2*mu + \
         2.6109148148e-7*T**2*V**2*sigma**2 - \
         1.4817938429e-6*T**2*V**2*sigma + \
         2.584581978913e-6*T**2*V**2 + \
         8.1820206167e-6*T**2*V*kappa**2 - \
         2.152989125e-5*T**2*V*kappa + \
         3.5087114657e-5*T**2*V*mu**2 - \
         7.62298625e-5*T**2*V*mu + \
         3.7441553065e-6*T**2*V*sigma**2 - \
         1.9480860556e-5*T**2*V*sigma + \
         8.078026266135e-5*T**2*V - \
         1.4079168102e-6*T**2*kappa**4 + \
         8.494577264e-7*T**2*kappa**3 - \
         1.9940069444e-6*T**2*kappa**2*mu**2 - \
         3.7981316228e-6*T**2*kappa**2*mu + \
         1.2927453704e-7*T**2*kappa**2*sigma**2 - \
         5.9931564037e-6*T**2*kappa**2*sigma + \
         2.918584228186e-5*T**2*kappa**2 - \
         5.4461049955e-7*T**2*kappa*mu**2 + \
         1.7234859722e-5*T**2*kappa*mu - \
         1.9437451044e-6*T**2*kappa*sigma**2 + \
         1.7031693056e-5*T**2*kappa*sigma - \
         6.0924394269614e-5*T**2*kappa - \
         1.3344139356e-5*T**2*mu**4 + \
         6.9440170146e-6*T**2*mu**3 - \
         6.7769083333e-6*T**2*mu**2*sigma**2 + \
         3.2599345224e-6*T**2*mu**2*sigma + \
         0.00022442210764479*T**2*mu**2 + \
         7.6319402846e-6*T**2*mu*sigma**2 + \
         7.2359722222e-6*T**2*mu*sigma - \
         0.0003491716663951*T**2*mu - \
         1.1631111786e-6*T**2*sigma**4 + \
         2.6724248621e-6*T**2*sigma**3 + \
         1.649502809964e-5*T**2*sigma**2 - \
         6.1536219106916e-5*T**2*sigma + \
         0.000140143705596224*T**2 + \
         2.7736623456e-7*T*V**4 + \
         1.5654708427e-5*T*V**3 + \
         2.1280663429e-6*T*V**2*kappa**2 - \
         3.8522365278e-6*T*V**2*kappa + \
         2.1368239286e-5*T*V**2*mu**2 - \
         3.6494785e-5*T*V**2*mu + \
         2.4667129876e-8*T*V**2*sigma**2 + \
         8.2984694444e-7*T*V**2*sigma - \
         1.16151114743199e-6*T*V**2 - \
         2.4386513472e-5*T*V*kappa**2 + \
         7.2626637568e-5*T*V*kappa - \
         8.6385334444e-5*T*V*mu**2 + \
         0.00022158505878*T*V*mu - \
         6.7549275e-6*T*V*sigma**2 + \
         6.7690826574e-5*T*V*sigma - \
         0.000319030815886*T*V + \
         5.7579921563e-6*T*kappa**4 - \
         5.9128382699e-6*T*kappa**3 + \
         5.8599504703e-6*T*kappa**2*mu**2 + \
         1.0919901389e-5*T*kappa**2*mu + \
         2.3570428983e-6*T*kappa**2*sigma**2 + \
         1.1143115278e-5*T*kappa**2*sigma - \
         9.05515382865e-5*T*kappa**2 - \
         5.1905069444e-6*T*kappa*mu**2 - \
         4.8460056021e-5*T*kappa*mu - \
         1.9954263889e-6*T*kappa*sigma**2 - \
         4.0364821013e-5*T*kappa*sigma + \
         0.0002112574003288*T*kappa + \
         3.6831785874e-5*T*mu**4 - \
         2.2679927293e-5*T*mu**3 + \
         1.3165213945e-5*T*mu**2*sigma**2 - \
         1.5411466667e-5*T*mu**2*sigma - \
         0.0005305662075903*T*mu**2 - \
         1.4521058333e-5*T*mu*sigma**2 - \
         2.36356423e-5*T*mu*sigma + \
         0.0008959089906671*T*mu + \
         3.1144363024e-6*T*sigma**4 - \
         1.1900713105e-5*T*sigma**3 - \
         3.6444491206927e-5*T*sigma**2 + \
         0.000221944115643271*T*sigma - \
         0.000541037097804642*T + \
         1.3872452626e-7*V**6 + \
         2.6280928855e-6*V**5 + \
         1.4116924747e-6*V**4*kappa**2 - \
         2.6532439544e-6*V**4*kappa + \
         4.0826945223e-6*V**4*mu**2 - \
         6.7305962741e-6*V**4*mu + \
         3.2377252949e-7*V**4*sigma**2 - \
         2.3853637883e-7*V**4*sigma - \
         2.12902319506e-6*V**4 - \
         3.4146924781e-6*V**3*kappa**2 + \
         1.1308634888e-5*V**3*kappa - \
         4.3203910556e-6*V**3*mu**2 + \
         2.0522142146e-5*V**3*mu - \
         1.1130033126e-7*V**3*sigma**2 + \
         9.6641202763e-6*V**3*sigma - \
         6.9565263098929e-5*V**3 + \
         9.124821437e-7*V**2*kappa**4 - \
         2.374870717e-7*V**2*kappa**3 + \
         2.0551150926e-6*V**2*kappa**2*mu**2 + \
         2.1737964134e-6*V**2*kappa**2*mu + \
         1.5400502778e-6*V**2*kappa**2*sigma**2 + \
         2.1762389258e-6*V**2*kappa**2*sigma - \
         1.908307985662e-5*V**2*kappa**2 - \
         2.2644466603e-6*V**2*kappa*mu**2 - \
         7.2844616667e-6*V**2*kappa*mu - \
         1.327008481e-6*V**2*kappa*sigma**2 - \
         7.3647294444e-6*V**2*kappa*sigma + \
         3.083805354799e-5*V**2*kappa + \
         6.376029485e-6*V**2*mu**4 - \
         3.4468156835e-6*V**2*mu**3 + \
         1.9027688889e-6*V**2*mu**2*sigma**2 - \
         1.9749133587e-6*V**2*mu**2*sigma - \
         8.992990807087e-5*V**2*mu**2 - \
         1.2826181036e-6*V**2*mu*sigma**2 - \
         8.6458333333e-7*V**2*mu*sigma + \
         0.000119137833700857*V**2*mu + \
         3.1330206897e-7*V**2*sigma**4 - \
         5.9611039059e-7*V**2*sigma**3 - \
         5.61848663091e-6*V**2*sigma**2 + \
         1.0076975521136e-5*V**2*sigma - \
         2.21975659993502e-6*V**2 - \
         5.5542905873e-6*V*kappa**4 + \
         8.2797491074e-6*V*kappa**3 - \
         5.5680226085e-6*V*kappa**2*mu**2 - \
         1.0037508333e-6*V*kappa**2*mu - \
         5.6383374563e-6*V*kappa**2*sigma**2 + \
         4.7589241667e-6*V*kappa**2*sigma + \
         8.232677423507e-5*V*kappa**2 + \
         2.0663045e-5*V*kappa*mu**2 + \
         7.3770892156e-6*V*kappa*mu + \
         1.2777554444e-5*V*kappa*sigma**2 + \
         4.4264848543e-6*V*kappa*sigma - \
         0.0002273892174654*V*kappa - \
         1.7755825532e-5*V*mu**4 + \
         2.2503273728e-5*V*mu**3 - \
         7.4071030901e-6*V*mu**2*sigma**2 + \
         5.0431555833e-5*V*mu**2*sigma + \
         0.00023575988653791*V*mu**2 + \
         1.6386676389e-5*V*mu*sigma**2 - \
         4.6608524981e-5*V*mu*sigma - \
         0.00059482989619126*V*mu - \
         4.6406148944e-7*V*sigma**4 + \
         1.2492365373e-5*V*sigma**3 + \
         6.64618520895e-6*V*sigma**2 - \
         0.00023138236574996*V*sigma + \
         0.000737904956621858*V + \
         8.9855895986e-7*kappa**6 + \
         2.5470308117e-9*kappa**5 + \
         1.2059968463e-6*kappa**4*mu**2 + \
         3.1793778702e-6*kappa**4*mu + \
         1.0401467472e-6*kappa**4*sigma**2 + \
         2.1906938947e-6*kappa**4*sigma - \
         2.433005812556e-5*kappa**4 - \
         4.7304299279e-7*kappa**3*mu**2 - \
         6.5083289391e-6*kappa**3*mu - \
         1.4235522045e-7*kappa**3*sigma**2 - \
         5.8027556768e-6*kappa**3*sigma + \
         1.4243394357223e-5*kappa**3 + \
         3.9738486622e-6*kappa**2*mu**4 - \
         2.1223005863e-6*kappa**2*mu**3 + \
         1.0183000926e-5*kappa**2*mu**2*sigma**2 - \
         1.2959512062e-5*kappa**2*mu**2*sigma - \
         3.411711272594e-5*kappa**2*mu**2 - \
         1.0317703172e-5*kappa**2*mu*sigma**2 + \
         1.7773383333e-5*kappa**2*mu*sigma - \
         3.0801900036594e-5*kappa**2*mu + \
         1.5939164141e-6*kappa**2*sigma**4 + \
         1.3527431493e-6*kappa**2*sigma**3 - \
         2.054649657405e-5*kappa**2*sigma**2 - \
         2.925616248782e-5*kappa**2*sigma + \
         0.00020667647366024*kappa**2 - \
         7.1129482287e-6*kappa*mu**4 + \
         2.347842395e-7*kappa*mu**3 - \
         2.3796144071e-5*kappa*mu**2*sigma**2 + \
         1.8951372222e-5*kappa*mu**2*sigma + \
         4.588667312192e-5*kappa*mu**2 + \
         2.7744291667e-5*kappa*mu*sigma**2 - \
         4.024641369e-5*kappa*mu*sigma + \
         0.0001409512704825*kappa*mu - \
         2.3923579072e-6*kappa*sigma**4 - \
         6.8064892728e-6*kappa*sigma**3 + \
         2.465996737684e-5*kappa*sigma**2 + \
         0.000123292481750089*kappa*sigma - \
         0.000440486811020821*kappa + \
         4.9302956951e-6*mu**6 + \
         2.366766686e-6*mu**5 + \
         8.8283463137e-6*mu**4*sigma**2 - \
         1.3070610726e-5*mu**4*sigma - \
         0.0001299209556219*mu**4 - \
         1.7778517447e-5*mu**3*sigma**2 + \
         2.0156346188e-5*mu**3*sigma + \
         7.22110904344e-6*mu**3 + \
         2.7912531806e-7*mu**2*sigma**4 - \
         1.6060547415e-6*mu**2*sigma**3 - \
         4.639438308883e-5*mu**2*sigma**2 + \
         5.69529527941e-5*mu**2*sigma + \
         0.00107865381644979*mu**2 + \
         1.6967924396e-6*mu*sigma**4 - \
         8.6474054636e-6*mu*sigma**3 + \
         5.0242308290501e-5*mu*sigma**2 + \
         0.00011747955359353*mu*sigma - \
         0.00160943569318205*mu + \
         7.81942578e-7*sigma**6 - \
         2.6973326406e-6*sigma**5 - \
         1.688778013916e-5*sigma**4 + \
         6.520497638123e-5*sigma**3 + \
         9.1407247080762e-5*sigma**2 - \
         0.00056464306868082*sigma + \
         0.00111457799998979

    return Smax
