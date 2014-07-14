"""
.. module:: parcel
    :synopsis: Collection of droplet activation routines

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""

import numpy as np
from numpy import min as nmin
from scipy.special import erfc, erf, erfinv

from thermo import *
from constants import *

def act_fraction(Smax, T, rs, kappa, r_drys, Nis):
    """Calculates the equilibrium activated fraction given the details of a population
    of aerosol sizes.

    NOTE - This works for a *single mode*. In order to study the entire aerosol
    population activated in the parcel model, this will need to be called however
    many times there are modes for each separate mode.

    """
    r_crits, s_crits = zip(*[kohler_crit(T, r_dry, kappa) for r_dry in r_drys])
    s_crits = np.array(s_crits)
    r_crits = np.array(r_crits)

    activated_eq = (Smax >= s_crits)
    activated_kn = (rs >= r_crits)

    N_tot = np.sum(Nis)

    eq_frac = np.sum(Nis[activated_eq])/N_tot
    kn_frac = np.sum(Nis[activated_kn])/N_tot

    return eq_frac, kn_frac

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

def arg2000(V, T, P, aerosols):

    ## Originally from Abdul-Razzak 1998 w/ Ma. Need kappa formulation
    alpha = (g*Mw*L)/(Cp*R*(T**2)) - (g*Ma)/(R*T)
    gamma = (R*T)/(es(T-273.15)*Mw) + (Mw*(L**2))/(Cp*Ma*T*P)

    ## Condensation effects
    G_a = (rho_w*R*T)/(es(T-273.15)*dv_cont(T, P)*Mw)
    G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka_cont(T)*T)
    G = 1./(G_a + G_b)

    Smis = []
    Sparts = []
    for aerosol in aerosols:

        sigma = aerosol.distribution.sigma
        am = aerosol.distribution.mu*1e-6
        N = aerosol.distribution.N*1e6
        kappa = aerosol.kappa

        fi = 0.5*np.exp(2.5*(np.log(sigma)**2))
        gi = 1.0 + 0.25*np.log(sigma)


        A = (2.*sigma_w(T)*Mw)/(rho_w*R*T)
        _, Smi2 = kohler_crit(T, am, kappa)

        zeta = (2./3.)*A*(np.sqrt(alpha*V/G))
        etai = ((alpha*V/G)**(3./2.))/(N*gamma*rho_w*2.*np.pi)

        ##
        Spa = fi*((zeta/etai)**(1.5))
        Spb = gi*(((Smi2**2)/(etai + 3.*zeta))**(0.75))
        S_part = (1./(Smi2**2))*(Spa + Spb)

        Smis.append(Smi2)
        Sparts.append(S_part)

    Smax = 1./np.sqrt(np.sum(Sparts))

    act_fracs = []
    for Smi, aerosol in zip(Smis, aerosols):
        ui = 2.*np.log(Smi/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.distribution.sigma))
        N_act = 0.5*aerosol.distribution.N*erfc(ui)
        act_fracs.append(N_act/aerosol.distribution.N)

    return Smax, act_fracs

def fn2005(V, T, P, aerosols, tol=1e-6, max_iters=100):
    """
    NS 2003 algorithm + FN 2005 corrections for diffusive growth rate
    """
    #aer = aerosols[0]

    A = (4.*Mw*sigma_w(T))/(R*T*rho_w)
    ## Compute rho_air by assuming air is saturated at given T, P
    # Petty (3.41)
    qsat = 0.622*(es(T-273.15)/P)
    Tv = T*(1.0 + 0.61*qsat)
    rho_air = P/Rd/Tv # air density
    #print "rho_air", rho_air

    Dp_big = 5e-6
    Dp_low = np.min([0.207683*(ac**-0.33048), 5.0])*1e-5
    Dp_B = 2.*dv_cont(T, P)*np.sqrt(2*np.pi*Mw/R/T)/ac
    Dp_diff = Dp_big - Dp_low
    Dv_ave = (dv_cont(T, P)/Dp_diff)*(Dp_diff - Dp_B*np.log((Dp_big + Dp_B)/(Dp_low+Dp_B)))

    G_a = (rho_w*R*T)/(es(T-273.15)*Dv_ave*Mw)
    G_b = (L*rho_w)*(L*Mw/R/T - 1.0)/(ka_cont(T)*T)
    G = 4./(G_a + G_b)

    alpha = (g*Mw*L)/(Cp*R*(T**2)) - (g*Ma)/(R*T)
    gamma = (P*Ma)/(Mw*es(T-273.15)) + (Mw*L*L)/(Cp*R*T*T)

    ## Compute sgi of each mode
    sgis = []
    for aerosol in aerosols:
        _, sgi = kohler_crit(T, aerosol.distribution.mu*1e-6, aerosol.kappa)
        sgis.append(sgi)
    #print "--"*20

    def S_integral(Smax):
        delta = Smax**4 - (16./9.)*(A*A*alpha*V/G)
        #delta = 1.0 - (16./9.)*alpha*V*(1./G)*((A/(Smax**2))**2)

        if delta > 0:
            sp_sm_sq = 0.5*(1.0 + np.sqrt(delta))
            sp_sm = np.sqrt(sp_sm_sq)
        else:
            arg = (2e7*A/3.)*Smax**(-0.3824)
            sp_sm = np.min([arg, 1.0])

        Spart = sp_sm*Smax
        #print "Spart", Smax, Spart, delta

        I1s, I2s = 0.0, 0.0
        for aerosol, sgi in zip(aerosols, sgis):

            log_sig = np.log(aerosol.distribution.sigma)
            Ni = aerosol.distribution.N*1e6

            upart = 2.*np.log(sgi/Spart)/(3.*np.sqrt(2)*log_sig)
            umax = 2.*np.log(sgi/Smax)/(3.*np.sqrt(2)*log_sig)

            def I1(Smax):
                A1 = (Ni/2.)*((G/alpha/V)**0.5)
                A2 = Smax
                C1 = erfc(upart)
                C2 = 0.5*((sgi/Smax)**2)
                C3a = np.exp(9.*(log_sig**2)/2.)
                C3b = erfc(upart + 3.*log_sig/np.sqrt(2.))
                return A1*A2*(C1 - C2*C3a*C3b)

            def I2(Smax):
                A1 = A*Ni/3./sgi
                A2 = np.exp(9.*(log_sig**2.)/8.)
                C1 = erf(upart - 3.*log_sig/(2.*np.sqrt(2.)))
                C2 = erf(umax - 3.*log_sig/(2.*np.sqrt(2.)))
                return A1*A2*(C1 - C2)

            beta = 0.5*np.pi*gamma*rho_w*G/alpha/V/rho_air
            #beta = 0.5*np.pi*gamma*rho_w/bet2_par/alpha/V/rho_air
            #print "++", Smax, I1(Smax), I2(Smax)

            I1s += I1(Smax)
            I2s += I2(Smax)

        return Smax*beta*(I1s + I2s) - 1.0

    x1 = 1e-5
    y1 = S_integral(x1)
    x2 = 1.0
    y2 = S_integral(x2)

    #print (x1, y1), (x2, y2)

    #print "BISECTION"
    #print "--"*20
    for i in xrange(max_iters):

        ## Bisection
        #y1, y2 = S_integral(x1), S_integral(x2)
        x3 = 0.5*(x1+x2)
        y3 = S_integral(x3)
        #print "--", x3, y3, "--"

        if np.sign(y1)*np.sign(y3) <= 0.:
            x2 = x3
            y2 = y3
        else:
            x1 = x3
            y1 = y3

        if np.abs(x2-x1) < tol*x1: break

        #print i, (x1, y1), (x2, y2)

    ## Converged ; return
    x3 = 0.5*(x1 + x2)

    Smax = x3
    #print "Smax = %f (%f)" % (Smax, 0.0)

    act_fracs = []
    for aerosol, sgi in zip(aerosols, sgis):
        ui = 2.*np.log(sgi/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.distribution.sigma))
        N_act = 0.5*aerosol.distribution.N*erfc(ui)
        act_fracs.append(N_act/aerosol.distribution.N)

    return Smax, act_fracs


def ns2003(V, T, P, aerosols, tol=1e-6, max_iters=100):
    """Sketch implementation of Nenes and Seinfeld (2003) parameterization
    """

    nmd_par = len(aerosols) # number of modes
    vhfi = 3.0 # van't hoff factor (ions/molecule)

    ## Setup constants
    akoh_par = 4.0*Mw*sigma_w(T)/R/T/rho_w
    ## Compute rho_air by assuming air is saturated at given T, P
    # Petty (3.41)
    qsat = 0.622*(es(T-273.15)/P)
    Tv = T*(1.0 + 0.61*qsat)
    rho_air = P/Rd/Tv # air density

    alpha = g*Mw*L/Cp/R/T/T - g*Ma/R/T
    gamma = P*Ma/es(T-273.15)/Mw + Mw*L*L/Cp/R/T/T

    bet2_par = R*T*rho_w/es(T-273.15)/Dv/Mw/4. + L*rho_w/4./Ka/T*(L*Mw/R/T - 1.0)
    beta = 0.5*np.pi*gamma*rho_w/bet2_par/alpha/V/rho_air

    cf1  = 0.5*(((1/bet2_par)/(alpha*V))**0.5)
    cf2  = akoh_par/3.0

    sgis = []
    for aerosol in aerosols:
        _, sgi = kohler_crit(T, aerosol.distribution.mu*1e-6, aerosol.kappa)
        sgis.append(sgi)

    def sintegral(spar):
        ## descriminant criterion
        descr = 1.0 - (16./9.)*alpha*V*bet2_par*((akoh_par/(spar**2))**2)

        if descr <= 0.0:
            crit2 = True
            ratio = (2e7/3.)*akoh_par*spar**(-0.3824)
            if ratio > 1.0: ratio = 1.0
            ssplt2 = spar*ratio
        else:
            crit2 = False
            ssplt1 = 0.5*(1.0 - np.sqrt(descr)) # min root of both
            ssplt2 = 0.5*(1.0 + np.sqrt(descr)) # max root of both
            ssplt1 = np.sqrt(ssplt1)*spar # multiply ratios with smax
            ssplt2 = np.sqrt(ssplt2)*spar
        ssplt = ssplt2 # store ssplit in common

        summ, summat, summa = 0, 0, 0

        sqtwo = np.sqrt(2.0)

        for aerosol, sgi in zip(aerosols, sgis):

            sg_par = sgi
            tpi = aerosol.distribution.N*1e6

            dlgsg = np.log(aerosol.distribution.sigma)
            dlgsp = np.log(sg_par/spar)

            orism1 = 2.0*np.log(sg_par/ssplt2)/(3.*sqtwo*dlgsg)
            orism2 = orism1 - 3.0*dlgsg/(2.0*sqtwo)

            orism3 = 2.0*dlgsp/(3.0*sqtwo*dlgsg) - 3.0*dlgsg/(2.0*sqtwo)
            orism4 = orism1 + 3.0*dlgsg/sqtwo

            orism5 = 2.0*dlgsp/(3*sqtwo*dlgsg)
            ekth = np.exp((9./2.)*dlgsg*dlgsg)

            integ1 = tpi*spar*((1-erf(orism1)) - 0.5*((sg_par/spar)**2)*ekth*(1.0-erf(orism4)))
            integ2 = (np.exp((9./8.)*dlgsg*dlgsg)*tpi/sg_par)*(erf(orism2) - erf(orism3))

            nd = (tpi/2.)*(1.0 - erf(orism5))

            summ += integ1
            summat += integ2
            summa += nd

        return summa, summ, summat

    ## Initial values for bisection
    x1 = 1e-5 # min cloud supersaturation
    ndrpl, sinteg1, sinteg2 = sintegral(x1)
    print ndrpl, sinteg1, sinteg2
    y1 = (sinteg1*cf1 + sinteg2*cf2)*beta*x1 - 1.0

    x2 = 1.0 # max cloud supersaturation
    ndrpl, sinteg1, sinteg2 = sintegral(x2)
    print ndrpl, sinteg1, sinteg2
    y2 = (sinteg1*cf1 + sinteg2*cf2)*beta*x2 - 1.0

    print (x1, y1), (x2, y2)

    print "BISECTION"
    print "--"*20
    ## Perform bisection
    for i in xrange(max_iters):
        x3 = 0.5*(x1 + x2)
        ndrpl, sinteg1, sinteg3 = sintegral(x3)
        y3 = (sinteg1*cf1 + sinteg2*cf2)*beta*x3 - 1.0

        if np.sign(y1)*np.sign(y3) <= 0.:
            y2 = y3
            x2 = x3
        else:
            y1 = y3
            x1 = x3

        if np.abs(x2-x1) <= tol*x1: break

        print i, (x1, y1), (x2, y2)

    ## Converged ; return
    x3 = 0.5*(x1 + x2)
    ndrpl, sinteg1, sinteg3 = sintegral(x3)
    y3 = (sinteg1*cf1 + sinteg2*cf2)*beta*x3 - 1.0

    Smax = x3
    print "Smax = %f (%f)" % (Smax, 0.0)

    act_fracs = []
    for aerosol, sgi in zip(aerosols, sgis):
        ui = 2.*np.log(sgi/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.distribution.sigma))
        N_act = 0.5*aerosol.distribution.N*erfc(ui)
        act_fracs.append(N_act/aerosol.distribution.N)

    return Smax, act_fracs

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
