"""
.. module:: parcel
    :synopsis: Collection of droplet activation routines

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""

import numpy as numpy
from scipy.special import erfc erf

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

def ming2006(V, T, P, aerosol):
    """Ming activation scheme.

    NOTE - right now, the variable names correspond to the FORTRAN implementation
    of the routine. Will change in the future.

    TODO: rename variables
    TODO: docstring
    TODO: extend for multiple modes.
    """

    Num = aerosol.Nis*1e-6

    RpDry = aerosol.mu*1e-6
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
    G_a = (rho_w*R*T)/(es(T-273.15)*Dv_T(T)*Mw)
    G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka_T(T)*T)
    G = 1./(G_a + G_b)

    Smis = []
    Sparts = []
    for aerosol in aerosols:

        sigma = aerosol.sigma
        am = aerosol.mu*1e-6
        N = aerosol.N*1e6
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
        ui = 2.*np.log(Smi/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.sigma))
        N_act = 0.5*aerosol.N*erfc(ui)
        act_fracs.append(N_act/aerosol.N)

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
    Dp_B = 2.*Dv_T(T)*np.sqrt(2*np.pi*Mw/R/T)/ac
    Dp_diff = Dp_big - Dp_low
    Dv_ave = (Dv_T(T)/Dp_diff)*(Dp_diff - Dp_B*np.log((Dp_big + Dp_B)/(Dp_low+Dp_B)))

    G_a = (rho_w*R*T)/(es(T-273.15)*Dv_ave*Mw)
    G_b = (L*rho_w)*(L*Mw/R/T - 1.0)/(ka_T(T)*T)
    G = 4./(G_a + G_b)

    alpha = (g*Mw*L)/(Cp*R*(T**2)) - (g*Ma)/(R*T)
    gamma = (P*Ma)/(Mw*es(T-273.15)) + (Mw*L*L)/(Cp*R*T*T)

    ## Compute sgi of each mode
    sgis = []
    for aerosol in aerosols:
        _, sgi = kohler_crit(T, aerosol.mu*1e-6, aerosol.kappa)
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

            log_sig = np.log(aerosol.sigma)
            Ni = aerosol.N*1e6

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

            beta = 0.5*np.pi*gamma*rho_w*G/alpha/V#/rho_air
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
        ui = 2.*np.log(sgi/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.sigma))
        N_act = 0.5*aerosol.N*erfc(ui)
        act_fracs.append(N_act/aerosol.N)

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
        _, sgi = kohler_crit(T, aerosol.mu*1e-6, aerosol.kappa)
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
            tpi = aerosol.N*1e6

            dlgsg = np.log(aerosol.sigma)
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
        ui = 2.*np.log(sgi/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.sigma))
        N_act = 0.5*aerosol.N*erfc(ui)
        act_fracs.append(N_act/aerosol.N)

    return Smax, act_fracs