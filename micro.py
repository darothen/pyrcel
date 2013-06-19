"""
.. module:: parcel
    :synopsis: Useful and commonly used microphysics constants and functions.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""
import numpy as np
from scipy.optimize import fminbound, bisect
from scipy.special import erfc, erf

## Microphysics constants
g = 9.81             #: Gravitational constant, m/s**2
Cp = 1004.0          #: Specific heat of dry air at constant pressure, J/kg
L = 2.5e6            #: Latent heat of condensation, J/kg
rho_w = 1e3          #: Density of water, kg/m**3
Rd = 287.0           #: Gas constant for dry air, J/(kg K)
R = 8.314            #: Universal gas constant, J/(mol K)
Mw = 18.0153/1e3     #: Molecular weight of water, kg/mol

# Additional constants from the auxiliary file
Ma = 28.9/1e3        #: Molecular weight of dry air, kg/mol
#Dv = 3.e-5          #: Diffusivity of water vapor in air, m^2/s
ac = 1.0             #: condensation constant
Ka = 2.e-2           #: Thermal conductivity of air, J/m/s/K
at = 0.96            #: thermal accomodation coefficient
epsilon = 0.622      #: molecular weight of water / molecular weight of dry air

## NOT CORRECTING FOR NON-CONTINUUM EFFECTS
Dv = 0.3/1e4 # Diffusivity of water vapor in air, m^2/s
Dv_T = lambda T, P=1.0 : 1e-4*(0.211/P)*((T/273.)**1.94) # Diffusivity of water vapor in air, m^2/s. T is Kelvin; P is atm, assumed 1
ka_T = lambda T: 1e-3*(4.39 + 0.071*T) # thermal conductivty of air, W/(m K) given T in Kelvin


## AUXILIARY FUNCTIONS
def dv(T, r):
    """Diffusivity of water vapor in air, modified for non-continuum effects

    Revise with equation 17.62, Seinfeld and Pandis?"""
    denom = 1.0 + (Dv_T(T)/(ac*r))*np.sqrt((2*np.pi*Mw)/(R*T))
    return Dv_T(T)/denom

def ka(T, rho, r):
    """Thermal conductivity of air, modified for non-continuum effects

    Revise with equation 17.71, Seinfeld and Pandis?"""
    denom = 1.0 + (ka_T(T)/(at*r*rho*Cp))*np.sqrt((2*np.pi*Ma)/(R*T))
    return ka_T(T)/denom

def activate_ming(V, T, P, aerosol):
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

def activate_new(V, T, P, aerosols, max_iters=30):

    ### Constants/coefficients
    A = (4.*Mw*sigma_w(T))/(R*T*rho_w)
    qsat = 0.622*(es(T-273.15)/P)
    Tv = T*(1.0 + 0.61*qsat)
    rho_air = P/Rd/Tv # air density

    #Dp_big = 5e-6
    #Dp_low = np.min([0.207683*(ac**-0.33048), 5.0])*1e-5
    #Dp_B = 2.*Dv_T(T)*np.sqrt(2*np.pi*Mw/R/T)/ac
    #Dp_diff = Dp_big - Dp_low
    #Dv_ave = (Dv_T(T)/Dp_diff)*(Dp_diff - Dp_B*np.log((Dp_big + Dp_B)/(Dp_low+Dp_B)))
    Dv_ave = Dv_T(T)

    G_a = (rho_w*R*T)/(es(T-273.15)*Dv_ave*Mw)
    G_b = (L*rho_w)*(L*Mw/R/T - 1.0)/(ka_T(T)*T)
    G = 1./(G_a + G_b)

    alpha = (g*Mw*L)/(Cp*R*(T**2)) - (g*Ma)/(R*T)
    gamma = (P*Ma)/(Mw*es(T-273.15)) + (Mw*L*L)/(Cp*R*T*T)
    gamma_star = 4.*np.pi*rho_w*gamma/rho_air

    root2 = np.sqrt(2.)

    def LHS(Smax):
        return alpha*V/gamma_star/G/Smax

    def RHS(T, Smax, aerosols):
        ## Compute sgi of each mode
        sgis = []
        for aerosol in aerosols:
            _, sgi = kohler_crit(T, aerosol.mu*1e-6, aerosol.kappa)
            sgis.append(sgi)

        ## Compute integral for each mode:
        integrals = []
        for sgi, aerosol in zip(sgis, aerosols):
            kappa = aerosol.kappa
            log_sigma = np.log(aerosol.sigma)
            Ni = aerosol.N*1e6

            u = lambda Smax, sgi=sgi, log_sigma=log_sigma: (root2/3./log_sigma)*np.log(sgi/Smax)

            u_smax = u(Smax)

            ## Integral 1 -
            I1_coeff = (Ni/2.)*np.sqrt(G/alpha/V)*Smax
            I1_left = erfc(u_smax)
            I1_right = 0.5*((sgi/Smax)**2)*np.exp((9./8.)*(log_sigma**2))* \
                       erfc(u_smax + (3./2./root2)*log_sigma)
            I1 = I1_coeff*(I1_left - I1_right)

            #print I1, (I1_coeff, I1_left, I1_right)

            ## Integral 2 -
            I2_coeff = (Ni*A/3./sgi)*np.exp((9./8.)*(log_sigma**2))
            I2_arg = (1.0 - erf(u_smax - (3./2./root2)*log_sigma))
            I2 = I2_coeff*I2_arg

            #print I2, (I2_coeff, I2_arg)

            integrals.append(I1 + I2)

        return np.sum(integrals)

    l = LHS(Smax)
    r = RHS(T, Smax, aerosols)

    #return LHS(Smax), RHS(T, Smax, aerosols), l-r

    f = lambda x: LHS(x) - RHS(T, x, aerosols)

    x1 = 1e-5
    y1 = f(x1)

    x2 = 0.5
    y2 = f(x2)

    print (x1, y1), (x2, y2)
    print "BISECTION"
    print "--"*20
    for i in xrange(max_iters):

        x3 = 0.5*(x1 + x2)
        y3 = f(x3)
        if np.sign(y1)*np.sign(y3) <= 0.:
            x2 = x3
            y2 = y3
        else:
            x1 = x3
            y1 = y3

        if np.abs(x2-x1) < 1e-6*x1: break

        print i, (x1, y1), (x2, y2)

    x3 = 0.5*(x1 + x2)
    Smax = x3

    act_fracs = []
    for aerosol, sgi in zip(aerosols, sgis):
        ui = 2.*np.log(sgi/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.sigma))
        N_act = 0.5*aerosol.N*erfc(ui)
        act_fracs.append(N_act/aerosol.N)

    return Smax, act_fracs

def activate_ARG(V, T, P, aerosols):

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

def activate_FN2(V, T, P, aerosols, tol=1e-6, max_iters=100):
    """
    NS 2003 algorithm + FN 2005 coorections for diffusive growth rate
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


def activate_FN(V, T, P, aerosols, tol=1e-6, max_iters=100):
    """Sketch implementation of Fountoukis and Nenes (2005) parameterization
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

def sigma_w(T):
    """Calculate the surface tension of water for a given temperature.

    .. math::
        \sigma_w = 0.0761 - 1.55\\times 10^{-4}(T - 273.15)

    :param T: temperature (deg C)
    :type T: float
    :returns: :math:`\sigma_w(T)` (J/m**2)

    """
    return 0.0761 - (1.55e-4)*(T-273.15)

def es(T):
    """Calculates the saturation vapor pressure over water for a given temperature.

    Implements Formula 2.17 from Rogers and Yau,

    .. math::

        e_s(T) = 6.112\exp\\frac{17.67T}{T+243.5}

    where *T* is the temperature in degrees Celsius.

    :param T: temperature (deg C)
    :type T: float

    :returns: :math:`e_s(T)` (Pa)

    """
    return 611.2*np.exp(17.67*T/(T+243.5))

def Seq(r, r_dry, T, kappa, neg=False):
    """Calculates the equilibrium supersaturation (relative to 100% RH) over an
    aerosol parcel of given dry/wet radius and of specified hygroscopicity in a
    parcel of a specific temperature.

    Following the parameterization by Petters and Kredenweis (2007), classical
    Kohler theory can be modified to parameterize the hygroscopicity of an aerosol
    particle using a parameter, :math:`\kappa`. The modified theory predicts that
    the supersaturation with respect to a given aerosol particle is,

    .. math::
        S_\\text{eq} &= a_w \exp \\left( \\frac{2\sigma_{w} M_w}{RT\\rho_w r} \\right) \\\\
        a_w &= \\left(1 + \kappa\\left(\\frac{r_d}{r}^3\\right) \\right)^{-1}

    with the relevant thermodynamic properties of water defined elsewhere in this
    module, :math:`r_d` is the particle dry radius (``r_dry``), :math:`r` is the
    radius of the droplet containing the particle (``r``), :math:`T` is the temperature
    of the environment (``T``), and :math:`\kappa` is the hygroscopicity parameter
    of the particle (``kappa``).

    This method has been extended to supply the *negative* of the supersaturation if
    specified using the argument ``neg``; this is useful when attempting to numerically
    estimate the particle's critical radius, as done in :func:`kohler_crit`. Otherwise,
    this method will return the supersaturation as a decimal with respect to 1.0,

    .. math::
        S_\\text{eq} = S - 1.0

    :param r: droplet radius (m)
    :type r: float

    :param r_dry: particle dry radius (m)
    :type r_dry: float

    :param T: environment temperature (K)
    :type T: float

    :param kappa: particle hygroscopicity
    :type kappa: float > 0.0

    :param neg: flag indicating to return the negative of the supersaturation, default ``False``
    :type neg: boolean -- optional

    :returns: :math:`S_\\text{eq}`

    """
    A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    if kappa == 0.0:
        B = 1.
    else:
        B = (r**3 - (r_dry**3))/(r**3 - (r_dry**3)*(1.-kappa))
    #print A, B, np.exp(A), np.exp(A)*B
    if neg:
        return 1.0 - np.exp(A)*B
    else:
        return np.exp(A)*B - 1.0

def critical_curve(T, r_a, r_b, kappa):
    """Calculates curves of critical wet radius and supersaturations for aerosol
    particles.

    Calls :func:`kohler_crit` for values of ``r_dry`` between ``r_a`` and ``r_b``
    to calculate how the critical supersaturation changes with the dry radius for a
    particle of specified :math:`\\kappa = ` ``kappa``.

    :param T: environment temperature (K)
    :type T: float

    :param r_a, r_b: parcel dry radii left and right bounds (m)
    :type r_a, r_b: float, float

    :param kappa: particle hygroscopicity
    :type kappa: float > 0.0

    :returns: arrays of dry radii, critical radii, and supersaturations

    """
    rs = np.logspace(np.log10(r_a), np.log10(r_b), 200)
    ss = np.array(map(lambda r_dry: kohler_crit(T, r_dry, kappa), rs))
    return rs, ss[:, 0], ss[:, 1]

def kohler_crit(T, r_dry, kappa):
    """Calculates the critical radius and supersaturation of an aerosol particle.

    Passes the negative or inverted supersaturation function from :func:`Seq` to
    a minimum-value optimization routine from SciPy to efficiently calculate the
    radius at which the curve achieves its maximum supersaturation, and that
    supersaturation value.

    :param T: environment temperature (K)
    :type T: float

    :param r_dry: parcel dry radius (m)
    :type r_dry: float

    :param kappa: particle hygroscopicity
    :type kappa: float

    :returns: r_crit - critical radius (value of `r` maximizing :func:`Seq`)
    :returns: s_crit - :func:`Seq` evaluated at ``r_crit``

    """
    out = fminbound(Seq, r_dry, r_dry*1e4, args=(r_dry, T, kappa, True),
                    xtol=1e-10, full_output=True, disp=0)
    r_crit, s_crit = out[:2]
    s_crit *= -1.0
    return r_crit, s_crit

def kappa_kohler_crit(T, r_dry, kappa):
    """Calculates the critical radius and supersaturation of an aerosol particle.

    Uses the approximation Seq = A/r - kappa*(r_dry/r)**3.
    """
    A = 4.*Mw*sigma_w(T)/R/T/rho_w
    rd3 = r_dry**3
    r_crit = np.sqrt(3.*kappa*rd3/A)
    s_crit = np.sqrt(4.*(A**3)/(27.*kappa*rd3))
    return r_crit, s_crit

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

def r_eff(rho, wc, Ni):
    """Calculates the cloud droplet effective radius given the parcel liquid
    water mixing ratio and the number of activated droplets, as well as the parcel
    air density

    Assuming the droplet population is monodisperse or close to it, the cloud droplet
    effective radius can be computed by

    .. math::
        r_{\\text{eff}} = \left(\\frac{3 \\rho_a w_c}{4 \pi N_i \\rho_w}\\right)^{1/3}

    :param rho: parcel air density (kg/m^3)
    :type rho: float

    :param wc: liquid water mixing ratio (kg/kg)
    :type wc: float

    :param Ni: droplet number concentration (1/m^3)
    :type Ni: float

    :returns: droplet effective radius (m^3)
    """
    return (3.*rho*wc/(4.*np.pi*rho_w*Ni))**(1./3.)
