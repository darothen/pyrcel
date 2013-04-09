from parcel_model.micro import *
import numpy as np

def activate_ming(V, T, P, aerosol):

    Num = aerosol.N

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
    Dpc = np.array(r_crits)
    Dp0 = r_crits/np.sqrt(3.)
    Sc = np.array(s_crits)+1.0
    DryDp = aerosol.r_drys*2.

    ## Begin algorithm
    Smax1 = 1.0
    Smax2 = 1.1

    iter_count = 1
    while (Smax2 - Smax1) > 1e-7:
        print iter_count, Smax2, Smax1
        Smax = 0.5*(Smax2 + Smax1)

        ## subroutine Grow()

        ## subroutine CalcG()
        # TODO: implement size-dependent effects on Dv, ka, using Dpc
        G_a = (rho_w*R*T)/(es(T-273.15)*Dv_T(T)*Mw)
        G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka_T(T)*T)
        G = 1./(G_a + G_b)

        Smax_large = (Smax > Sc) # if(Smax>Sc(count1,count2))
        WetDp = Dpc[Smax_large]**2. + 1e12*(G/(alpha*V))*np.sqrt((Smax-1.0)**2.5 - (Sc[Smax_large]-1.0)**2.5)

        ## subroutine Activity()
        # Is this just the Kohler curve?
        Act = np.ones_like(WetDp)
        WetDp_large = (WetDp > 1e-5) # if(WetDp(i,j)>1e-5)
        Act[WetDp_large] = Seq(WetDp[WetDp_large]/2., DryDp[WetDp_large]/2., T, kappa) + 1.0

        ## subroutine Conden()

        ## subroutine CalcG()
        # TODO: implement size-dependent effects on Dv, ka, using WetDp
        # note - before the above TODO is implemented, this is a redundant calc and should be skipped
        #G_a = (rho_w*R*T)/(es(T-273.15)*Dv_T(T)*Mw)
        #G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka_T(T)*T)
        #G = 1./(G_a + G_b)

        WetDp_large = (WetDp > Dpc) # (WetDp(count1,count2)>Dpc(count1,count2))
        CondenRate = np.sum((np.pi/2.)*1e3*G*(WetDp[WetDp_large]*1e-6)*Num[WetDp_large]*1e6*
                             (Smax-Act[WetDp_large]))
        DropletNum = np.sum(Num[WetDp_large])
        ActDp = 0.0
        for i in xrange(1, len(WetDp)):
            if (WetDp[i] > Dpc[i]) and (WetDp[i-1] < Dpc[i]):
                ActDp = DryDp[i]

        ## Iteration logic
        if CondenRate < (alpha*V/gamma):
            Smax1 = Smax
        else:
            Smax2 = Smax

        iter_count += 1

    return Smax, None

if __name__ == "__main__":

    from parcel_model.parcel import AerosolSpecies
    from parcel_model.lognorm import Lognorm
    from parcel_model.micro import activate_ARG, activate_FN2

    P0 = 80000. # Pressure, Pa
    T0 = 283.15 # Temperature, K
    S0 = -0.00 # Supersaturation. 1-RH from wv term
    V = 0.5 # m/s

    aerosol1 = AerosolSpecies('(NH4)2SO4', Lognorm(mu=0.01, sigma=2.5, N=200.),
                          bins=40, kappa=0.71, r_min=0.001, r_max=5.0)
    initial_aerosols = [aerosol1, ]
    #ming_Smax, ming_ratio = activate_ming(V, T0, P0, initial_aerosols[0])

    V = V
    T = T0
    P = P0
    aerosol = initial_aerosols[0]

    Num = aerosol.Nis

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
    while (Smax2 - Smax1) > 1e-5:
        print "\t", iter_count, Smax1, Smax2
        Smax = 0.5*(Smax2 + Smax1)
        print "---", Smax-1.0

        ## subroutine Grow()

        ## subroutine CalcG()
        # TODO: implement size-dependent effects on Dv, ka, using Dpc
        #G_a = (rho_w*R*T)/(es(T-273.15)*Dv_T(T)*Mw)
        G_a = (rho_w*R*T)/(es(T-273.15)*dv(T, Dpc/2.)*Mw)
        #G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka_T(T)*T)
        G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka(T, 1.007e3, Dpc/2.)*T)
        G = 4./(G_a + G_b) # multiply by four since we're doing diameter this time

        Smax_large = (Smax > Sc) # if(Smax>Sc(count1,count2))
        WetDp = np.zeros_like(Dpc)
        WetDp[Smax_large] = np.sqrt(Dpc[Smax_large]**2. + 1e12*(G[Smax_large]/(alpha*V))*((Smax-1.0)**2.5 - (Sc[Smax_large]-1.0)**2.5))

        #print Dpc
        #print WetDp/DryDp
        #print WetDp

        ## subroutine Activity()
        # Is this just the Kohler curve?
        Act = np.ones_like(WetDp)
        WetDp_large = (WetDp > 1e-5) # if(WetDp(i,j)>1e-5)
        Act[WetDp_large] = Seq((WetDp[WetDp_large]*1e-6)/2., DryDp[WetDp_large]/2., T, kappa) + 1.0

        #print Act

        ## subroutine Conden()

        ## subroutine CalcG()
        # TODO: implement size-dependent effects on Dv, ka, using WetDp
        #G_a = (rho_w*R*T)/(es(T-273.15)*Dv_T(T)*Mw)
        G_a = (rho_w*R*T)/(es(T-273.15)*dv(T, WetDp/2.)*Mw)
        #G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka_T(T)*T)
        G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka(T, 1.007e3, WetDp/2.)*T)
        G = 4./(G_a + G_b) # multiply by four since we're doing diameter this time

        WetDp_large = (WetDp > Dpc) # (WetDp(count1,count2)>Dpc(count1,count2))
        #WetDp_large = (WetDp > 0)
        print (Smax-Act[WetDp_large])
        CondenRate = np.sum((np.pi/2.)*1e3*G[WetDp_large]*(WetDp[WetDp_large]*1e-6)*Num[WetDp_large]*1e6*
                             (Smax-Act[WetDp_large]))

        print "%r" % (CondenRate, )
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

    print "---"*20
    print "MING", Smax-1.0
    print " ARG", activate_ARG(V, T, P, initial_aerosols)[0]
    print " FN2", activate_FN2(V, T, P, initial_aerosols)[0]
