from parcel_model.micro import *
import numpy as np

if __name__ == "__main__":

    from parcel_model.parcel import AerosolSpecies
    from parcel_model.lognorm import Lognorm
    from parcel_model.micro import activate_ARG, activate_FN2

    import sys

    P0 = 80000. # Pressure, Pa
    T0 = 283.15 # Temperature, K
    S0 = -0.00 # Supersaturation. 1-RH from wv term
    #V = 0.1 # m/s
    V = float(sys.argv[1])

    aerosol1 = AerosolSpecies('(NH4)2SO4', Lognorm(mu=0.05, sigma=2., N=1000.),
                          bins=200, kappa=0.71)#, r_min=0.001, r_max=5.0)
    initial_aerosols = [aerosol1, ]
    #ming_Smax, ming_ratio = activate_ming(V, T0, P0, initial_aerosols[0])

    V = V
    T = T0
    P = P0
    aerosol = initial_aerosols[0]

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
        G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka(T, 1.3e3, (Dpc*1e-6)/2.)*T)
        G = 1./(G_a + G_b) # multiply by four since we're doing diameter this time

        Smax_large = (Smax > Sc) # if(Smax>Sc(count1,count2))
        WetDp = np.zeros_like(Dpc)
        #WetDp[Smax_large] = np.sqrt(Dpc[Smax_large]**2. + 1e12*(G[Smax_large]/(alpha*V))*((Smax-.0)**2.4 - (Sc[Smax_large]-.0)**2.4))
        WetDp[Smax_large] = 1e6*np.sqrt((Dpc[Smax_large]*1e-6)**2. + (G[Smax_large]/(alpha*V))*((Smax-1.0)**2.5 - (Sc[Smax_large]-1.0)**2.5))

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

        print iter_count, "%r %r %r" % (Smax, CondenRate, alpha*V/gamma)
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
        #if iter_count == 9: break

    print "---"*20

    print "    FN2", activate_FN2(V, T, P, initial_aerosols)[0]
    print "   MING", (Smax-1.0)
    print "    ARG", activate_ARG(V, T, P, initial_aerosols)[0]

    parcel = False
    if parcel:
        from parcel_model.parcel import ParcelModel
        pm = ParcelModel(initial_aerosols, V, T0, S0, P0, console=False)
        parcel, aerosols = pm.run(100.0, 0.001)
        print " PARCEL", parcel.S.max()

