from libc.math cimport exp, sqrt, abs
import numpy as np
import parcel_model.constants as c

cimport cython

cdef extern from "math.h":
    bint isnan(double x)

cdef:
    double Mw = c.Mw        #: Molecular weight of water, kg/mol
    double Ma = c.Ma        #: Molecular weight of dry air, kg/mol
    double R = c.R          #: Universal gas constant, J/(mol K)
    double rho_w = c.rho_w  #: Density of water, kg/m**3
    double Rd = c.Rd        #: Gas constant for dry air, J/(kg K)
    double g = c.g          #: Gravitational constant, m/s**2
    double Dv = c.Dv        #: Diffusivity of water vapor in air, m^2/s
    double ac = c.ac        #: Condensation Constant
    double Ka = c.Ka        #: Thermal conductivity of air, J/m/s/K
    double at = c.at        #: thermal accomodation coefficient
    double L = c.L          #: Latent heat of condensation, J/kg
    double Cp = c.Cp        #: Specific heat of dry air at constant pressure, J/kg
    double PI = 3.14159265358979323846264338328 #: Pi, constant
    double a = 0.5          #: weighting coefficient for shift factor (Simmel and Wurzler, 2006)
    int N_TRAP = 1000        #: number of steps to use in internal trapezoid rule

cdef inline double cmax(double a, double b): return a if a >= b else b
cdef inline double cmin(double a, double b): return a if a <= b else b
    
cdef inline double sigma_w(double T) nogil:
    return 0.0761 - (1.55e-4)*(T-273.15)

cdef inline double ka(double T, double r, double rho) nogil:
    cdef double denom, ka_cont
    ka_cont = 1e-3*(4.39 + 0.071*T)
    denom = 1.0 + (ka_cont/(at*r*rho*Cp))*sqrt((2*PI*Ma)/(R*T))
    return ka_cont/denom

cdef inline double dv(double T, double r, double P) nogil:
    cdef double denom, dv_cont, P_atm
    P_atm = P*1.01325e-5 # Pa -> atm
    dv_cont = 1e-4*(0.211/P_atm)*((T/273.)**1.94)
    denom = 1.0 + (dv_cont/(ac*r))*sqrt((2*PI*Mw)/(R*T))
    return dv_cont/denom

cdef inline double es(double T):
    return 611.2*exp(17.67*T/(T+243.5))

cdef inline double m_to_r(double m, double rho):
    if m <= 0.: 
        return 0.
    else:
        return (m*0.75/(rho*PI))**(1./3.)

cdef inline double r_to_m(double r, double rho):
    if r <= 0.:
        return 0.
    else:
        return (4./3.)*PI*rho_w*(r**3)

cpdef double Seq(double r, double r_dry, double T, double kappa) nogil:
    cdef double A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    cdef double B = 1.0
    if kappa > 0.0:
        B = (r**3 - (r_dry**3))/(r**3 - (r_dry**3)*(1.-kappa))
    cdef double returnval = exp(A)*B - 1.0
    return returnval

cdef double growth_rate(double r, double r_dry, double T, double P, double S,
                        double kappa, double rho_air):
    cdef double pv_sat = es(T - 273.15)
    cdef:
        double G_a, G_b, G
        double delta_S, dr_dt
        double dv_r, ka_r, P_atm, Seq_r
        double dm_dyn, dm_eq
    
    ## Non-continuum diffusivity/thermal conductivity of air near
    ## near particle
    dv_r = dv(T, r, P)
    ka_r = ka(T, r, rho_air)

    ## Condensation coefficient
    G_a = (rho_w*R*T)/(pv_sat*dv_r*Mw)
    G_b = (L*rho_w*((L*Mw/(R*T))-1.))/(ka_r*T)
    G = 1./(G_a + G_b)

    ## Difference between ambient and particle equilibrium supersaturation
    Seq_r = Seq(r, r_dry, T, kappa)
    S = S
    Seq_r = Seq_r
    delta_S = S - Seq_r
    #print "        ->", r, r_dry, S, Seq_r

    ## Size and liquid water tendencies
    dr_dt = (G/r)*delta_S    
    dm_dyn = 4.*PI*rho_w*r*r*dr_dt
    
    ## Equilibrium value
    #dm_eq =  4.*PI*rho_w*r*G*S
    
    ## Check if the droplet is shrinking; if it is, we can't go lower than the dry 
    ## radius
    #if delta_S < 0:
    #    dr_max = (r_dry - r) # negative because it's shrinking
    #    dm_max = (4.*PI/3.)*rho_w*r*r*dr_max
    #    return cmax(dm_max, dm_dyn)

    #if delta_S >= 0:
    #    return cmin(dm_eq, dm_dyn)
    #elif delta_S < 0:
    #    return cmax(dm_eq, dm_dyn)
    
    return dm_dyn

@cython.cdivision(True)
cdef double cn_x_linear(double x, double xk, double xkp1, double Nk, double Mk, 
                       double moment=0.0, int flag=0):
    cdef double xk_hat = 0.5*(xk + xkp1)
    cdef double n0 = Nk/(xkp1 - xk)
    cdef double a = 12.*(Mk - xk_hat*Nk)/((xkp1 - xk)**3)
    cdef double xi = xk_hat - n0/a
    
    cdef double n_xk = n0 + a*(xk - xk_hat)
    cdef double n_xkp1 = n0 + a*(xkp1 - xk_hat)
    
    #if flag == 1:
    #    print "   n_xk", n_xk
    #    print " n_xkp1", n_xkp1

    cdef double x_star, a_star

    """
    There is a problem with the (xi - x_star)**2 calcution for a_star. If the two values
    are very close to each other, then we'll overflow and produce a garbage answer.

    As of 8/7/2013 I have no good solution to this problem. So the workaround for now will
    be to merely exclude the endpoints when I integrate over the linear approximations...
    The issue only ever crops up on the edges anyways, and having a high N_TRAP should reduce
    the numerical error from this.
    """
    
    if n_xk < 0:
        x_star = 3.*Mk/Nk - 2.*xkp1
        a_star = 2.*Nk/((xkp1 - x_star)**2)
        if x <= x_star : 
            return 0.
        else:
            return (x**moment)*a_star*(x - x_star)
        
    if n_xkp1 < 0:
        x_star = 3.*Mk/Nk - 2.*xk
        a_star = -2.*Nk/((xk - x_star)**2)
        if x >= x_star :
            return 0.
        else:
            if flag == 1: 
                print "HERE", x, x_star, x-x_star
                print "    ", x**moment, a_star, xk, x_star, x
            return (x**moment)*a_star*(x - x_star)
        
    return (x**moment)*(n0 + a*(x - xk_hat))

cdef inline double moment_linear(xl, xr, l, n_x):
    return 0.5*(xr - xl)*((xl**l)*n_x(xl) + (xr**l)*n_x(xr))

@cython.cdivision(True)
cdef double trap_integrate(double xleft, double xright, # bounds to integrate
                           double x_low, double x_high, # params for new linear func
                           double Nk, double Mk, int moment):
    cdef double Tot_add, dx
    cdef double xi
    cdef unsigned int i
    cdef double mom = moment*1.0
    
    Tot_add = 0.0
    #print xright, xleft

    dx = (xright - xleft)/(1.0*N_TRAP)
    ## Trapezoidal rule integration of linear function, N = N_TRAP
    """
    Workaround 8/7/2013: Don't add the edges to the integration! see note in cn_x_linear
    """
    #Tot_add = 0.5*cn_x_linear(xleft, x_low, x_high, Nk, Mk, mom) + \
    #          0.5*cn_x_linear(xright, x_low, x_high, Nk, Mk, mom)
    for i in range(1, N_TRAP):
        xi = xleft + dx*i
        Tot_add += cn_x_linear(xi, x_low, x_high, Nk, Mk, mom)
    Tot_add *= dx

    #if mom == 0.:
    #    if Tot_add >= Nk:
    #        print "TRAP BAD", Tot_add, Nk, xleft, xright
    #        print "        ", cn_x_linear(xleft, x_low, x_high, Nk, Mk, mom, 1)
    #        print "        ", cn_x_linear(xright, x_low, x_high, Nk, Mk, mom, 1)
    
    return Tot_add

@cython.cdivision(True)
def adjust_bins(double[::1] xks, double[::1] dms,
                double[::1] Nks, double[::1] Mks, double[::1] Mks_dry,
                int nk, int output_log):
    new_Nks = np.zeros(shape=(nk), dtype='d')
    new_Mks = np.zeros(shape=(nk), dtype='d')
    new_Mks_dry = np.zeros(shape=(nk), dtype='d')

    cdef:
        double xk, xkp1, Nk, Mk, Nk_dry, Mk_dry
        double x_low, x_high, xm, xmp1
        double N_acc, M_acc
        double Nks_add, Mks_add
        double f_N, f_M, f_aer
        unsigned int m, k

    for k in range(nk):
        xk = xks[k]
        xkp1 = xks[k+1]

        Nk = Nks[k]
        Mk = Mks[k]
        Nk_dry = Nks[k]
        Mk_dry = Mks_dry[k]

        if (Nk == 0.) or (Mk == 0.):
            if output_log > 0:
                print "k", k, "no mass in bin"
            continue
        
        dm = dms[k]

        Nk_new = Nk*1.0
        Mk_new = Mk + Nk*dm
        
        x_low = xk + dm
        x_high = xkp1 + dm

        if output_log > 0:
            print "k", k, "(%1.5e, %1.5e)->(%1.5e, %1.5e)" % (xk, xkp1, x_low, x_high)
            print "   dm", dm
            print "   Mk", Mk, Mk_new
            print "   Nk", Nk, Nk_new      
        
        ## Loop over bins; if a bin is partially contained with [x_low, x_high],
        ## then integrate over the linear approximation in that bin to compute the
        ## shift in moments
        N_acc, M_acc = 0., 0.
        for m in range(nk):
            xm = xks[m]
            xmp1 = xks[m+1]

            if (xm <= x_low < xmp1 <= x_high):
                xleft, xright = x_low, xmp1
                if output_log > 1:
                    print "   CASE 1 (xm <= x_low < xmp1 <= x_high)"
            elif (x_low <= xm < x_high < xmp1):
                xleft, xright = xm, x_high
                if output_log > 1:
                    print "   CASE 2 (x_low <= xm < x_high < xmp1)"
            elif (x_low <= xm) and (xmp1 <= x_high):
                xleft, xright = xm, xmp1
                if output_log > 1:
                    print "   CASE 3 (x_low <= xm and xmp1 <= x_high)"
            elif (xm < x_low) and (x_high < xmp1):
                xleft, xright = x_low, x_high
                if output_log > 1:
                    print "   CASE 4 (xm < x_low and x_high < xmp1)"
            else: 
                continue

            """
            if (xm <= x_low) and (x_low < xmp1 <= x_high):
                ## Contribute xlow->xmp1 to bin m
                xleft, xright = x_low, xmp1
                if output_log > 1:
                    print "   CASE 1 (xm <= x_low, x_low < xmp1 <= x_high)"
            elif (x_low <= xm < xmp1 <= x_high):
                ## Contribute xm->xmp1 to bin m
                xleft, xright = xm, xmp1
                if output_log > 1:
                    print "   CASE 2 (x_low <= xm < xmp1 <= x_high)"
            elif (xm <= x_high) and (x_high <= xmp1):
                ## Contribute xmp1->xhigh to bin m
                xleft, xright = np.max([xm, x_low]), x_high
                if output_log > 1:
                    print "   CASE 3 (xm <= x_high <= xmp1) "
            else:
                continue
            """

            ## Wet Droplets            
            Nks_add = trap_integrate(xleft, xright, x_low, x_high, Nk_new, Mk_new, 0)
            Mks_add = trap_integrate(xleft, xright, x_low, x_high, Nk_new, Mk_new, 1)
            new_Nks[m] += Nks_add
            new_Mks[m] += Mks_add
        
            N_acc += Nks_add
            M_acc += Mks_add
            
            if output_log > 1:
                print "   water -"
                print "   added %1.5e/cc to bin %d" % (Nks_add, m)
                print "   added %1.5e kg to bin %d, (%1.3e, %1.3e)" % (Mks_add, m, xm, xmp1)
                print "      after integrating from (%1.3e, %1.3e)" % (xleft, xright)
        
            ## Dry particles
            f_N = Nks_add/Nk_new
            f_M = Mks_add/Mk_new
            f_aer = a*f_M + (1. - a)*f_N
            new_Mks_dry[m] += Mk_dry*f_aer            
            
            if output_log > 1:
                print "   aerosol -"
                print "   added %1.5e/cc to bin %d" % (Nks_add, m)
                print "   added %1.5e kg to bin %d" % (Mk_dry*f_aer, m)
            
        if output_log > 0:
            print "+++ %1.5e (%1.5e/%1.5e) total N acc" % (N_acc/Nk_new, N_acc, Nk_new)
            print "+++ %1.5e (%1.5e/%1.5e) total M acc" % (M_acc/Mk_new, M_acc, Mk_new)
            print "---------------------------------"

    ## Final sanity check is something negative? Shouldn't be!
    cdef double N_tot = 0.0
    cdef double new_N_tot = 0.0
    if output_log > 0:
        for k in range(nk):
            print k, "%1.3e -> %1.3e" % (Nks[k], new_Nks[k])
            N_tot += Nks[k]
            new_N_tot += new_Nks[k]
        print "  ---------------  "
        print "%1.3e -> %1.3e" % (N_tot, new_N_tot)

    for k in range(nk):
        if new_Nks[k] < 0: new_Nks[k] = 0
        if new_Mks[k] < 0: new_Mks[k] = 0
        if new_Mks_dry[k] < 0: new_Mks_dry[k] = 0

    return new_Nks, new_Mks, new_Mks_dry

@cython.cdivision(True)
def der(double t, double[::1] y,
        double[::1] Nks, double[::1] Mks, double[::1] Mks_dry, 
        double V, double kappa, double rho, int nk,
        int output_log=0):

    cdef double z = y[0]
    cdef double P = y[1]
    cdef double T = y[2]
    cdef double wv = y[3]
    cdef double wc = y[4]
    cdef double S = y[5]

    cdef double T_c = T-273.15 # convert temperature to Celsius
    cdef double pv_sat = es(T_c) # saturation vapor pressure
    cdef double wv_sat = wv/(S+1.) # saturation mixing ratio
    cdef double Tv = (1.+0.61*wv)*T
    cdef double rho_air = P/(Rd*Tv)

    cdef:
        double dz_dt, dP_dt, dwc_dt, dwv_dt, dT_dt, dS_dt

    dP_dt = (-g*P*V)/(Rd*Tv)
    dwc_dt = 0.0

    cdef:
        double Nk, Mk, Nk_dry, Mk_dry
        double mean_m, mean_m_dry, mean_r, mean_r_dry
        double dm_dt

    dms_dt = np.zeros(shape=(nk), dtype='d')

    for k in range(nk):
        Nk = Nks[k]
        Mk = Mks[k]
        Nk_dry = Nks[k]
        Mk_dry = Mks_dry[k]
        
        mean_m = Mk/Nk # kg
        mean_m_dry = Mk_dry/Nk_dry
        mean_r = m_to_r(mean_m, rho_w) # m
        mean_r_dry = m_to_r(mean_m_dry, rho)
        
        if isnan(mean_m) or (mean_m == 0.0):
            if output_log > 0:
                print "k", k, "no mass in bin"
            continue

        dm_dt = growth_rate(mean_r, mean_r_dry, T, P, S, kappa, rho_air)
        dms_dt[k] = dm_dt
        dwc_dt += dm_dt*Nk
    dwc_dt *= 1./rho_air

    dwv_dt = -dwc_dt
    dT_dt = -g*V/Cp - L*dwv_dt/Cp
    dz_dt = V

    ## GHAN (2011)
    cdef double alpha, gamma
    alpha = (g*Mw*L)/(Cp*R*(T**2)) - (g*Ma)/(R*T)
    gamma = (P*Ma)/(Mw*pv_sat) + (Mw*L*L)/(Cp*R*T*T)
    dS_dt = alpha*V - gamma*dwc_dt

    x = np.empty(shape=(nk+6), dtype='d')
    x[0] = dz_dt
    x[1] = dP_dt
    x[2] = dT_dt
    x[3] = dwv_dt
    x[4] = dwc_dt
    x[5] = dS_dt
    x[6:nk+6] = dms_dt
    
    return x