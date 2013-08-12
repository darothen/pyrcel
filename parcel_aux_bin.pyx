#cython: cdivision=True

from libc.math cimport exp, sqrt, abs
from libc.stdio cimport printf
import numpy as np
import parcel_model.constants as c

from cython_gsl cimport *

cimport cython

## CythonGSL typedefs
ctypedef double * double_ptr
ctypedef void * void_ptr


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
    int N_TRAP = 1100        #: number of steps to use in internal trapezoid rule

cdef inline double cmax(double a, double b) nogil: return a if a >= b else b
cdef inline double cmin(double a, double b) nogil: return a if a <= b else b
    
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

cdef inline double es(double T) nogil:
    return 611.2*exp(17.67*T/(T+243.5))

cdef inline double m_to_r(double m, double rho) nogil:
    if m <= 0.: 
        return 0.
    else:
        return (m*0.75/(rho*PI))**(1./3.)

cdef inline double r_to_m(double r, double rho) nogil:
    if r <= 0.:
        return 0.
    else:
        return (4./3.)*PI*rho_w*(r**3)

cdef double Seq_wrap(double x, void * params) nogil:
    cdef double r_dry, T, kappa, offset

    r_dry = (<double_ptr> params)[0]
    T = (<double_ptr> params)[1]
    kappa = (<double_ptr> params)[2]
    offset = (<double_ptr> params)[3]
    prefac = (<double_ptr> params)[4]

    return Seq(x, r_dry, T, kappa, offset, prefac)

@cython.cdivision(True)
cpdef double Seq(double r, double r_dry, double T, double kappa, double offset=0.0,
                 double prefactor=1.) nogil:
    cdef double A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    cdef double B = 1.0
    if kappa > 0.0:
        B = (r**3 - (r_dry**3))/(r**3 - (r_dry**3)*(1.-kappa))
    cdef double returnval = exp(A)*B - 1.0
    return prefactor*returnval - offset

def kohler_crit_public(double T, double r_dry, double kappa):
    cdef double x = 0.
    with nogil:
        x = kohler_crit(T, r_dry, kappa)
    return x

cdef double kohler_crit(double T, double r_dry, double kappa) nogil:
    cdef int status, it, maxit
    it = 0
    maxit = 100

    cdef double params[6]
    params[0] = r_dry
    params[1] = T
    params[2] = kappa
    params[3] = 0. # offset
    params[4] = -1. # negate

    cdef gsl_min_fminimizer_type * fmin_type = gsl_min_fminimizer_brent
    cdef gsl_min_fminimizer * s = gsl_min_fminimizer_alloc(fmin_type)

    cdef gsl_function F
    F.function = &Seq_wrap
    F.params = &params

    cdef double a = r_dry+1e-10
    cdef double b = r_dry*1e4

    ## nudge over until find a satisfactory initial condition
    cdef double m = a
    cdef double fm = Seq_wrap(m, params)
    cdef double fa = Seq_wrap(a, params)
    cdef double fb = Seq_wrap(b, params)
    cdef int max_space = 1000
    cdef double dx = (b/a)**(1./max_space)
    cdef double counter = 1.
    while (fm >= fa or fm >= fb) and (counter < max_space):
        m = a*(dx**counter)
        fm = Seq_wrap(m, params)
        counter += 1.
    #printf("total moves - %f", counter)

    '''
    #cdef int nr = 100
    #cdef double ri = 1.
    #cdef double dr = (b-a)/nr
    #while ri < nr:
    #    printf("%e %e\n", a+ri*dr, Seq_wrap(a+ri*dr, params))
    #    ri += 1.

    printf("\n>> import numpy as np; from pylab import *; ion()\n")
    printf(">> from parcel_aux_bin import Seq\n")
    printf(">> a, b, m= %e, %e, %e\n", a, b, m)
    printf(">> rs = np.logspace(np.log10(a), np.log10(b), 200)\n")
    printf(">> plot(rs, [Seq(r, %e, %e, %e, 0., -1.) for r in rs])\n", r_dry, T, kappa)
    printf(">> semilogx(); ylim(-1e-1, 1e-1)\n\n")
    printf("%e %e %e\n", fa,fm, fb)
    '''


    status = gsl_min_fminimizer_set(s, &F, m,  a,  b)

    cdef double GSL_CONTINUE = -2
    status = -2
    while (it < maxit and status == GSL_CONTINUE):
        status = gsl_min_fminimizer_iterate(s)

        m = gsl_min_fminimizer_x_minimum(s)
        a = gsl_min_fminimizer_x_lower(s)
        b = gsl_min_fminimizer_x_upper(s)

        status = gsl_min_test_interval(a, b, 0., 1e-15)
        if (status == GSL_SUCCESS):
            break

        it += 1

    gsl_min_fminimizer_free(s)

    return m

cpdef double r_eq_find(double S, double low, double high,
                       double mean_r_dry, double T, double kappa,
                       int maxit=100, double tol=1e-15) nogil:

    cdef int it = 0

    cdef double params[4]
    params[0] = mean_r_dry
    params[1] = T
    params[2] = kappa
    params[3] = S
    params[4] = 1.

    cdef gsl_function F
    F.function = &Seq_wrap
    F.params = &params

    cdef gsl_root_fsolver_type * fsol_type = gsl_root_fsolver_brent
    cdef gsl_root_fsolver * s = gsl_root_fsolver_alloc(fsol_type)

    cdef int status
    status = gsl_root_fsolver_set(s, &F, low, high)

    cdef double r = 0.0

    #printf("Using %s method\n", gsl_root_fsolver_name(s))
    #printf("%5s [%9s, %9s] %9s %10s %9s\n",
    #        "iter", "lower", "upper", "root", "err", "err(est)")

    #printf("initial sanity check\n")
    #printf("%e, %e | %e, %e\n", low, high, Seq(low, mean_r_dry, T, kappa, S), Seq(high, mean_r_dry, T, kappa, S))

    cdef double GSL_CONTINUE = -2
    status = -2
    while (it < maxit and status == GSL_CONTINUE) :
        status = gsl_root_fsolver_iterate(s)
        r = gsl_root_fsolver_root(s)
        low = gsl_root_fsolver_x_lower(s)
        high = gsl_root_fsolver_x_upper(s)
        status = gsl_root_test_interval(low, high, 0., tol)

        if status == GSL_SUCCESS:
        #    print "Converged!"
            break

        #printf("%5d [%1.5e, %1.5e] %1.5e %+1.5e %1.5e\n", \
        #        it, low, high, r, 0.,  high-low)

        it += 1

    return r

cdef double growth_rate(double r, double r_dry, double T, double P, double S,
                        double kappa, double rho_air, double dt=0.) nogil:
    cdef double pv_sat = es(T - 273.15)
    cdef:
        double G_a, G_b, G
        double delta_S, dr_dt
        double dv_r, ka_r, P_atm, Seq_r
        double dm_dyn, dm_eq
        double equil_r
        double r_low, r_high, r_crit
    
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
    # only do this if the droplet is large and dilute
    if r > 1e-6:
        return dm_dyn 

    if dt > 0:
        if delta_S >= 0:
            r_crit = kohler_crit(T, r_dry, kappa) 
            r_low, r_high = r, r_crit

            if Seq(r_high, r_dry, T, kappa, S)*Seq(r_high,r_dry,T,kappa,S) > 0:
                ## same sign, won't be able to find an equilibrium value.
                ## FUDGE FACTOR
                return dm_dyn

            #printf("[%1.3e, %1.3e] - (%1.3e %1.3e)\n", r_low, r_high,
            #        Seq(r_low, r_dry, T, kappa, S), Seq(r_high, r_dry, T, kappa, S))     
            #printf("   Seq(r, %e, %e, %e, %e)", r_dry, T, kappa, S)       
            equil_r = r_eq_find(S, r_low, r_high, r_dry, T, kappa)

            dr_dt = (equil_r - r)/dt
            dm_eq = 4.*PI*rho_w*r*r*dr_dt

            #printf("gr %1.5e | (%1.5e) %1.5e %1.5e\n", r, equil_r, dm_eq, dm_dyn)
            return cmin(dm_eq, dm_dyn)

        else:
            #r_low, r_high = r_dry, r

            #printf("[%1.3e, %1.3e] - (%1.3e %1.3e)\n", r_low, r_high,
            #       Seq(r_low, r_dry, T, kappa), Seq(r_high, r_dry, T, kappa))  
            #equil_r = r_eq_find(S, r_low, r_high, r_dry, T, kappa)
            #dr_dt = (equil_r - r)/dt
            #dm_eq = 4.*PI*rho_w*r*r*dr_dt

            equil_r = r_dry
            dr_dt =  (r_dry - r)/dt
            dm_eq = 4.*PI*rho_w*r*r*dr_dt

            #printf("ev %1.5e | (%1.5e) %1.5e %1.5e\n", r, equil_r, dm_eq, dm_dyn)
            return cmax(dm_eq, dm_dyn)

        #else
        #    r_low, r_high
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

cdef double bin_growth_integrate(double mean_r_dry, double T, double P, double S, 
                                 double kappa, double rho_air,
                                 double xk, double xkp1, double Nk, double Mk):

    cdef gsl_integration_workspace * w
    cdef double result, error
    cdef double params[9]
    w = gsl_integration_workspace_alloc(1000)

    cdef gsl_function F
    F.function = &bin_growth

    params[0] = mean_r_dry
    params[1] = T
    params[2] = P
    params[3] = S
    params[4] = kappa
    params[5] = rho_air
    params[6] = xk
    params[7] = xkp1
    params[8] = Nk
    params[9] = Mk

    F.params = params

    gsl_integration_qags(&F, xk, xkp1, 0., 0., 1000, w, &result, &error)
    return result    

cdef double bin_growth(double x, void * params) nogil:

    cdef double mean_r_dry, T, P, S, kappa, rho_air
    cdef double xk, xkp1, Nk, Mk
    mean_r_dry = (<double_ptr> params)[0]
    T = (<double_ptr> params)[1]
    P = (<double_ptr> params)[2]
    S = (<double_ptr> params)[3]
    kappa = (<double_ptr> params)[4]
    rho_air = (<double_ptr> params)[5]
    xk = (<double_ptr> params)[6]
    xkp1 = (<double_ptr> params)[7]
    Nk = (<double_ptr> params)[8]
    Mk = (<double_ptr> params)[9]

    cdef double Tot_add = 0.0
    cdef double r = m_to_r(x, rho_w)
    cdef double G_r, N_x

    r = m_to_r(x, rho_w)
    if r <= mean_r_dry: 
        return 0.
    else:
        G_r = growth_rate(r, mean_r_dry, T, P, S, kappa, rho_air)
        N_x = cn_x_linear(x, xk, xkp1, Nk, Mk, 0)
        return G_r*N_x

cdef double cn_x_linear_wrap(double x, void * params) nogil:
    cdef double xk, xkp1, Nk,  Mk, moment
    xk = (<double_ptr> params)[0]
    xkp1 = (<double_ptr> params)[1]
    Nk = (<double_ptr> params)[2]
    Mk = (<double_ptr> params)[3]
    moment = (<double_ptr> params)[4]

    return cn_x_linear(x, xk, xkp1, Nk, Mk, moment, 0)

@cython.cdivision(True)
cdef double cn_x_linear(double x, double xk, double xkp1, double Nk, double Mk, 
                         double moment=0.0, int flag=0) nogil:
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
            #if flag == 1: 
            #    print "HERE (n_xk<0)", x, x_star, x-x_star
            #    print "    ", x**moment, a_star, xk, x_star, x
            return (x**moment)*a_star*(x - x_star)
        
    if n_xkp1 < 0:
        x_star = 3.*Mk/Nk - 2.*xk
        a_star = -2.*Nk/((xk - x_star)**2)

        if x >= x_star :
            return 0.
        else:
            #if flag == 1: 
            #    print "HERE (n_xkp1<0)", x, x_star, x-x_star
            #    print "    ", x**moment, a_star, xk, x_star, x
            return (x**moment)*a_star*(x - x_star)    
        
    return (x**moment)*(n0 + a*(x - xk_hat))

cdef double quad_integrate(double xleft, double xright, # bounds to integrate
                           double x_low, double x_high, # params for new linear func
                           double Nk, double Mk, int moment):

    cdef gsl_integration_workspace * w
    cdef double result, error
    cdef double mom = moment*1.0
    cdef double params[4]
    w = gsl_integration_workspace_alloc(1000)

    cdef gsl_function F
    F.function = &cn_x_linear_wrap

    params[0] = x_low
    params[1] = x_high
    params[2] = Nk
    params[3] = Mk
    params[4] = mom
    F.params = params

    gsl_integration_qags(&F, x_low, x_high, 0., 0., 1000, w, &result, &error)
    return result    

## Deprecated, 8/8/2013 after introduction of GSL quadrature routines
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
    '''
    Workaround 8/7/2013: Don't add the edges to the integration! see note in cn_x_linear
    '''
    #Tot_add = 0.5*cn_x_linear(xleft, x_low, x_high, Nk, Mk, mom) + \
    #          0.5*cn_x_linear(xright, x_low, x_high, Nk, Mk, mom)
    for i in range(1, N_TRAP):
        xi = xleft + dx*i
        Tot_add += cn_x_linear(xi, x_low, x_high, Nk, Mk, mom)
    Tot_add *= dx

    #if mom == 0.:
    #    if abs((Tot_add - Nk)/Nk) > 0.1:
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
        double N_acc, M_acc, Md_acc
        double Nks_add, Mks_add, Mks_dry_add
        double f_N, f_M, f_aer
        int last_bin_added
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

        ## Check if things were made too small and adjust
        # - should no longer be triggered 8/8/2013
        if Mk_new < 0:
            dm = (Mk_dry - Mk)/Nk
            Mk_new = Mk_dry

            x_low = xk + dm
            if x_low < 0:
                x_low = xks[0]
            x_high = xkp1 + dm

            if output_log > 0:
                print "   (too much evaporation occurred)"

        if output_log > 0:
            print "k", k, "(%1.5e, %1.5e)->(%1.5e, %1.5e)" % (xk, xkp1, x_low, x_high)
            print "   dm", dm
            print "   Mk", Mk, Mk_new
            print "   Nk", Nk, Nk_new      
        
        ## Loop over bins; if a bin is partially contained with [x_low, x_high],
        ## then integrate over the linear approximation in that bin to compute the
        ## shift in moments
        N_acc, M_acc, Md_acc = 0., 0., 0.
        last_bin_added = -1
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

            last_bin_added = m

            ## Wet Droplets            
            #Nks_add = quad_integrate(xleft, xright, x_low, x_high, Nk_new, Mk_new, 0)
            #Mks_add = quad_integrate(xleft, xright, x_low, x_high, Nk_new, Mk_new, 1)            
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
            Mks_dry_add = Mk_dry*f_aer     

            new_Mks_dry[m] += Mks_dry_add

            Md_acc += Mks_dry_add  
            
            if output_log > 1:
                print "   aerosol -"
                print "   added %1.5e/cc to bin %d" % (Nks_add, m)
                print "   added %1.5e kg to bin %d" % (Mk_dry*f_aer, m)

        ## Adjust for numerical spillage from bins
        
        if Nk_new - N_acc > 0:
            new_Nks[last_bin_added] += Nk_new - N_acc
            N_acc += Nk_new - N_acc
        if Mk_new - M_acc > 0:
            new_Mks[last_bin_added] += Mk_new - M_acc
            M_acc += Mk_new - M_acc
        if Mk_dry - Md_acc > 0:
            new_Mks_dry[last_bin_added] += Mk_dry - Md_acc
            Md_acc += Mk_dry - Md_acc
        
            
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
def der(double[::1] y, double t, double[::1] xks,
        double[::1] Nks, double[::1] Mks, double[::1] Mks_dry, 
        double V, double kappa, double rho, int nk, double dt,
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
        double dm_dt, xk, xkp1, dx, Tot_add, dmk_dt, xi, ri
        unsigned int i, k

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

        dm_dt = growth_rate(mean_r, mean_r_dry, T, P, S, kappa, rho_air, dt)
        dms_dt[k] = dm_dt 

        
        ## To more accurately estimate the change in cloud droplet water, 
        ## integrate the growth rate over the bin
        xk = xks[k]
        xkp1 = xks[k+1]
        #Tot_add = bin_growth_integrate(mean_r_dry, T, P, S, kappa, rho_air,
        #                               xk, xkp1, Nk, Mk)
        
        Tot_add = 0.0
        dx = (xkp1 - xk)/(100.)
        for i in range(1, 100):
            xi = xk + dx*i
            ri = m_to_r(xi, rho_w)
            if ri <= mean_r_dry: continue

            dmk_dt = growth_rate(ri, mean_r_dry, T, P, S, kappa, rho_air)
            Tot_add += dmk_dt*cn_x_linear(xi, xk, xkp1, Nk, Mk, 0)
        Tot_add *= dx
        
        dwc_dt += Tot_add
        
        

        #dwc_dt += dm_dt*Nk
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