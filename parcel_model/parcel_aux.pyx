#cython: embedsignature:True
#cython: profile=False
#cython: boundscheck=False
#cython: cdivision=True
# (embedsignature adds doc-strings accessible by Sphinx)
""" Parcel model derivative calculations implemented in Cython.

"""

from cython.parallel cimport prange, parallel
from libc.math cimport exp, sqrt

import numpy as np
cimport numpy as np
cimport cython
#cimport openmp

import constants as c

## Thermodynamic/chemistry constants
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

## Auxiliary, single-value calculations with GIL released for derivative
## calculations
cdef inline double sigma_w(double T) nogil:
    """See :func:`parcel_model.micro.sigma_w` for full documentation
    """
    return 0.0761 - (1.55e-4)*(T-273.15)

cdef inline double ka(double T, double r, double rho) nogil:
    """See :func:`parcel_model.micro.ka` for full documentation
    """
    cdef double denom, ka_cont
    ka_cont = 1e-3*(4.39 + 0.071*T)
    denom = 1.0 + (ka_cont/(at*r*rho*Cp))*sqrt((2*PI*Ma)/(R*T))
    return ka_cont/denom

cdef inline double dv(double T, double r, double P, double accom) nogil:
    """See :func:`parcel_model.micro.dv` for full documentation
    """
    cdef double denom, dv_cont, P_atm
    P_atm = P*1.01325e-5 # Pa -> atm
    dv_cont = 1e-4*(0.211/P_atm)*((T/273.)**1.94)
    denom = 1.0 + (dv_cont/(accom*r))*sqrt((2*PI*Mw)/(R*T))
    return dv_cont/denom

cdef inline double es(double T):
    """See :func:`parcel_model.micro.es` for full documentation
    """
    return 611.2*exp(17.67*T/(T+243.5))

cdef double Seq(double r, double r_dry, double T, double kappa) nogil:
    """See :func:`parcel_model.micro.Seq` for full documentation.
    """
    cdef double A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    cdef double B = 1.0
    if kappa > 0.0:
        B = (r**3 - (r_dry**3))/(r**3 - (r_dry**3)*(1.-kappa))
    cdef double returnval = exp(A)*B - 1.0
    return returnval

## Jacobian of derivative
def jac(double[::1] y, double t,
        int nr, double[::1] r_drys, double[::1] Nis, double V, double[::1] kappas):

    cdef double z = y[0]
    cdef double P = y[1]
    cdef double T = y[2]
    cdef double wv = y[3]
    cdef double wc = y[4]
    cdef double S = y[5]
    #cdef np.ndarray[double, ndim=1] rs = y[5:]
    cdef double[::1] rs = y[6:]

    cdef double T_c = T-273.15 # convert temperature to Celsius
    cdef double pv_sat = es(T_c) # saturation vapor pressure
    cdef double wv_sat = wv/(S+1.) # saturation mixing ratio
    cdef double Tv = (1.+0.61*wv)*T
    cdef double rho_air = P/(Rd*Tv)

    return 

## RHS Derivative callback function
def der(double[::1] y, double t,
        int nr, double[::1] r_drys, double[::1] Nis, double V, double[::1] kappas,
        double accom):
    """ Calculates the instantaneous time-derivate of the parcel model system.

    Given a current state vector `y` of the parcel model, computes the tendency
    of each term including thermodynamic (pressure, temperature, etc) and aerosol
    terms. The basic aerosol properties used in the model must be passed along
    with the state vector (i.e. if being used as the callback function in an ODE
    solver).

    This function is implemented in NumPy and Python, and is likely *very* slow
    compared to the available Cython version.

    Parameters
    ----------
    y : array_like
        Current state of the parcel model system,
            * y[0] = altitude, m
            * y[1] = pressure, Pa
            * y[2] = temperature, K
            * y[3] = water vapor mass mixing ratio, kg/kg
            * y[4] = droplet liquid water mass mixing ratio, kg/kg
            * y[5] = parcel supersaturation
            * y[`nr`:] = aerosol bin sizes (radii), m
    t : float
        Current simulation time, in seconds.
    nr : Integer
        Number of aerosol radii being tracked.
    r_drys : array_like
        Array recording original aerosol dry radii, m.
    Nis : array_like
        Array recording aerosol number concentrations, 1/(m**3).
    V : float
        Updraft velocity, m/s.
    kappas : array_like
        Array recording aerosol hygroscopicities.

    Returns
    -------
    x : array_like
        Array of shape (`nr`+6, ) containing the evaluated parcel model
        instaneous derivative.

    Notes
    -----
    This Python sketch of the derivative function shouldn't really be used for
    any computational purposes. Instead, see the cythonized version in the auxiliary
    file, **parcel_aux.pyx**. In the default configuration, once the code has been
    built, you can set the environmental variable **OMP_NUM_THREADS** to control
    the parallel for loop which calculates the condensational growth rate for each
    bin.

    """
    cdef double z = y[0]
    cdef double P = y[1]
    cdef double T = y[2]
    cdef double wv = y[3]
    cdef double wc = y[4]
    cdef double S = y[5]
    #cdef np.ndarray[double, ndim=1] rs = y[5:]
    cdef double[::1] rs = y[6:]

    cdef double T_c = T-273.15 # convert temperature to Celsius
    cdef double pv_sat = es(T_c) # saturation vapor pressure
    cdef double wv_sat = wv/(S+1.) # saturation mixing ratio
    cdef double Tv = (1.+0.61*wv)*T
    cdef double rho_air = P/(Rd*Tv)

    ## Begin computing tendencies
    cdef:
        double dP_dt, dwc_dt, dwv_dt, dT_dt, dS_dt
        double[::1] drs_dt, x

    dP_dt = (-g*P*V)/(Rd*Tv)
    dwc_dt = 0.0
    drs_dt = np.empty(shape=(nr), dtype="d")

    cdef: # variables set in parallel loop
        unsigned int i
        double G_a, G_b, G
        double r, r_dry, delta_S, kappa, dr_dt, Ni
        double dv_r, ka_r, P_atm, A, B, Seq_r

    for i in prange(nr, nogil=True):#, schedule='static'):#, num_threads=40):
    #for i in range(nr):
        r = rs[i]
        r_dry = r_drys[i]
        kappa = kappas[i]
        Ni = Nis[i]

        ## Non-continuum diffusivity/thermal conductivity of air near
        ## near particle
        dv_r = dv(T, r, P, accom)
        ka_r = ka(T, r, rho_air)

        ## Condensation coefficient
        G_a = (rho_w*R*T)/(pv_sat*dv_r*Mw)
        G_b = (L*rho_w*((L*Mw/(R*T))-1.))/(ka_r*T)
        G = 1./(G_a + G_b)

        ## Difference between ambient and particle equilibrium supersaturation
        Seq_r = Seq(r, r_dry, T, kappa)
        delta_S = S - Seq_r

        ## Size and liquid water tendencies
        dr_dt = (G/r)*delta_S
        dwc_dt += Ni*r*r*dr_dt # Contribution to liq. water tendency due to growth
        drs_dt[i] = dr_dt
    dwc_dt *= (4.*PI*rho_w/rho_air) # Hydrated aerosol size -> water mass

    dwv_dt = -dwc_dt

    dT_dt = -g*V/Cp - L*dwv_dt/Cp

    dz_dt = V

    ''' Alternative methods for calculation supersaturation tendency
    # Used eq 12.28 from Pruppacher and Klett in stead of (9) from Nenes et al, 2001
    #cdef double S_a, S_b, S_c, dS_dt
    #cdef double S_b_old, S_c_old, dS_dt_old
    #S_a = (S+1.0)

    ## NENES (2001)
    #S_b_old = dT_dt*wv_sat*(17.67*243.5)/((243.5+(Tv-273.15))**2.)
    #S_c_old = (rho_air*g*V)*(wv_sat/P)*((0.622*L)/(Cp*Tv) - 1.0)
    #dS_dt_old = (1./wv_sat)*(dwv_dt - S_a*(S_b_old-S_c_old))

    ## PRUPPACHER (PK 1997)
    #S_b = dT_dt*0.622*L/(Rd*T**2.)
    #S_c = g*V/(Rd*T)
    #dS_dt = P*dwv_dt/(0.622*es(T-273.15)) - S_a*(S_b + S_c)

    ## SEINFELD (SP 1998)
    #S_b = L*Mw*dT_dt/(R*T**2.)
    #S_c = V*g*Ma/(R*T)
    #dS_dt = dwv_dt*(Ma*P)/(Mw*es(T-273.15)) - S_a*(S_b + S_c)
    '''

    ## GHAN (2011)
    cdef double alpha, gamma
    alpha = (g*Mw*L)/(Cp*R*(T**2)) - (g*Ma)/(R*T)
    gamma = (P*Ma)/(Mw*pv_sat) + (Mw*L*L)/(Cp*R*T*T)
    dS_dt = alpha*V - gamma*dwc_dt

    x = np.empty(shape=(nr+6), dtype='d')
    x[0] = dz_dt
    x[1] = dP_dt
    x[2] = dT_dt
    x[3] = dwv_dt
    x[4] = dwc_dt
    x[5] = dS_dt
    x[6:] = drs_dt[:]

    return x
