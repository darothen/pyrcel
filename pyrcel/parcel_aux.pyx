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

## Define double DTYPE 
DTYPE = np.float

## Thermodynamic/chemistry constants
cdef:
    double Mw = c.Mw        #: Molecular weight of water, kg/mol
    double Ma = c.Ma        #: Molecular weight of dry air, kg/mol
    double R = c.R          #: Universal gas constant, J/(mol K)
    double rho_w = c.rho_w  #: Density of water, kg/m**3
    double Rd = c.Rd        #: Gas constant for dry air, J/(kg K)
    double Rv = c.Rv        #: Gas constant for water vapor, J/(kg K)
    double g = c.g          #: Gravitational constant, m/s**2
    double Dv = c.Dv        #: Diffusivity of water vapor in air, m^2/s
    double ac = c.ac        #: Condensation Constant
    double Ka = c.Ka        #: Thermal conductivity of air, J/m/s/K
    double at = c.at        #: thermal accomodation coefficient
    double L = c.L          #: Latent heat of condensation, J/kg
    double Cp = c.Cp        #: Specific heat of dry air at constant pressure, J/kg
    double PI = 3.14159265358979323846264338328 #: Pi, constant

    int N_STATE_VARS = c.N_STATE_VARS #: number of thermodynamic state variables

## Auxiliary, single-value calculations with GIL released for derivative
## calculations
cdef inline double sigma_w(double T) nogil:
    """See :func:`pyrcel.micro.sigma_w` for full documentation
    """
    return 0.0761 - (1.55e-4)*(T-273.15)

cdef inline double ka(double T, double r, double rho) nogil:
    """See :func:`pyrcel.micro.ka` for full documentation
    """
    cdef double denom, ka_cont
    ka_cont = 1e-3*(4.39 + 0.071*T)
    denom = 1.0 + (ka_cont/(at*r*rho*Cp))*sqrt((2*PI*Ma)/(R*T))
    return ka_cont/denom

cdef inline double dv(double T, double r, double P, double accom) nogil:
    """See :func:`pyrcel.micro.dv` for full documentation
    """
    cdef double denom, dv_cont, P_atm
    P_atm = P*1.01325e-5 # Pa -> atm
    dv_cont = 1e-4*(0.211/P_atm)*((T/273.)**1.94)
    denom = 1.0 + (dv_cont/(accom*r))*sqrt((2*PI*Mw)/(R*T))
    return dv_cont/denom

cdef inline double es(double T):
    """See :func:`pyrcel.micro.es` for full documentation
    """
    return 611.2*exp(17.67*T/(T+243.5))

cdef double Seq(double r, double r_dry, double T, double kappa) nogil:
    """See :func:`pyrcel.micro.Seq` for full documentation.
    """
    cdef double A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    cdef double B = 1.0
    if kappa > 0.0:
        B = (r**3 - (r_dry**3))/(r**3 - (r_dry**3)*(1.-kappa))
    cdef double returnval = exp(A)*B - 1.0
    return returnval

## RHS Derivative callback function
def der(double[::1] y, double t,
        int nr, double[::1] r_drys, double[::1] Nis, double V, double[::1] kappas,
        double accom):
    """See :func:`pyrcel.parcel.der` for full documentation

    """
    cdef double z  = y[0]
    cdef double P  = y[1]
    cdef double T  = y[2]
    cdef double wv = y[3]
    cdef double wc = y[4]
    cdef double wi = y[5]
    cdef double S  = y[6]
    cdef double[::1] rs = y[N_STATE_VARS:]

    cdef double T_c = T-273.15 # convert temperature to Celsius
    cdef double pv_sat = es(T_c) # saturation vapor pressure
    cdef double wv_sat = wv/(S+1.) # saturation mixing ratio
    cdef double Tv = (1.+0.61*wv)*T
    cdef double e = (1. + S)*pv_sat # water vapor pressure

    ## Compute air densities from current state
    cdef double rho_air     = P/Rd/Tv
    cdef double rho_air_dry = (P-e)/Rd/T #: TODO - port to parcel.py

    ## Begin computing tendencies
    cdef:
        double dwc_dt, dwv_dt, dwi_dt, dT_dt, dS_dt
        double[::1] drs_dt, x

    dP_dt = -1.*rho_air*g*V

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
    dwc_dt *= (4.*PI*rho_w/rho_air_dry) # Hydrated aerosol size -> water mass
                        # use rho_air_dry for mixing ratio definition consistency
    # No freezing implemented yet
    dwi_dt = 0.0 

    ## MASS BALANCE CONSTRAINT
    dwv_dt = -1.*(dwc_dt + dwi_dt)

    ## ADIABATIC COOLING
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

    x = np.empty(shape=(nr+N_STATE_VARS), dtype='d')
    x[0] = dz_dt
    x[1] = dP_dt
    x[2] = dT_dt
    x[3] = dwv_dt
    x[4] = dwc_dt
    x[5] = dwi_dt
    x[6] = dS_dt
    x[N_STATE_VARS:] = drs_dt[:]

    return x
