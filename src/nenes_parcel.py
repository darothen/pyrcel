'''
Cloud Parcel Model based on Nenes et al, 2001
'''

import sys
import numpy as np
import pandas
from scipy.optimize import fmin

## Constants
V = 0.01 # Updraft velocity, m/s
g = 9.81 # Gravitational constant, m/s**2
Cp = 1004.0 # Specific heat of dry air at constant pressure, J/kg 
L = 2.5e6 # Latent heat of condensation, J/kg
rho_w = 1e3 # Density of water, kg/m**3
Rd = 287.0 # Gas constant for dry air, J/(kg K)
R = 8.314 # Universal gas constant, J/(mol K)
Mw = 18.0153/1e3 # Molecular weight of water, kg/mol
sigma_w = lambda T: 0.0761 - 1.55e-4*(T-273.15) # surface tension of water, J/m^2 given T in Kelvin

## NOT CORRECTING FOR NON-CONTINUUM EFFECTS
Dv = 0.3/1e4 # Diffusivity of water vapor in air, m^2/s
ka = lambda T: 419.*(5.69 + 0.017*(T-273.15))*1e-5 # thermal conductivty of air, W/(m K) given T in Kelvin

## Aerosol properties
Ni = 3000.*1e6 # m**-3
ri = 0.1*1e-6 # m

# Use Sulfate particles, like in Steele thesis
#r_dry = 0.01*1e-6 # m
#r_dry = 0.1*1e-6 # m
r_dry = ri
d_dry = 2.*r_dry
Ms = 0.13214 # Molecular weight, kg/mol
rho_s = 2.16*1e-3*1e6
rho_u = 2.17*1e-3*1e6

##1
epsilon = 0.01 # mass fraction of soluble material in the dry particle
rho_p = rho_u/(1.-epsilon*(1.-(rho_u/rho_s)))
##2
#rho_p = 2.1699*1e3
#epsilon = (1. - (rho_u/rho_p))/(1.-(rho_u/rho_s))

nu = 2.0 # number of ions a solute dissolves into
ns = (4.*rho_p*np.pi*epsilon*(d_dry**3.))/(3.*Ms)
# Start with a single bulk aerosol - all same size. Pick clean continental Nuclei mode


## AUXILIARY FUNCTIONS
def es(T):
    """Returns saturation vapor pressure (Pascal) at temperature T (Celsius)
    Formula 2.17 in Rogers&Yau"""
    return 611.2*np.exp(17.67*T/(T+243.5))

def Seq(T, r):
    '''Equilibrium supersaturation predicted by Kohler theory'''
    A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    B = (3.*ns*Mw*nu)/(4.*np.pi*rho_w*(r**3 - r_dry**3))
    return np.exp(A - B) - 1.0

## DERIVATIVE
def der(t, y):
    #print t
    P, T, wv, wc, r, S = y
    
    pv_sat = es(T-273.15) # saturation vapor pressure
    wv_sat = wv/(S+1) # saturation mixing ratio
    Tv = (1.+0.61*wv)*T
    rho_air = P/(Rd*Tv)
    
    # 1) dP_dt 
    dP_dt = (-g*P*V)/(Rd*T)
    
    # 2) dr_dt
    G_a = (rho_w*R*T)/(pv_sat*Dv*Mw)
    G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka(T)*T)
    G = 1./(G_a + G_b)
    
    dr_dt = (G/r)*(S - Seq(T, r))
        
    # 3) dwc_dt
    dwc_dt = (4.*np.pi*rho_w/rho_air)*(Ni *( r**2) * dr_dt)
    
    # 4) dwv_dt
    dwv_dt = -dwc_dt
    
    # 5) dT_dt 
    dT_dt = -g*V/Cp - L*dwv_dt/Cp
    
    # 6) dS_dt
    # Used eq 12.28 from Pruppacher and Klett in stead of (9) from Nenes et al, 2001
    S_a = P*dwv_dt/(0.622*pv_sat)
    S_b = 1.+S
    S_c = 0.622*L*dT_dt/(R*T**2) + g*V/(R*T)
    dS_dt = S_a - S_b*S_c
    
    return [dP_dt, dT_dt, dwv_dt, dwc_dt, dr_dt, dS_dt]

    
if __name__ == "__main__":
    
    ## Initial conditions
    P0 = 100000. # Pressure, Pa
    T0 = 293.15 # Temperature, K
    wv0 = 0.98*0.622*es(T0-273.15)/P0 # Water Vapor mixing ratio, kg/kg. Change the leading term, RH
    S0 = -0.02 # Supersaturation. 1-RH from wv term.
    # Need to compute r and wc
    f = lambda r: Seq(T0, r) - S0
    r0 = fmin(f, ri+ri/100.)[0]
    wc0 = ((4.*np.pi/3.)*(r0**3 - ri**3))*rho_w*Ni

    y0 = [P0, T0, wv0, wc0, r0, S0]
    print y0

    t0 = 0.
    #t_end = 10./V
    t_end = 200./V
    dt = 1.
    t = np.arange(t0, t_end+dt, dt)
    print len(t), t_end

    raw_input("continue?")
    
    from scipy.integrate import ode
    r = ode(der).set_integrator("vode", method="bdf", nsteps=5000, order=5)
    r.set_initial_value(y0)

    out = []
    while r.successful() and r.t < t_end:
        print "1", r.t, r.y
        r.integrate(r.t+dt)
        #print "11", r.t, r.y
        out.append(r.y)
        if r.t%60 == 0: print r.t
        if not r.successful: print "FUCK"; exit()
    x = np.array(out)
    
    if len(x) == len(t)-1: sys.exit()
    '''
    from scipy.integrate import odeint
    x = odeint(der, y0, t)
    '''

    heights = t*V
    df = pandas.DataFrame( {'P':x[:,0], 'T':x[:,1], 'wv':x[:,2],
                            'wc':x[:,3], 'r':x[:,4], 'S':x[:,5]}, index=heights )


    fig = figure()    
    df.plot(subplots=True, figure=fig)
    
    
    