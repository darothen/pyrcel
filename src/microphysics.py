"""
"""

import numpy as np
from aerocalc.std_atm import sat_press

Cp = 1004.0   # Specific heat at constant pressure, J/(kg*K)
Rd = 287.0    # Universal gas constant, J/(kg*K)
Rv = 461.5    # Specific gas constant, Water Vapor, J(kg*K)
L = 2.5e6     # Latent heat of vaporization, J/kg
g = 9.8       # gravity constant, m/(s**2)
gamma_d = 9.8 # dry adiabatic lapse rate, C/km

##
epsilon = Rd/Rv

def calc_sat_mixing_ratio(temperature, pressure):
    """
    Compute the saturation mixing ratio given temperature (in K) and pressure (in Pa)
    """
    es = sat_press(T=temperature, RH=1.0, temp_units="K", press_units="pa")
    ws = (epsilon*es)/temperature
    return ws

def calc_virtual_temp(temperature, q):
    """
    Calculate virtual temperature from temperature (K) and specific humidity (kg/kg)
    """
    return (1.0 + 0.61*q)*temperature

    


    
    