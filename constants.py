"""
.. module:: parcel
    :synopsis: Useful and commonly used microphysics/thermodynamics constants.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""

## Therymodynamic/chemistry constants
g = 9.81             #: Gravitational constant, m/s^2
Cp = 1004.0          #: Specific heat of dry air at constant pressure, J/kg
L = 2.5e6            #: Latent heat of condensation, J/kg
rho_w = 1e3          #: Density of water, kg/m^3
Rd = 287.0           #: Gas constant for dry air, J/(kg K)
Rv = 461.5           #: Gas constant for water vapor, J/(kg K)
R = 8.314            #: Universal gas constant, J/(mol K)
Mw = 18.0153/1e3     #: Molecular weight of water, kg/mol
Ma = 28.9/1e3        #: Molecular weight of dry air, kg/mol
Dv = 3.e-5           #: Diffusivity of water vapor in air, m^2/s
ac = 1.0             #: condensation constant
Ka = 2.e-2           #: Thermal conductivity of air, J/m/s/K
at = 0.96            #: thermal accomodation coefficient
epsilon = 0.622      #: molecular weight of water / molecular weight of dry air

import pandas as pd
std_atm = pd.read_csv("std_atm.txt", delim_whitespace=True)