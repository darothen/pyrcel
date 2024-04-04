""" Commonly used constants in microphysics and aerosol thermodynamics equations as
well as important model parameters.

================= ============= ========== ==========        ======================
Symbol            Variable      Value      Units             Description
================= ============= ========== ==========        ======================
:math:`g`         ``g``         9.8        m s**-2           gravitational constant
:math:`C_p`       ``Cp``        1004.0     J/kg              specific heat of dry air
                                                             at constant pressure
:math:`\\rho_w`    ``rho_w``     1000.0     kg m**-3          density of water at STP
:math:`R_d`       ``Rd``        287.0      J/kg/K            gas constant for dry air
:math:`R_v`       ``Rv``        461.5      J/kg/K            gas constant for water vapor
:math:`R`         ``R``         8.314      J/mol/K           universal gas constant
:math:`M_w`       ``Mw``        0.018      kg/mol            molecular weight of water
:math:`M_a`       ``Ma``        0.0289     kg/mol            molecular weight of dry air
:math:`D_v`       ``Dv``        3e-5       m**2/s            diffusivity of water vapor
                                                             in air
:math:`L_v`       ``L``         2.25e6     J/kg/K            latent heat of vaporization
                                                             of water
:math:`\\alpha_c`  ``ac``        1.0        unitless          condensation coefficient
:math:`K_a`       ``Ka``        0.02       J/m/s/K           thermal conductivity of air
:math:`a_T`       ``at``        0.96       unitless          thermal accommodation
                                                             coefficient
:math:`\epsilon`  ``epsilon``   0.622      unitless          ratio of :math:`M_w/M_a`
================= ============= ========== ==========        ======================

Additionally, a reference table containing the
`1976 US Standard Atmosphere <http://www.pdas.com/atmos.html>`_ is implemented in the
constant ``std_atm``, which is a pandas DataFrame with the fields

- ``alt``, altitude in km
- ``sigma``, ratio of density to sea-level density
- ``delta``, ratio of pressure to sea-level pressure
- ``theta``, ratio of temperature to sea-level temperature
- ``temp``, temperature in K
- ``press``, pressure in Pa
- ``dens``, air density in kg/m**3
- ``k.visc``, air kinematic viscosity
- ``ratio``, ratio of speed of sound to kinematic viscosity in m**-1

Using default pandas functons, you can interpolate to any reference pressure or
height level.

"""

import pandas as pd
import pkg_resources

g = 9.81  #: Gravitational constant, m/s^2
Cp = 1004.0  #: Specific heat of dry air at constant pressure, J/kg
L = 2.25e6  #: Latent heat of condensation, J/kg
rho_w = 1e3  #: Density of water, kg/m^3
R = 8.314  #: Universal gas constant, J/(mol K)
Mw = 18.0 / 1e3  #: Molecular weight of water, kg/mol
Ma = 28.9 / 1e3  #: Molecular weight of dry air, kg/mol
Rd = R / Ma  #: Gas constant for dry air, J/(kg K)
Rv = R / Mw  #: Gas constant for water vapor, J/(kg K)
Dv = 3.0e-5  #: Diffusivity of water vapor in air, m^2/s
ac = 1.0  #: condensation constant
Ka = 2.0e-2  #: Thermal conductivity of air, J/m/s/K
at = 0.96  #: thermal accomodation coefficient
epsilon = 0.622  #: molecular weight of water / molecular weight of dry air

# Additional fixed model parameters
N_STATE_VARS = 7
STATE_VARS = ["z", "P", "T", "wv", "wc", "wi", "S"]
STATE_VAR_MAP = {var: i for i, var in enumerate(STATE_VARS)}

# Read the standard atmosphere CSV file
_std_atm_fn = pkg_resources.resource_filename("pyrcel", "data/std_atm.csv")
std_atm = pd.read_csv(_std_atm_fn, sep="\s+")
