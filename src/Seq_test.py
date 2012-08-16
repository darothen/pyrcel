import numpy as np

## Constants
V = 0.01 # Updraft velocity, m/s
g = 9.81 # Gravitational constant, m/s**2
Cp = 1004.0 # Specific heat of dry air at constant pressure, J/kg 
L = 2.5e6 # Latent heat of condensation, J/kg
rho_w = 1e3 # Density of water, kg/m**3
Rd = 287.0 # Gas constant for dry air, J/(kg K)
R = 8.314 # Universal gas constant, J/(mol K)
Mw = 18.0153*1e-3 # Molecular weight of water, kg/mol
sigma_w = lambda T: 0.0761 - 1.55e-4*(T-273.15) # surface tension of water, J/m^2 given T in Kelvin

## Aerosol properties
# Use Sulfate particles, like in Steele thesis
#r_dry = 0.01*1e-6 # m
r_dry = 0.1*1e-6 # m
d_dry = 2.*r_dry
Ms = 0.13214*10. # Molecular weight, kg/mol
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

def Seq(T, r):
    '''Equilibrium supersaturation predicted by Kohler theory'''
    A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    B = (3.*ns*Mw*nu)/(4.*np.pi*rho_w*(r**3. - r_dry**3.))
    return np.exp(A - B) - 1.0

def Seq2(T, r):
	A = (2.*Mw*sigma_w(T))/(T*R*rho_w)
	B = 4.3*nu*epsilon*(rho_p*(4.*np.pi*(r_dry**3)/3.))/Ms
	return A/r - B/(r**3)

def Seqa(T, r):
	return np.exp((2.*Mw*sigma_w(T))/(R*T*rho_w*r))

def Seqb(T, r):
	return np.exp((3.*ns*Mw*nu)/(4.*np.pi*rho_w*(r**3. - r_dry**3.)))

T = 293.


from pylab import *
clf()
rs = np.logspace(-7,-5,1000)[1:]
subplot(1,3,1)
xyz = np.array([Seq(T, rrr) for rrr in rs])
xyz2 = np.array([Seq2(T, rrr) for rrr in rs])
plot(rs, xyz)
plot(rs, xyz2)
semilogx()
ylim(-0.2, .05)
hlines([0.,], 1e-8, 1e-5, linestyle='dotted')

subplot(1,3,2)
xyz_a = np.array([Seqa(T, rrr) for rrr in rs])
xyz2_a = np.array([1. + (2.*Mw*sigma_w(T))/(T*R*rho_w)/(rrr) for rrr in rs])
plot(rs, xyz_a)
plot(rs, xyz2_a)
semilogx()
hlines([1., ], 1e-8, 1e-5, linestyle='dotted')
ylim(0.8, 1.05)

subplot(1,3,3)
xyz_b = np.array([Seqb(T, rrr) for rrr in rs])
plot(rs, xyz_b)
semilogx()
hlines([1., ], 1e-8, 1e-5, linestyle='dotted')
#ylim(0.8, 1.05)

show()