import numpy as np
from pylab import *

from parcel_model import *

M, N = 15, 36
N_equ = (M/2)*(M+1) + M*(N-M) + 3
diag_plots = True

P0 = 80000.
T0 = 283.15
V = 0.5
S0 = -0.00
wv0 = (1.-S0)*0.622*es(T0-273.15)/(P0-es(T0-273.15))

print "M (aerosol) = ", M
print "N (droplet) = ", N
print "Total Number of Equations: ", N_equ

aerosol_rho = 1760. # kg/m^3
aerosol_kappa = 0.6

r_aerosol_min = 0.025 # micron
r_droplet_min = 0.025 # micron

x_aerosol_min = (4./3.)*np.pi*((r_aerosol_min*1e-6)**3)*aerosol_rho # kg
x_droplet_min = (4./3.)*np.pi*((r_droplet_min*1e-6)**3)*rho_w

## Calculate bin edges in radius and mass space
xks_aerosol = [x_aerosol_min, ]
for i in xrange(M):
	xks_aerosol.append(2.*xks_aerosol[-1])
xks_aerosol = np.array(xks_aerosol) # kg
rks_aerosol = (xks_aerosol*0.75/(aerosol_rho*np.pi))**(1./3.) # meters

xks_droplet = [x_droplet_min, ]
for i in xrange(N):
	xks_droplet.append(2.*xks_droplet[-1])
xks_droplet = np.array(xks_droplet) # kg
rks_droplet = (xks_droplet*0.75/(rho_w*np.pi))**(1./3.) # meters

print "\nAerosol bin ranges"
print "   radius - [%2.2e, %2.2e] micron" % (rks_aerosol[0]*1e-6, rks_aerosol[-1]*1e-6)
print "     mass - [%2.2e, %2.2e] kg" % (xks_aerosol[0], xks_aerosol[-1])

print "\nDroplet bin ranges"
print "   radius - [%2.2e, %2.2e] micron" % (rks_droplet[0]*1e-6, rks_droplet[-1]*1e-6)
print "     mass - [%2.2e, %2.2e] kg" % (xks_droplet[0], xks_droplet[-1])

## Setup the 2D grid
mass_grid = np.zeros((M, N))
number_grid = np.zeros((M, N))

## Build a simple aerosol size_distribution to place onto the 2D grid
size_dist = Lognorm(mu=0.05, sigma=2.0, N=500)
aerosol = AerosolSpecies("sulfate", size_dist, 
						 kappa=aerosol_kappa, rho=aerosol_rho, bins=rks_aerosol*1e6)
Nis = aerosol.Nis # m^3

if diag_plots:
	figure(1); clf()
	bar(rks_aerosol[:-1], Nis,
		width=rks_aerosol[1:]-rks_aerosol[:-1])
	semilogx()

## Equilibrate the aerosol population
