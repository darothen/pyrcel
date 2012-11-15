from pylab import *
ion()

import numpy as np
from micro import *


## Sulfate
Ms = 0.13214
nu = 3.0
epsilon = 0.7
rho_p = 1.769e-3*1e6
kappa = 0.15
r_dry = 0.02 # micron

# Original supersaturation
Seq1 = Seq
# New supersaturation w/ kohler theory
def Seq2(T, r, r_dry, kappa):
    '''Equilibrium supersaturation predicted by Kohler theory'''
    #A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    #B = 1. + kappa*(r_dry/r)**3.
    #return np.exp(A)/B - 1.0
    B = (r**3 - (r_dry**3))/(r**3 - (r_dry**3)*(1.-kappa))
    return np.exp(A)*B - 1.

l, r = np.log10(r_dry), np.log10(5.0)
rs = np.logspace(l, r, 200)

#ss_orig  = np.array(map(lambda r: Seq1(274., r*1e-6, r_dry*1e-6, epsilon, rho_p, Ms, nu), rs))
ss_kappa = np.array(map(lambda r: Seq2(274., r*1e-6, r_dry*1e-6, kappa), rs))

figure(1)
#plot(rs, ss_orig, label="molality")
plot(rs, ss_kappa, label=r"$\kappa = $ %3.2f" % kappa)

legend(loc="upper right")
semilogx()
xlim(1e-3, 1e1)
ylim(-.1, 0.2)

xlabel(r"wet radius ($\mu$m)"); ylabel("supersaturation (%)")

r_crit, s_crit = kohler_crit(274., r_dry*1e-6, kappa)
print r_crit, r_crit/(r_dry*1e-6), s_crit
