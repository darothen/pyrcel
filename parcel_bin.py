from parcel_model import *
import numpy as np
from scipy.integrate import quad, ode
from scipy.optimize import bisect

from pandas import *

import parcel_aux_bin as bin_methods

EQUIL_ITER = 2

s = 2
p = 2**(1./s)
nx = 100*s
r0 = 0.001*1e-6 ## meters
rmax = 100.0*1e-6

P0 = 80000.
T0 = 283.15
V = 0.5
S0 = -0.001
wv0 = (1.-S0)*0.622*es(T0-273.15)/(P0-es(T0-273.15))
aerosol_rho = 1760.

size_dist = Lognorm(mu=0.05*1e-6, sigma=2.0, N=1000.*1e6)
size_dist2 = Lognorm(mu=0.05, sigma=2.0, N=1000.)

#############################

@np.vectorize
def moment_bin(xl, xr, moment, n_x, method='quad'):
    l = moment
    if method == "quad":
        return quad(lambda x: (x**l)*n_x(x), xl, xr)[0]
    else:
        return 0.5*(xr - xl)*((xl**l)*n_x(xl) + (xr**l)*n_x(xr))


x0 = (4./3.)*np.pi*(r0**3)*aerosol_rho ## kg/m^3
xmax = (4./3.)*np.pi*(rmax**3)*aerosol_rho

rks = [r0, ]
xks = [x0, ]
while rks[-1] <= rmax:
	xk = xks[-1]*p
	rk = (xk*0.75/(aerosol_rho*np.pi))**(1./3.) 

	rks.append(rk)
	xks.append(xk)
## Add a giant catch-all bin at the end
for k in xrange(1):
	xks.append(xks[-1]*20.)
	rks.append((xks[-1]*0.75/(aerosol_rho*np.pi))**(1./3.))

xks_edges = np.array(xks)
rks_edges = np.array(rks)
xks = xks_edges[:-1]
rks = rks_edges[:-1]

nk = len(xks)
print "nk =", nk
aerosol = AerosolSpecies("sulfate", size_dist, bins=rks_edges, kappa=0.6, rho=aerosol_rho)
aerosol2 = AerosolSpecies("sulfate", size_dist2, bins=rks_edges*1e6, kappa=0.6, rho=aerosol_rho)

#############################

def n_x_approx(xk, xkp1, Nk, Mk):
    #if Mk < 1e-10: Mk = 0
    #print "(n_x_approx)", "(%e, %e)" % (xk, xkp1), Nk, Mk
    xk_hat = 0.5*(xk + xkp1)
    n0 = Nk/(xkp1 - xk)
    a = 12.*(Mk - xk_hat*Nk)/((xkp1 - xk)**3)
    xi = xk_hat - n0/a
    
    nx_linear = lambda x: n0 + a*(x - xk_hat)
        
    if nx_linear(xk) < 0 and nx_linear(xkp1) < 0: 
        raise ValueError("both endpoints were negative")
        
    if nx_linear(xk) < 0:
        #print "left endpoint negative"
        x_star = 3.*Mk/Nk - 2.*xkp1
        #xstar = xi
        a_star = 2.*Nk/((xkp1 - x_star)**2)
        def nx_linear(x):
            if x < x_star: 
                return 0.
            else:
                return a_star*(x - x_star)        
        #print "(n_x_approx)", xk, "...0...", x_star, ".......", xkp1
        if not (xk < x_star < xkp1): 
            raise ValueError("n(xk) < 0 - x_star out of bounds (%1.3e | %1.3e, %1.3e)" %
                                                      (x_star, xk, xkp1))
        
    if nx_linear(xkp1) < 0:
        #print "right endpoint negative"
        x_star = 3.*Mk/Nk - 2.*xk
        #x_star = xi
        a_star = -2.*Nk/((xk - x_star)**2)
        def nx_linear(x):
            if x > x_star:
                return 0.
            else:
                return a_star*(x - x_star)        
        #print "(n_x_approx)", xk, ".......", x_star, "...0...", xkp1
        if not (xk < x_star < xkp1):
            raise ValueError("n(xkp1) < 0 - x_star out of bounds (%1.3e | %1.3e, %1.3e)" %
                                                      (x_star, xk, xkp1))
        
        
    return np.vectorize(nx_linear)

#############################

T_c = T0-273.15 # convert temperature to Celsius
pv_sat = es(T_c) # saturation vapor pressure
wv_sat = wv0/(S0+1.) # saturation mixing ratio
Tv = (1.+0.61*wv0)*T0
rho_air = P0/(Rd*Tv)

t_end = 20
print aerosol2
pm = ParcelModel([aerosol2,], V,  T0, S0, P0, console=False)
pp, a = pm.run(t_end=t_end, dt=0.01, solver='odeint')
#pp.S.plot()

seed = 20.

if __name__ == "__main__":

	Nks_dry0, Mks_dry0 = np.zeros_like(xks), np.zeros_like(xks)
	xs_to_rs = np.vectorize(lambda x: (x*0.75/(aerosol_rho*np.pi))**(1./3.)) # meters
	pdf_func = np.vectorize(lambda x: size_dist.pdf(xs_to_rs(x)))
	print "INITIAL PASS THROUGH DRY MASS"
	for k in xrange(nk):
		rk = rks_edges[k]
		rkp1 = rks_edges[k+1]
		Nk = moment_bin(rk, rkp1, 0, size_dist.pdf)
		Mk = moment_bin(rk, rkp1, 3, size_dist.pdf)*(4.*np.pi*aerosol_rho/3.)
		print "(%1.3e, %1.3e) %1.3e %1.3e -> %1.3e" % (rk, rkp1, Nk, Mk, xs_to_rs(Mk/Nk))
		Nks_dry0[k] = Nk
		Mks_dry0[k] = Mk
	Nks0 = Nks_dry0[:]

	## Equilibrate the droplets in each bin. To do this, we want to set up a scenario such
	## that when we calculate dm_dt in the timestep function, we retrieve 0 if we hold S
	## steady. Thus, we want to predict Mk in each bin such that the mean droplet size
	## is in steady state
	f = lambda r, r_dry: (Seq(r, r_dry, T0, aerosol.kappa) - S0)
	Mks0 = Mks_dry0[:]
	dms = np.zeros_like(xks)
	wc0 = 0.

	print "EQUILIBRATION CALCULATION"
	for it in xrange(EQUIL_ITER):
		for k in xrange(nk):
			# Attempt 2 - predict dm_dt using method from before, and re-bin
			if Nks0[k] == 0 or Mks_dry0[k] == 0:
				dms[k] = 0.
				continue

			mean_x_wet = Mks0[k]/Nks0[k]
			mean_r_wet = xs_to_rs(mean_x_wet)
			mean_x_dry = Mks_dry0[k]/Nks0[k]
			mean_r_dry = xs_to_rs(mean_x_dry)

			r_b, _ = kohler_crit(T0, mean_r_dry, aerosol.kappa)
			r_a = mean_r_dry

			mean_r_new = bisect(f, r_a, r_b, args=(mean_r_dry, ), xtol=1e-30)
			mean_x_new = (4.*np.pi*aerosol_rho/3.)*(mean_r_dry**3) + \
			   		 	 (4.*np.pi*rho_w/3.)*((mean_r_new**3) - (mean_r_dry**3))

			Mk_new = Nks0[k]*mean_x_new
			#dms[k] = (Mk_new - Mks0[k])
			dms[k] = mean_x_new - mean_x_wet

			#print "k %3d" % k, "%1.3e->%1.3e | %1.3e" % (mean_r_dry, mean_r_new, mean_r_wet)
			#print "      (%1.3e, %1.3e) %1.3e" % (xks[k], p*xks[k], mean_x_new)
			#print "      mass: %1.3e to %1.3e (%1.3e kg/drop)" % (Mks0[k], Mk_new, dms[k])
		Nks1, Mks1, Mks_dry1 = bin_methods.adjust_bins(xks_edges, dms,
													   Nks0, Mks0, Mks_dry0,
													   nk, 0)
		wc_add = np.sum(dms*Nks0)/rho_air
		wc0 += wc_add

		print "Adjustment summary (%d):" % it
		print "       Nks0/Nks_dry0", np.sum(Nks1)/np.sum(Nks_dry0)
		print "           Mks0/Mks0", np.sum(Mks1)/np.sum(Mks0)
		print "   Mks_dry1/Mks_dry0", np.sum(Mks_dry1)/np.sum(Mks_dry0)
		print "   wc0 - ", wc0

		Nks0, Mks0, Mks_dry0 = Nks1, Mks1, Mks_dry1

		#if np.abs(wc_add < EQUIL_WC_TOL):
		#	print "FINISHED EQUIL AFTER %d ITERATIONS" % it
		#	break

	print "  ---------- " 


	rhs = bin_methods.der

	dt = 1.
	y0 = [0.0, P0, T0, wv0, wc0, S0]
	y0.extend(np.zeros_like(xks))
	t0 = 0.
	#t_end = 0.4

	r = ode(rhs)
	log = 0
	r.set_integrator('vode', method='adams', order=5, max_step=dt/100.)
	r.set_initial_value(y0)
	r.set_f_params(*[xks_edges, Nks0, Mks0, Mks_dry0, V, aerosol.kappa, aerosol_rho, nk, log])

	Nks = Nks0[:]
	Mks = Mks0[:]
	Mks_dry = Mks_dry0[:]

	print t0, y0

	t0s = [t0, ]
	y0s = [y0[:6], ]
	while r.successful() and r.t < t_end:
		r.integrate(r.t + dt)
		print r.t, r.y[:6]
		y0s.append(r.y[:6])
		t0s.append(r.t)

		dms = r.y[6:nk+6]*1.0

		Nks, Mks, Mks_dry = bin_methods.adjust_bins(xks_edges, dms,
													Nks, Mks, Mks_dry,
													nk, 0)
		#for rk, N, M in zip(rks, Nks, Mks):
		#	print rk, N, M
		print "     Nk", np.ma.sum(Nks)/np.ma.sum(Nks0)
		print "     Mk", np.ma.sum(Mks)/np.ma.sum(Mks0)
		print " Mk_dry", np.ma.sum(Mks_dry)/np.ma.sum(Mks_dry0)
		#if np.ma.sum(Mks_dry)/np.ma.sum(Mks_dry0) > 1.1 or \
		#	np.ma.sum(Mks_dry)/np.ma.sum(Mks_dry0) < 0.9:
		#	raise ValueError("Too much dry mass change")

		r.set_f_params(*[xks_edges, Nks, Mks, Mks_dry, V, aerosol.kappa, aerosol_rho, nk, log])

		r.y[6:nk+6] = 0.
		#r.y[:6] = y0[:6]

	x = np.array(y0s)

	x = pandas.DataFrame( {'P':x[:,1], 'T':x[:,2], 'wv':x[:,3],
	                       'wc':x[:,4], 'S':x[:,5], 'z':x[:,0]},
	                       index=t0s)

	from pylab import *
	ion()

	figure(1); clf();
	bar(rks, Nks0, edgecolor='k', color='b',
	    width=(rks_edges[1:] - rks), alpha=1., label='before')
	bar(rks, Nks, edgecolor='k', color='y',
	    width=(rks_edges[1:] - rks), alpha=.5, label='after')
	bar(rks, Nks_dry0, edgecolor='k', color='w', hatch="//",
		width=(rks_edges[1:] - rks), alpha=0.5, label="dry")
	semilogx()

	figure(2); clf();
	step(Nks0, 'b',  marker='x', label='before')
	step(Nks,  'y',  marker='o', label='after')
	step(Nks_dry0, 'k', marker='^', label='dry')
	semilogy()
	ylim(1e6, 1e9)
	
	
	figure(3); clf();
	fig, axes = subplots(2, 3, num=3, figsize=(20,8))
	cols = { 'P': (79900, 80001),
			 'T': (283., 284.),
			 'S': (S0, 0.002),
			 'wv': (pp.wv.min(), pp.wv.max()),
			 'wc': (pp.wc.min(), pp.wc.max()),

	}
	for ax, p_key in zip(axes.flatten(), cols.keys()):
		y_lims = cols[p_key]
		ax.set_ylim(y_lims)
		ax.set_title(p_key)
		pp[p_key].plot(ax=ax)
		x[p_key].plot(ax=ax, linewidth=0, marker='x')
		ax.set_xlim(0, t_end)



		




