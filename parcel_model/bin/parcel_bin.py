from parcel_model import *
import numpy as np
from scipy.integrate import quad, ode
from scipy.optimize import bisect

from pandas import *

import parcel_aux_bin as bin_methods
import odespy

EQUIL_ITER = 2

s = 2
p = 2**(1./s)
nx = 100*s
r0 = 0.001*1e-6 ## meters
rmax = 100.0*1e-6

P0 = 80000.
T0 = 283.15
V = 0.5
S0 = -0.00
wv0 = (1.-S0)*0.622*es(T0-273.15)/(P0-es(T0-273.15))
aerosol_rho = 1760.

size_dist = Lognorm(mu=0.05*1e-6, sigma=2.0, N=500.*1e6)
size_dist2 = Lognorm(mu=0.05, sigma=2.0, N=500.)

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
	xks.append(xks[-1]*p)
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

			#r_b, _ = kohler_crit(T0, mean_r_dry, aerosol.kappa)
			r_b = bin_methods.kohler_crit_public(T0, mean_r_dry, aerosol.kappa)
			r_a = mean_r_dry

			mean_r_new = bisect(f, r_a, r_b, args=(mean_r_dry, ), xtol=1e-30)
			#mean_r2 = bin_methods.r_eq_find(S0, r_a, r_b, mean_r_dry, T0, aerosol.kappa, 
			#								100, 1e-12)
			#print "|--| k %d|" % (k, ), mean_r_new, mean_r2
			mean_x_new = (4.*np.pi*aerosol_rho/3.)*(mean_r_dry**3) + \
			   		 	 (4.*np.pi*rho_w/3.)*((mean_r_new**3) - (mean_r_dry**3))

			Mk_new = Nks0[k]*mean_x_new
			#dms[k] = (Mk_new - Mks0[k])
			dms[k] = mean_x_new - mean_x_wet

			#print "k %3d" % k, "%1.3e->%1.3e | %1.3e" % (mean_r_dry, mean_r_new, mean_r_wet)
			#print "      (%1.3e, %1.3e) %1.3e" % (xks[k], p*xks[k], mean_x_new)
			#print "      mass: %1.3e to %1.3e (%1.3e kg/drop)" % (Mks0[k], Mk_new, dms[k])
		Nks1, Mks1, Mks_dry1 = bin_methods.adjust_bins(xks_edges, dms, dms,
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

	dt = 0.05
	y0 = [0.0, P0, T0, wv0, wc0, S0]
	y0.extend(np.zeros_like(xks))
	t0 = 0.
	#t_end = 2.*dt
	log = 0

	solver = odespy.Euler(rhs)
	solver.set_initial_condition(y0)
	solver.set(f_args=[xks_edges, Nks0, Mks0, Mks_dry0, 
			   		   V, aerosol.kappa, aerosol_rho, nk, dt, log])

	#r = ode(rhs)
	#r.set_integrator('vode', method='bdf', order=5, max_step=dt/100.)
	#r.set_initial_value(y0)
	#r.set_f_params(*[xks_edges, Nks0, Mks0, Mks_dry0, 
	#				 V, aerosol.kappa, aerosol_rho, nk, dt, log])

	Nks = Nks0[:]
	Mks = Mks0[:]
	Mks_dry = Mks_dry0[:]

	print t0, y0

	t0s = [t0, ]
	y0s = [y0[:6], ]
	Nks_all, Mks_all = [Nks0, ], [Mks0, ]
	#while r.successful() and r.t < t_end:
	#	r.integrate(r.t + dt)
	#	print r.t, r.y[:6]
	#	y0s.append(r.y[:6])
	#	t0s.append(r.t)
	t = t0
	y = y0
	while t < t_end:
		time_points = [t, t+dt]
		solver.set_initial_condition(y)
		solver.set(f_args=[xks_edges, Nks, Mks, Mks_dry, 
			 			   V, aerosol.kappa, aerosol_rho, nk, dt, log])
		y, _ = solver.solve(time_points)
		y = y[-1, :]

		t += dt

		y0s.append(y[:6])
		t0s.append(t)

		#dms = r.y[6:nk+6]*1.0
		#dms = y[6:nk+6]*1.0
		dms_low = r.y[6:nk+6]*1.0
		dms_high = r.y[nk+6:2*nk+6]*1.0

		Nks, Mks, Mks_dry = bin_methods.adjust_bins(xks_edges, dms_low, dms_high,
													Nks, Mks, Mks_dry,
													nk, 0)
		Nks_all.append(Nks)
		Mks_all.append(Mks)
		print dms
		#for rk, N, M in zip(rks, Nks, Mks):
		#	print rk, N, M
		print "     Nk", np.ma.sum(Nks)/np.ma.sum(Nks0)
		print "     Mk", np.ma.sum(Mks)/np.ma.sum(Mks0)
		print " Mk_dry", np.ma.sum(Mks_dry)/np.ma.sum(Mks_dry0)
		#if np.ma.sum(Mks_dry)/np.ma.sum(Mks_dry0) > 1.1 or \
		#	np.ma.sum(Mks_dry)/np.ma.sum(Mks_dry0) < 0.9:
		#	raise ValueError("Too much dry mass change")

	#	r.set_f_params(*[xks_edges, Nks, Mks, Mks_dry, 
	#					 V, aerosol.kappa, aerosol_rho, nk, dt, log])

	#	r.y[6:nk+6] = 0.
	#	y[6:nk+6] = 0.
		r.y[6:2*nk+6] = 0.
		#r.y[:6] = y0[:6]

	Nks_all = np.array(Nks_all)
	Mks_all = np.array(Mks_all)
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
	fig, (ax1, ax2) = subplots(2, 1, sharex=True, num=2, figsize=(10, 4))
	ax1.plot(xks_edges[:-1], Nks0, 'b',  marker='x', label='before')
	ax1.plot(xks_edges[:-1], Nks,  'y',  marker='o', label='after')
	ax1.plot(xks_edges[:-1], Nks_dry0, 'k', marker='^', label='dry')
	ax1.semilogy()
	ax1.semilogx()
	ax_top = ax1.twiny()
	ax_top.grid()
	ax_top.set_xlim(rks_edges[0], rks_edges[-2])
	ax_top.semilogx()
	ax1.set_xlim(xks_edges[0], xks_edges[-2])
	ax1.set_ylim(1e5, 1e9)	
	ax1.set_ylabel("Number in bin (m$^-3$)")
	ax_top.set_xlabel("Droplet radius (m)")

	ax2.plot(xks_edges[:-1], Mks0, 'b',  marker='x', label='before')
	ax2.plot(xks_edges[:-1], Mks,  'y',  marker='o', label='after')
	ax2.plot(xks_edges[:-1], Mks_dry0, 'k', marker='^', label='dry')
	ax2.semilogy()
	ax2.semilogx()
	ax2.set_xlim(xks_edges[0], xks_edges[-2])
	ax2.set_ylim(1e-12, 1e-5)
	ax2.set_xlabel("Droplet Mass (kg)")
	ax2.set_ylabel("Mass in bin (kg/m$^3$)")
	ax2.legend(loc='upper left')

	ax2.grid(False, 'minor')
	ax1.grid(False, 'minor')
	ax_top.grid(False, 'major')
	ax_top.grid(False, 'minor')

	fig.subplots_adjust(hspace=0.05)
	draw()

	
	figure(3); clf();
	fig, axes = subplots(2, 3, num=3, figsize=(20,8))
	cols = { 'P': (79900, 80001),
			 'T': (283., 284.),
			 'S': (S0, pp.S.max()*1.01),
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

	from matplotlib.colors import LogNorm
	fig = figure(5); clf();
	ax = fig.add_subplot(111)
	yy = t0s
	xx = range(nk)
	xx, yy = np.meshgrid(xx, yy)
	pc = ax.pcolor(xx, yy, Nks_all, norm=LogNorm(vmin=1e-1, vmax=Nks_all.max(), clip=False),
					   cmap='PuBu', rasterized=True)
	#im = ax.imshow(Nks_all, norm=LogNorm(vmin=1e-1, vmax=Nks_all.max(), clip=True),
	#				   cmap='PuBu', interpolation='None', origin='lower')
	#im.set_extent([0, nk, 0, t0s[-1]])
	ax.grid()
	cb = fig.colorbar(pc)
	cb.set_label("Number Concentration (m$^-3$)")
	ax.set_xlim(0, nk); 
	ax.set_xlabel("Bin Number")
	ax.set_ylabel("time")

	ax_S = ax.twiny()
	ax_S.plot(pp.S, pp.index, color='k', linewidth=1.5)
	ax_S.plot(x.S, t0s, color='r', marker='x', markersize=2)
	ax_S.set_xlim(S0, 1.01*np.max([pp.S.max(), x.S.max()]))
	ax_S.set_xlabel("Supersaturation")
	ax_S.grid()
	ax.set_ylim(0, t0s[-1])

	draw()

	#from mpl_toolkits.mplot3d import Axes3D
	#fig = figure(6); clf()
	#ax = fig.add_subplot(111, projection='3d')


"""
from parcel_aux_bin import Seq
from pylab import *
ion()

import numpy as np

rs = np.logspace(np.log10(r_low), np.log10(r_high), 1000)
plot(rs, [Seq(r, 1.053431e-08, 2.831451e+02, 6.000000e-01) for r in rs])
"""

