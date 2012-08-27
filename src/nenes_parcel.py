'''
Cloud Parcel Model based on Nenes et al, 2001
'''

import sys
import numpy as np
import pandas
from lognorm import Lognorm
from pylab import *
ion()
from scipy.optimize import fsolve, broyden1, broyden2
from scipy.integrate import ode

import nenes_parcel_aux as npa
from micro import *

class ParcelModel(object):
    
    def __init__(self, aerosol, V):
        """
        Initialize the model with information about updraft speed and the aerosol
        distribution. To keep things as general as possible, the aerosol distribution
        should be a dictionary with two keys - "r_drys" and "Nis", which are the 
        aerosol dry radii (meters) and number concentartions (per m**3) given as
        lists/arrays with the dry radii and concentrations divided into discrete
        bins.
        """
        self.aerosol = aerosol
        self.V = V
        
    def _setup_run(self, P0, T0, S0, make_plot=False):
        """
        Setup the initial parcel state for the run, given the initial
        temperature (K), pressure (Pa), and supersaturation, as well as
        the number of bins in which to divide the aerosol distribution.        
        """        
        out = dict()
        r_drys = self.aerosol['r_drys']
        Nis = self.aerosol['Nis']
        
        ## 1) Setup aerosols
        # a - aerosol droplet solutes
        d_drys = 2.*r_drys
        nss = np.array([(rho_p*np.pi*epsilon*(d_dry**3.))/(6.*Ms) for d_dry in d_drys])
        out['nss'] = nss
        print "AEROSOL DISTRIBUTION"
        print "%8s %6s" % ("r", "N")
        for r, N in zip(r_drys, Nis):
            print "%2.2e %4.1f" % (r, N)
        print "\n"+"-"*44
            
        ## 2) Setup parcel initial conditions
        # a) water vapor
        wv0 = (1.-S0)*0.622*es(T0-273.15)/(P0-es(T0-273.15)) # Water Vapor mixing ratio, kg/kg
        
        # b) find equilibrium wet particle radius
        # need to seed the solver with a good guess
        r_guesses = []
        print "calculating seeds for equilibrium solver..."
        r_guesses = npa.guesses(T0, S0, r_drys, nss)
        print " done"
        '''
        for i in xrange(len(r_drys)):
            print i
            rdi = r_drys[i]
            #xi = np.arange(rdi+1e-9, 1e-5, 1e-9)
            xi = np.linspace(rdi+1e-9, 1e-5, 10000.)
            yi = np.array([Seq(T0, r, rdi, nss[i]) for r in xi])
            idx_min = (np.abs(yi - S0)).argmin()
            r_guesses.append(xi[idx_min])
        '''
                
        f = lambda r, r_dry, ns: (Seq(T0, r, r_dry, ns) - S0)
        r0s = np.array([fsolve(f, guess, args=(r, ns), xtol=1e-10)[0] for guess, r, ns in zip(r_guesses, r_drys, nss)])
        #print r0s
        print "equilib r test: ", r0s[0], Seq(T0, r0s[0], r_drys[0], nss[0]) ,"\n"
        for i, (r, r_dry, ns) in enumerate(zip(r0s, r_drys, nss)):
            ss = Seq(T0, r, r_dry, ns)
            #print r, ss
            if r < 0 or r > 1e-3: 
                print "Found bad r", r, r_dry
                #r0s[i] = fsolve(f, r_dry+1e-9, args=(r_dry, ns), xtol=1e-10)[0]
            if np.abs(ss-S0)/S0 > 0.02: print "Found S discrepancy", ss, r_dry
        #print r0s
        
        #while np.any(r0s < 0):
        #    r0s = np.sign(r0s)*r0s
        #    r0s_new = np.array([fsolve(f, r0, args=(r, ns), xtol=1e-20)[0] for r0, r, ns in zip(r0s, r_drys, nss)])
        #    print r0s_new[0], Seq(T0, r0s_new[0], r_drys[0], nss[0])
        #    r0s = r0s_new[:]

        #fig = figure(2)
        #x = np.arange(r_drys[0]+1e-9, 1e-5, 1e-9)
        #y = np.array([Seq(T0, r, r_drys[0], nss[0]) for r in x])
        #plot(x, y)
        ##for xi, yi in zip(x,y): print xi, yi
        #semilogx(); ylim(-.2, 0.1)
        
        # c) compute equilibrium droplet water content
        wc0 = np.sum([(4.*np.pi/3.)*rho_w*Ni*(r0**3 - r_dry**3) for r0, r_dry, Ni in zip(r0s, r_drys, Nis)])
        
        # d) concat into initial conditions arrays
        y0 = [P0, T0, wv0, wc0, S0]
        print "PARCEL INITIAL CONDITIONS"
        print "    {:>8} {:>8} {:>8} {:>8} {:>8}".format("P (hPa)", "T (K)", "wv", "wc", "S")
        print "      %4.1f   %3.2f  %3.1e   %3.1e   %01.2f" % (P0/100., T0, wv0, wc0, S0)
        y0.extend(r0s)
        y0 = np.array(y0)
        out['y0'] = y0
        
        return out
        
    def run(self, P0, T0, S0, z_top, dt=0.1, max_steps=20000, integrator="vode"):
        
        setup_results = self._setup_run(P0, T0, S0)
        nss = setup_results['nss']
        y0 = setup_results['y0']
        r_drys = self.aerosol['r_drys']
        Nis = self.aerosol['Nis']
        nr = len(r_drys)
        
        ## Setup run time conditions        
        t0 = 0.
        t_end = z_top/self.V
        t = np.arange(t0, t_end+dt, dt)
        print "\n"+"n_steps = %d" % (len(t))+"\n"
        
        raw_input("Continue run?")
        
        ## Setup integrator
        if integrator == "vode":
            r = ode(npa.der).set_integrator("vode", method="bdf", nsteps=max_steps, order=5)
        else:
            r = ode(npa.der).set_integrator("dop853", nsteps=max_steps)
        r.set_initial_value(y0)
        
        nss, r_drys, Nis = np.array(nss), np.array(r_drys), np.array(Nis)
        r.set_f_params(nr, nss, r_drys, Nis, self.V)
    
        out = []
        while r.successful() and r.t < t_end:
            print "1", r.t, r.t*self.V, r.y[:5]
            r.integrate(r.t+dt)
            out.append(r.y)
            if not r.successful: print "failure"; exit()
        x = np.array(out)
    
        heights = t*self.V
        print len(heights), x.shape
        offset = 0
        if len(heights) > x.shape[0]:
            offset = 1
        
        df1 = pandas.DataFrame( {'P':x[:,0], 'T':x[:,1], 'wv':x[:,2],
                                'wc':x[:,3], 'S':x[:,4]} , index=heights[offset:])
        df1.to_csv("output_parcel.csv", index_label="height")
        
        labels = ["r%03d" % i for i in xrange(nr)]
        radii_dict = dict()
        for i, label in enumerate(labels):
            radii_dict[label] = x[:,5+i]
        df2 = pandas.DataFrame( radii_dict, index=heights[offset:])
        df2.to_csv("output_aerosol.csv", index_label="height")
        
        return df1, df2

if __name__ == "__main__":
    
    ## Initial conditions
    P0 = 100000. # Pressure, Pa
    T0 = 298.15 # Temperature, K
    S0 = -0.02 # Supersaturation. 1-RH from wv term
    V = 0.1 # m/s
    
    ## Aerosol properties
    ## RS SHOULD BE MONOTONICALLY INCREASING!!!!!!
    mu, sigma, N, bins = 0.1, 2.0, 100., 200
    l = 0
    r = bins
    aerosol_dist = Lognorm(mu=mu, sigma=sigma, N=N)
    rs = np.logspace(-2., 0., num=bins+1)[:]
    mids = np.array([np.sqrt(a*b) for a, b in zip(rs[:-1], rs[1:])])[l:r]
    Nis = np.array([0.5*(b-a)*(aerosol_dist.pdf(a) + aerosol_dist.pdf(b)) for a, b in zip(rs[:-1], rs[1:])])[l:r]
    r_drys = mids*1e-6
    ##r_drys = np.array([0.01e-6])
    ##Nis = np.array([100.])
    figure(10)
    clf()
    bar(rs[:-1], Nis, diff(rs))
    vlines(mids, 0, ylim()[1], color='red', linestyle='dotted')
    semilogx()
    Nis = Nis*1e6
    Ntot = np.sum(Nis)
    aerosol = {"r_drys":r_drys, "Nis":Nis}    
    pm = ParcelModel(aerosol, V)
    
    ## Run model
    
    dt = np.max([V/100., 0.01])
    dt = 0.005
    parcel, aerosols = pm.run(P0, T0, S0, z_top=50.0, 
                              dt=dt, max_steps=6000, integrator="vode")
    
    figure(1)
    p = parcel.S.plot(logx=False)
    max_idx = np.argmax(parcel.S)
    max_z = parcel.index[max_idx]
    vlines([max_z], ylim()[0], ylim()[1], color='k', linestyle='dashed')
    #figure(2)
    #p2 = aerosols.r000.plot(logx=False)    
    
    #Smax = parcel.S.max()
    def s_crit(d_s, T):
        #A = (0.66/T)*1e-6
        A = (4.*Mw*sigma_w(T))/(R*T*rho_w)        
        B = (4.*(A**3.)*rho_w*Ms)/(27.*nu*rho_p*Mw*(d_s**3))
        return np.exp(B**0.5) - 1.
    
    def d_crit(d_s, T):
        A = (4.*Mw*sigma_w(T))/(R*T*rho_w)
        ns = (rho_p*np.pi*epsilon*(d_s**3.))/(6.*Ms) 
        B = (6.*ns*Mw)/(np.pi*rho_w)
        return np.sqrt(3.*B/A)
    '''
    Neq = []
    Nkn = []
    Nunact = []
    S_max = S0
    print "N analysis"
    for S, T, i in zip(parcel.S, parcel['T'], xrange(len(parcel.S))):
        print parcel.index[i]
        if S > S_max: S_max = S
        Neq_step = 0.
        s_crits = [s_crit(2.*r, T) for r in r_drys]
        for sc, Ni in zip(s_crits, Nis):
            if S_max > sc: Neq_step += Ni
        Neq.append(Neq_step)
        
        Nkn_step = 0.
        r_crits = [d_crit(2.*r, T) for r in r_drys]
        rs = np.array(aerosols.ix[i])
        min_found = False
        for j, (r, rc) in enumerate(zip(rs, r_crits)):
            if r > rc: min_found = True
            if min_found: Nkn_step += Nis[j]
        Nkn.append(Nkn_step)
        Nunact.append(Ntot - Nkn_step)    
            
    parcel['Neq'] = Neq
    parcel['Nkn'] = Nkn
    parcel['Nunact'] = Nunact
    
    alphaz = parcel.Nkn/parcel.Neq
    alphaz[isnan(alphaz)] = 0.
    phiz = Nunact/parcel.Nkn
    phiz[phiz == inf] = 1.
    
    parcel['alpha'] = alphaz
    parcel['phi'] = phiz
    '''
