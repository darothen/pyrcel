"""Cloud Parcel Model based on Nenes et al, 2001; implemented by 
Daniel Rothenberg (darothen@mit.edu) as part of research undertaken for the
General examination in the Program in Atmospheres, Oceans, and Climate"""

__docformat__ = 'reStructuredText'
import pandas
from lognorm import Lognorm, MultiModeLognorm

from pylab import *
ion()

from scipy.optimize import fsolve
from scipy.integrate import odeint

import nenes_parcel_aux as npa
from micro import Seq, es, rho_w, Mw, sigma_w, R, kohler_crit

class AerosolSpecies(object):
    """This class organizes metadata about an aerosol species"""
    
    def __init__(self, **values):
        self.species = values['species'] # Species molecular formula
        self.Ms = values['Ms'] # Molecular weight, kg/mol
        self.epsilon = values['epsilon'] # mass fraction of soluble material
        self.nu = values['nu'] # number of ions into which a solute dissolves
        
        if not "rho_p" in values:
            self.rho_s = values['rho_s']
            self.rho_u = values['rho_u']
            self.rho_p = self.rho_u/(1.-self.epsilon*(1.-(self.rho_u/self.rho_s)))
        else:
            self.rho_p = values['rho_p']
        
        self.distribution = values['distribution']
        self.r_drys = values['r_drys']
        self.Nis = values['Nis']
        self.nr = len(self.r_drys)
        self.rs = values['rs']
        
    def __repr__(self):
        return "%s - %r" % (self.species, self.distribution)

class ParcelModel(object):
    """This class wraps logic for initializing and running the Nenes et al, 2001
    cloud parcel model
    
    The model here is implemented in order to give an object-oriented approach to
    running the model. Instead of hard-coding aerosol size distributions and parcel
    initial conditions, the user can independently setup these parameters and pass
    them to the model to calculate everything necessary for running the model.
    
    :ivar aerosols: 
    :ivar V: Parcel updraft velocity (m/s)
    :ivar T0: Parcel initial temperature (K)
    :ivar S0: Parcel initial supersaturation, relative to 100% RH 
    :ivar P0: Parcel initial pressure (Pa)        
    """
    
    def __init__(self, aerosols, V, T0, S0, P0, console=False):
        """
        Initialize the model with information about updraft speed and the aerosol
        distribution. 
        """
        self.aerosols = aerosols
        self.V = V
        self.T0 = T0
        self.S0 = S0
        self.P0 = P0
        
        self.console = console
        
    def _setup_run(self, P0, T0, S0, make_plot=False):
        """
        Setup the initial parcel state for the run, given the initial
        temperature (K), pressure (Pa), and supersaturation, as well as
        the number of bins in which to divide the aerosol distribution.        
        """        
        out = dict()
        
        ## 1) Setup aerosols
        # a - aerosol droplet solutes
        nss_chunks, species = [], []
        r_drys, Nis = [], []
        for aerosol in self.aerosols:
            rho_p, epsilon, Ms = aerosol.rho_p, aerosol.epsilon, aerosol.Ms
            r_drys.extend(aerosol.r_drys)
            Nis.extend(aerosol.Nis)
            d_drys = 2.*aerosol.r_drys
            species.extend([aerosol.species]*len(d_drys))
            nss_chunks.extend([(rho_p*np.pi*epsilon*(d_dry**3.))/(6.*Ms) for d_dry in d_drys])
        
        nss = np.array(nss_chunks)
        r_drys = np.array(r_drys)
        Nis = np.array(Nis)
        
        out['nss'] = nss
        out['r_drys'] = r_drys
        out['Nis'] = Nis
        
        if self.console:
            print "AEROSOL DISTRIBUTION"
            print "%8s %6s" % ("r", "N")
            for sp, r, N in zip(species, r_drys, Nis):
                print "%10s %2.2e %4.1f" % (sp, r, N)
            print "\n"+"-"*44
            
        ## 2) Setup parcel initial conditions
        # a) water vapor
        wv0 = (1.-S0)*0.622*es(T0-273.15)/(P0-es(T0-273.15)) # Water Vapor mixing ratio, kg/kg
        
        # b) find equilibrium wet particle radius
        if self.console:
            print "calculating seeds for equilibrium solver..."
            r_guesses = npa.guesses(T0, S0, r_drys, aerosol.epsilon, aerosol.rho_p, aerosol.Ms, aerosol.nu) #@UndefinedVariable
            print " done"
        else:
            r_guesses = npa.guesses(T0, S0, r_drys, aerosol.epsilon, aerosol.rho_p, aerosol.Ms, aerosol.nu) #@UndefinedVariable
                
        f = lambda r, r_dry: (Seq(T0, r, r_dry, aerosol.epsilon, aerosol.rho_p, aerosol.Ms, aerosol.nu) - S0)
        r0s = np.array([fsolve(f, (rd+guess)/2., args=(rd, ), xtol=1e-10)[0] for guess, rd in zip(r_guesses, r_drys)])
        if self.console: 
            #print "equilib r test: ", r0s[0], Seq(T0, r0s[0], r_drys[0], nss[0]) ,"\n"
            for (r, rg, r_dry) in zip(r0s, r_guesses, r_drys,):
                ss = Seq(T0, r, r_dry, aerosol.epsilon, aerosol.rho_p, aerosol.Ms, aerosol.nu)
                rc, sc = kohler_crit(T0, r_dry, aerosol.epsilon, aerosol.rho_p, aerosol.Ms, aerosol.nu)
                #print r, rg, r_dry, ss, rc, sc
                if r < 0 or r > 1e-3: print "Found bad r", r, r_dry
                if np.abs(ss-S0)/S0 > 0.02: print "Found S discrepancy", ss, r_dry
        out['r0s'] = r0s

        # c) compute equilibrium droplet water content
        wc0 = np.sum([(4.*np.pi/3.)*rho_w*Ni*(r0**3 - r_dry**3) for r0, r_dry, Ni in zip(r0s, r_drys, Nis)])
        
        # d) concat into initial conditions arrays
        y0 = [P0, T0, wv0, wc0, S0]
        if self.console:
            print "PARCEL INITIAL CONDITIONS"
            print "    {:>8} {:>8} {:>8} {:>8} {:>8}".format("P (hPa)", "T (K)", "wv", "wc", "S")
            print "      %4.1f   %3.2f  %3.1e   %3.1e   %01.2f" % (P0/100., T0, wv0, wc0, S0)
        y0.extend(r0s)
        y0 = np.array(y0)
        out['y0'] = y0
        self.y0 = y0
        
        return out
        
    def run(self, P0, T0, S0, z_top, dt=0.1, max_steps=1000):
        
        setup_results = self._setup_run(P0, T0, S0)
        nss = setup_results['nss']
        y0 = setup_results['y0']
        r_drys = setup_results['r_drys']
        Nis = setup_results['Nis']
        nr = len(r_drys)
        
        aerosol = self.aerosols[0]
        
        ## Setup run time conditions        
        t0 = 0.
        if self.V:
            t_end = z_top/self.V
        else:
            t_end = dt*1000
        t = np.arange(t0, t_end+dt, dt)
        if self.console:
            print "\n"+"n_steps = %d" % (len(t))+"\n"
            raw_input("Continue run?")
        
        ## Setup integrator
        if self.console:
            x, info = odeint(npa.der, y0, t, args=(nr, r_drys, Nis, self.V, 
                                                   aerosol.epsilon, aerosol.rho_p, aerosol.Ms, aerosol.nu),
                             full_output=1, printmessg=1, ixpr=1, mxstep=max_steps,
                             mxhnil=0, atol=1e-15, rtol=1e-12)
        else:
            x = odeint(npa.der, y0, t, args=(nr, r_drys, Nis, self.V,  
                                             aerosol.epsilon, aerosol.rho_p, aerosol.Ms, aerosol.nu),
                             mxstep=max_steps,
                             mxhnil=0, atol=1e-15, rtol=1e-12)
        #s = info['message']
        #l = len(s)
        #print "#"*(l+2)
        #print "#"+s+"#"
        #print "#"*(l+2)
    
        heights = t*self.V
        #print len(heights), x.shape
        offset = 0
        if len(heights) > x.shape[0]:
            offset = 1
        
        df1 = pandas.DataFrame( {'P':x[:,0], 'T':x[:,1], 'wv':x[:,2],
                                'wc':x[:,3], 'S':x[:,4]} , index=heights[offset:])
        
        aerosol_dfs = {}
        species_shift = 0 # increment by nr to select the next aerosol's radii
        for aerosol in self.aerosols:
            nr = aerosol.nr
            species = aerosol.species
            
            labels = ["r%03d" % i for i in xrange(nr)]
            radii_dict = dict()
            for i, label in enumerate(labels):
                radii_dict[label] = x[:,5+species_shift+i]
                
            aerosol_dfs[species] = pandas.DataFrame( radii_dict, index=heights[offset:])
            species_shift += nr
        
        
        return df1, pandas.Panel(aerosol_dfs)

if __name__ == "__main__":
    
    ## Initial conditions
    P0 = 100000. # Pressure, Pa
    T0 = 294. # Temperature, K
    S0 = -0.02 # Supersaturation. 1-RH from wv term
    V = 3.0 # m/s
    
    ## Aerosol properties
    ## AEROSOL 1 - (NH4)2SO4
    # Chemistry
    ammonium_sulfate = { 'Ms': 0.13214, # Molecular weight, kg/mol
                         'rho_s': 1.769*1e-3*1e6,
                         'rho_u': 1.700*1e-3*1e6,
                         'nu': 3.0, # number of ions into which a solute dissolve
                         'epsilon': 0.3 # mass fraction of soluble material in the dry particle
                       }
    #ammonium_sulfate['rho_p'] = ammonium_sulfate['rho_u']/(1.-ammonium_sulfate['epsilon']*(1.-(ammonium_sulfate['rho_u']/ammonium_sulfate['rho_s']))) # total wet particle density
    #ammonium_sulfate['rho_p'] = 1.769*1e-3*1e6

    # Size Distribution
    mu, sigma, N, bins = 0.05, 2.0, 10000., 200
    l = 0
    r = bins
    aerosol_dist = Lognorm(mu=mu, sigma=sigma, N=N)
    #rs = np.logspace(-2.5, 1.5, num=bins+1)[:]
    lr, rr = np.log10(mu/(10.*sigma)), np.log10(mu*10.*sigma)
    #lr, rr = np.log10(1e-3), np.log10(5.0)
    rs = np.logspace(lr, rr, num=bins+1)[:]
    mids = np.array([np.sqrt(a*b) for a, b in zip(rs[:-1], rs[1:])])[l:r]
    Nis = np.array([0.5*(b-a)*(aerosol_dist.pdf(a) + aerosol_dist.pdf(b)) for a, b in zip(rs[:-1], rs[1:])])[l:r]
    r_drys = mids*1e-6
    
    ammonium_sulfate['distribution'] = aerosol_dist
    ammonium_sulfate['r_drys'] = r_drys
    ammonium_sulfate['rs'] = rs
    ammonium_sulfate['Nis'] = Nis*1e6
    ammonium_sulfate['species'] = '(NH4)2SO4'
    ammonium_sulfate['nr'] = len(r_drys)
    
    ## AEROSOL 2 - NaCl
    # Chemistry
    NaCl = { 'Ms': 0.06, # Molecular weight, kg/mol
                         'rho_s': 1.001*1e-3*1e6,
                         'rho_u': 1.001*1e-3*1e6,
                         'nu': 2.0, # number of ions into which a solute dissolve
                         'epsilon': 0.05 # mass fraction of soluble material in the dry particle
                       }
    #NaCl['rho_p'] = NaCl['rho_u']/(1.-NaCl['epsilon']*(1.-(NaCl['rho_u']/NaCl['rho_s']))) # total wet particle density
    NaCl['rho_p'] = 2.16*1e-3*1e6

    # Size Distribution
    '''
    mu, sigma, N, bins = 0.1, 2., 10., 11
    l = 0
    r = bins
    aerosol_dist = Lognorm(mu=mu, sigma=sigma, N=N)
    #rs = np.logspace(-2.5, -1.6, num=bins+1)[:]
    lr, rr = np.log10(mu/(10.*sigma)), np.log10(mu*10.*sigma)
    rs = np.logspace(lr, rr, num=bins+1)[:]
    mids = np.array([np.sqrt(a*b) for a, b in zip(rs[:-1], rs[1:])])[l:r]
    Nis = np.array([0.5*(b-a)*(aerosol_dist.pdf(a) + aerosol_dist.pdf(b)) for a, b in zip(rs[:-1], rs[1:])])[l:r]
    r_drys = mids*1e-6
    '''
    r_drys = np.array([10.5*1e-6, ])
    rs = np.array([r_drys[0]*0.9, r_drys[0]*1.1])*1e6
    Nis = np.array([10, ])
    
    
    NaCl['distribution'] = aerosol_dist
    NaCl['r_drys'] = r_drys
    NaCl['rs'] = rs
    NaCl['Nis'] = Nis*1e6
    NaCl['species'] = 'NaCl'
    NaCl['nr'] = len(r_drys)
    
    ######
    
    #initial_aerosols = [AerosolSpecies(**ammonium_sulfate), AerosolSpecies(**NaCl)]
    #initial_aerosols = [AerosolSpecies(**NaCl)]
    initial_aerosols = [AerosolSpecies(**ammonium_sulfate)]
    print initial_aerosols
    
    for aerosol in initial_aerosols: print np.sum(aerosol.Nis) *1e-6
    
    figure(2, figsize=(12,10))
    #clf()
    subplot(3,2,1)
    colors = 'bgrcmyk'
    for i, aerosol in enumerate(initial_aerosols):
        rs, Nis = aerosol.rs, aerosol.Nis
        bar(rs[:-1], Nis/np.diff(rs)*1e-6, diff(rs), color=colors[i], alpha=0.5)
        #semilogy()
    #vlines(mids, 0, ylim()[1], color='red', linestyle='dotted')
    semilogx()
        
        
    pm = ParcelModel(initial_aerosols, V, T0, S0, P0, console=True)
    
    ## Run model    
    dt = np.max([V/100., 0.01])
    dt = 0.05
    parcel, aerosols = pm.run(P0, T0, S0, z_top=200.0, 
                              dt=dt, max_steps=500)
    
    xs = np.arange(501)
    parcel = parcel.ix[parcel.index % 1 == 0]
    aero_subset = {}
    for key in aerosols:
        aerosol = aerosols[key]
        subset = aerosol.ix[aerosol.index % 1 == 0]
        aero_subset[key] = subset
    aerosols = pandas.Panel(aero_subset)
    
    
    subplot(3,2,4)
    p = parcel.S.plot(logx=False)
    max_idx = np.argmax(parcel.S)
    max_z = parcel.index[max_idx]
    vlines([max_z], ylim()[0], ylim()[1], color='k', linestyle='dashed')
    xlabel("Height"); ylabel("Supersaturation")
    
    #subplot(3,2,2)
    #step = 1 if len(aerosols) < 20000 else 50
    #for r in aerosols[::step]:
    #    (aerosols[r]*1e6).plot(logy=True)
    #vlines([max_z], ylim()[0], ylim()[1], color='k', linestyle='dashed')
    #xlabel("Height"); ylabel("Wet Radius, micron")
        
    
    subplot(3,2,3)
    plot(parcel['T'], parcel['P']*1e-2)
    ylim(ylim()[::-1])
    hlines([parcel.P[max_idx]*1e-2], xlim()[0], xlim()[1], color='k', linestyle='dashed')
    xlabel("Temperature"); ylabel("Pressure (hPa)")
    
    ## PLOT AEROSOLS!!
    
    n_species = len(initial_aerosols)
    fig = figure(1, figsize=(9, 5*n_species))
    clf()
    for n, key in enumerate(aerosols):
        subplot(2, n_species, n+1)
        aerosol = aerosols[key]
        for r in aerosol:
            (aerosol[r]*1e6).plot(logy=True)
        vlines([max_z], ylim()[0], ylim()[1], color='k', linestyle='dashed')
        xlabel("Height"); ylabel("Wet Radius, micron")
        title(key)
    
        subplot(2, n_species, n+1 + n_species)
        initial_aerosol = initial_aerosols[n]
        rs, Nis = initial_aerosol.rs, initial_aerosol.Nis
        bar(rs[:-1], Nis, diff(rs), color=colors[n], alpha=0.2)
        
        rs, Nis = aerosol.ix[-1]*1e6, initial_aerosol.Nis
        plot(rs, Nis, color='k', alpha=0.5)
        semilogx(); #semilogy()
        
        rs, Nis = aerosol.ix[0]*1e6, initial_aerosol.Nis
        plot(rs, Nis, color='r', alpha=0.5)
        
        #rs, Nis = pm.y0[-initial_aerosol.nr:]*1e6, initial_aerosol.Nis
        #plot(rs, Nis, color='b', alpha=0.5)
        
    
    def sd_crits(d_s, T, aerosol):
        #A = (0.66/T)*1e-6
        A = (4.*Mw*sigma_w(T))/(R*T*rho_w)  
        ns = (aerosol.rho_p*np.pi*aerosol.epsilon*(d_s**3.))/(6.*aerosol.Ms) 
        B = (6.*ns*Mw)/(np.pi*rho_w)
        s_crit = np.exp(np.sqrt(4.*(A**3)/(27.*B))) - 1.0
        d_crit = np.sqrt(3.*B/A)
        return s_crit, d_crit


    Neq = []
    Nkn = []
    Nunact = []    
    S_max = S0
    aerosol = initial_aerosols[0]
    '''
    s_crits, d_crits = zip(*[sd_crits(2.*r, T0, aerosol) for r in aerosol.r_drys])
    s_crits = np.array(s_crits)
    d_drys = 2.*aerosol.r_drys
    nss = np.array([(aerosol.rho_p*np.pi*aerosol.epsilon*(d_dry**3.))/(6.*aerosol.Ms) for d_dry in d_drys])
    #r_crits = [npa.guesses(T0, s_crit, np.array([r_dry]), np.array([ns]))[0] for
    #           s_crit, r_dry, ns in zip(s_crits, aerosol.r_drys, nss)]
    r_crits = np.array(d_crits)/2.
    for a, b, s in zip(r_drys, r_crits, s_crits):
        print a, b, a > b, s
    '''

    raw_input("N analysis...")
    aerosols = aerosols[aerosol.species]
    for S, T, i in zip(parcel.S, parcel['T'], xrange(len(parcel.S))):
        '''
        print parcel.index[i],
        s_crits, d_crits = zip(*[sd_crits(2.*r, T, aerosol) for r in aerosol.r_drys])
        s_crits = np.array(s_crits)
        r_crits = [npa.guesses(T, s_crit, np.array([r_dry]), np.array([ns]))[0] for
               s_crit, r_dry, ns in zip(s_crits, aerosol.r_drys, nss)]
        '''
        r_crits, s_crits = zip(*[kohler_crit(T, r_dry, aerosol.epsilon, aerosol.rho_p, aerosol.Ms, aerosol.nu) for r_dry in aerosol.r_drys])
        s_crits = np.array(s_crits)
        r_crits = np.array(r_crits)
        if S > S_max: S_max = S

        big_s =  S_max >= s_crits
        Neq.append(np.sum(Nis[big_s]))

        rstep = np.array(aerosols.ix[i])
        #active_radii = (S > s_crits) & (rstep > r_crits)
        active_radii = (rstep > r_crits)
        #sar = np.min(active_radii) if len(active_radii) > 0 else 1e99
        if len(active_radii) > 0:
            Nkn.append(np.sum(Nis[active_radii]))
            Nunact.append(np.sum(Nis[(rstep < r_crits)]))
        else:
            Nkn.append(0.0)           
            Nunact.append(np.sum(Nis))
        
        '''    
        unactivated = Nis[rstep < r_crits]
        #unactivated = Nis[(rstep < r_crits) & (S > s_crits)]
        if len(unactivated) > 0:
            Nunact.append(np.sum(unactivated))
        else:
            Nunact.append(0.0)
        '''
        
        print parcel.index[i], Neq[i], Nkn[i], Nunact[i], S_max, S


    parcel['Neq'] = Neq
    parcel['Nkn'] = Nkn
    parcel['Nunact'] = Nunact
    
    alphaz = parcel.Nkn/parcel.Neq
    alphaz[isnan(alphaz)] = 0.
    phiz = Nunact/parcel.Nkn
    phiz[phiz == inf] = 1.
    
    parcel['alpha'] = alphaz
    parcel['phi'] = phiz
    
    figure(2)
    ax = subplot(3,2,5)
    parcel[['Neq', 'Nkn']].plot(ax=ax, grid=True)
    xlabel("Height")
    
    subplot(3,2,6)
    parcel.alpha.plot()
    parcel.phi.plot()
    ylim(0, 1)
    xlabel("Height"); ylabel(r'$\alpha(z),\quad\phi(z)$')
    print parcel.alpha.ix[-1]
