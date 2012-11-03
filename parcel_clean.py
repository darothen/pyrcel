"""Cloud Parcel Model based on Nenes et al, 2001; implemented by 
Daniel Rothenberg (darothen@mit.edu) as part of research undertaken for the
General examination in the Program in Atmospheres, Oceans, and Climate TEST"""

__docformat__ = 'reStructuredText'
import pandas
from lognorm import Lognorm

from pylab import *
ion()

from scipy.optimize import fsolve
from scipy.integrate import odeint

from nenes_parcel_aux import der, guesses
from micro import Seq, es, rho_w, kohler_crit

class AerosolSpecies(object):
    """This class organizes metadata about an aerosol species"""
    
    def __init__(self, **values):
        self.species = values['species'] # Species molecular formula
        self.kappa = values['kappa'] # Kappa hygroscopicity parameter
        
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
    
    def __init__(self, aerosols, V, T0, S0, P0, console=False, write_output=False):
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
        self.write_output = write_output
        
    def _setup_run(self, P0, T0, S0, make_plot=False):
        """
        Setup the initial parcel state for the run, given the initial
        temperature (K), pressure (Pa), and supersaturation, as well as
        the number of bins in which to divide the aerosol distribution.        
        """        
        out = dict()
        
        ## 1) Setup aerosols
        # a) grab all the initial aerosol size/concentrations
        species = []
        r_drys, Nis, kappas = [], [], []
        for aerosol in self.aerosols:
            r_drys.extend(aerosol.r_drys)
            kappas.extend([aerosol.kappa]*aerosol.nr)
            Nis.extend(aerosol.Nis)
            species.extend([aerosol.species]*aerosol.nr)
        
        r_drys = np.array(r_drys)
        kappas = np.array(kappas)
        Nis = np.array(Nis)
        
        out['r_drys'] = r_drys
        out['kappas'] = kappas
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
        if self.console: print "calculating seeds for equilibrium solver..."
        r_guesses = np.array(guesses(T0, S0, r_drys, kappas))
        #r_guesses = []
        #for aerosol in self.aerosols:
        #    r_guesses.extend(guesses(T0, S0, r_drys, kappas))
        #r_guesses = np.array(r_guesses)
        if self.console: print " done"
                            
        # wrapper function for quickly computing deviation from chosen equilibrium supersaturation given
        # current size and dry size
        f = lambda r, r_dry, kappa: (Seq(r, r_dry, T0, kappa) - S0)
        ## Compute the equilibrium wet particle radii
        r0s = np.array([fsolve(f, (rd+guess)/2., args=(rd, kappa), xtol=1e-10)[0] for guess, rd, kappa in zip(r_guesses, r_drys, kappas)])
        ## Console logging output, if requested, of the equilibrium calcuations. Useful for 
        ## checking if the computations worked
        if self.console: 
            for (r,  r_dry, sp, kappa) in zip(r0s, r_drys, species, kappas):
                ss = Seq(r, r_dry, T0, kappa)
                rc, _ = kohler_crit(T0, r_dry, kappa)
                if r < 0 or r > 1e-3: print "Found bad r", r, r_dry, sp
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
        y0 = setup_results['y0']
        r_drys = setup_results['r_drys']
        kappas = setup_results['kappas']
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
            x, info = odeint(der, y0, t, args=(nr, r_drys, Nis, self.V, kappas),
                             full_output=1, printmessg=1, ixpr=1, mxstep=max_steps,
                             mxhnil=0, atol=1e-15, rtol=1e-12)
        else:
            x = odeint(der, y0, t, args=(nr, r_drys, Nis, self.V, kappas),
                       mxstep=max_steps, mxhnil=0, atol=1e-15, rtol=1e-12)
    
        heights = t*self.V
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
        
        if self.write_output:
            df1.to_csv("parcel_trace.csv")
            for aerosol in self.aerosols:
                species = aerosol.species
                df = aerosol_dfs[species]
                df.to_csv("%s.csv" % species)
        
        return df1, pandas.Panel(aerosol_dfs)

if __name__ == "__main__":
    
    ## Initial conditions
    P0 = 95000. # Pressure, Pa
    T0 = 285.2 # Temperature, K
    S0 = -0.05 # Supersaturation. 1-RH from wv term
    V = 0.5 # m/s
    
    ## Aerosol properties
    ## AEROSOL 1 - (NH4)2SO4
    # Chemistry
    ammonium_sulfate = { 
        'kappa': 0.6, # Hygroscopicity parameter
    }
    # Size Distribution
    mu, sigma, N, bins = 0.05, 2.0, 300., 200
    l = 0
    r = bins
    aerosol_dist = Lognorm(mu=mu, sigma=sigma, N=N)
    lr, rr = np.log10(mu/(10.*sigma)), np.log10(mu*10.*sigma)
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
    NaCl = { 'kappa': 0.01, # Hygroscopicity parameter
    }
    # Size Distribution
    r_drys = np.array([0.25*1e-6, ])
    rs = np.array([r_drys[0]*0.9, r_drys[0]*1.1])*1e6
    Nis = np.array([10000., ])
    
    
    NaCl['distribution'] = aerosol_dist
    NaCl['r_drys'] = r_drys
    NaCl['rs'] = rs
    NaCl['Nis'] = Nis*1e6
    NaCl['species'] = 'NaCl'
    NaCl['nr'] = len(r_drys)
    
    ######
 
    initial_aerosols = [AerosolSpecies(**ammonium_sulfate), AerosolSpecies(**NaCl)]
    print initial_aerosols
    
    aer_species = [a.species for a in initial_aerosols]
    aer_dict = dict()
    for aerosol in initial_aerosols:
        aer_dict[aerosol.species] = aerosol
    
    for aerosol in initial_aerosols: print np.sum(aerosol.Nis) *1e-6
        
    pm = ParcelModel(initial_aerosols, V, T0, S0, P0, console=True)
    
    ## Run model    
    dt = 0.01
    parcel, aerosols = pm.run(P0, T0, S0, z_top=200.0, dt=dt, max_steps=500)

