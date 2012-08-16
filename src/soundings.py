"""
Utilities for reading and visualizing soundings
"""

import numpy as np
from pandas import *
from pylab import *

from microphysics import gamma_d, Rd, g, epsilon
from pywrfplotUtils import gamma_s, es
from scipy.interpolate import interp1d

def read_PACS(filename):
    """
    Read the ASCII data from a PACS sounding file, and save the resulting sounding in a pandas
    DataFrame. There are more than one sounding in each PACS file, but this method will split up
    individual soundings and save them separately as a .csv generated from the DataFrame which can
    easily be re-imported.
    """
    f = open(filename, "r")
    
    ## First, split out all the individual soundings from the file
    lines = f.readlines()
    sounding_lines = []
    current_sounding = []
    for i, line in enumerate(lines):
        if line.startswith("Data Type"):
            if current_sounding:
                sounding_lines.append(current_sounding)
            current_sounding = []
        current_sounding.append(line)
    f.close()  
    
    ## Process each sounding
    soundings = []
    for lines in sounding_lines:
        s = {}
        # Top metadata
        for line in lines[:5]:
            left, right = line.strip().split(": ")
            s[left] = right
        
        # Sounding data
        pressure, temp, dewpoint, rh, alt = [], [], [], [], []
        for line in lines[15:]:
            bits = line.strip().split()
            p, t, dt, r, a = bits[1], bits[2], bits[3], bits[4], bits[14]
            pressure.append(p)
            temp.append(t)
            dewpoint.append(dt)
            rh.append(r)
            alt.append(a)
        data = {}
        data['pressure'] = pressure
        data['temperature'] = temp
        data['dewpoint'] = dewpoint
        data['relative humidity'] = rh
        data['height'] = alt
        s['data'] = data
        soundings.append(s)
    
    for s in soundings:
        sounding = DataFrame(s['data'])
        site = s['Release Site Type/Site ID'].strip()
        site = site.split(" ")[0]
        date = s['UTC Release Time (y,m,d,h,m,s)'].strip()
        date = date.replace(", ", "-")
        print "data/%s_%s" % (site, date)
        sounding.to_csv("data/%s_%s.csv" % (site, date))

def read_sounding(filename):
    """Retreive a previously generated sounding made with method read_PACS, and return
    the sounding as a new DataFrame"""
    return read_csv(filename, index_col=0, na_values=['999.0', '9999.0'])

def plot_sounding(sounding, parcel_sounding=None):
    """
    Skew-T plot of sounding, adapted from pywrfplot
    """
    from pywrfplotParams import T_zero,kappa,P_bot
    from pywrfplotUtils import gamma_s
    
    skewness = 37.5
    P_b = 105000. # Bottom pressure, Pa
    P_t = 10000. # Top pressure, Pa
    dp = 100. # Pressure increment, Pa
    plevs = np.arange(P_b, P_t-1, -dp)
    
    def _skewnessTerm(P):
        return skewness * np.log(P_bot/P)

    def _isotherms():
        for temp in np.arange(-140,50,10):
            plt.semilogy(temp + _skewnessTerm(plevs), plevs,  basey=math.e, \
                         color = ('blue' if temp <= 0 else 'red'), \
                         linestyle=('solid' if temp == 0 else 'dashed'), linewidth = .5)
            
    def _isobars():
        for n in np.arange(P_bot,P_t-1,-10**4):
            plt.plot([-40,50], [n,n], color = 'black', linewidth = .5)
            
    def _dry_adiabats():
        for tk in T_zero+np.arange(-30,210,10):
            dry_adiabat = tk * (plevs/P_bot)**kappa - T_zero + _skewnessTerm(plevs)
            plt.semilogy(dry_adiabat, plevs, basey=math.e, color = 'brown', \
                         linestyle='dashed', linewidth = .5)
    
    def _moist_adiabats():
        ps = [p for p in plevs if p<=P_bot]
        for temp in np.concatenate((np.arange(-40.,10.1,5.),np.arange(12.5,45.1,2.5))):
            moist_adiabat = []
            for p in ps:
                temp -= dp*gamma_s(temp,p)
                moist_adiabat.append(temp + _skewnessTerm(p))
            plt.semilogy(moist_adiabat, ps, basey=math.e, color = 'green', \
                         linestyle = 'dotted', linewidth = .5)
        
    def _temperature(temperature, pressure):
        plt.semilogy(temperature + _skewnessTerm(pressure), pressure, basey=math.e, color = 'black', \
                     linestyle='solid', linewidth = 1.5)
    
    def _dewpoint(dewpoints, pressure):
        plt.semilogy(dewpoints + _skewnessTerm(pressure), pressure, basey=math.e, color = 'red', \
                     linestyle='solid', linewidth = 1.5)
    
    ## Start making figure
    fig = figure()
    _isotherms()
    _isobars()
    _dry_adiabats()
    _moist_adiabats()
    
    pressure = sounding['pressure']*100. # hPa -> Pa
    temperature = sounding['temperature'] # C
    dewpoint = sounding['dewpoint']
    
    _temperature(temperature, pressure)
    _dewpoint(dewpoint, pressure)
    
    if parcel_sounding is not None:
        p_pressure, p_temperature = parcel_sounding['pressure']*100., parcel_sounding['temperature']
        plt.semilogy(p_temperature + _skewnessTerm(p_pressure), p_pressure, basey=math.e, color='purple',
                     linestyle='dashed', linewidth=1.5)
    
    axis([40, 50, P_b, P_t])
    xlabel('Temperature ($^{\circ}\! C$) at 1000 hPa')
    xtick = np.arange(-40, 51, 5)
    xticks(xtick,['' if tick%10 != 0 else str(tick) for tick in xtick])
    ylabel('Pressure (hPa)')
    ytick = np.arange(P_bot, P_t-1, -10**4)
    yticks(ytick, ytick/100.)
    
def virtual_temperature(sounding):
    temperature = sounding['temperature']
    rh = sounding['relative humidity']/100. # percent -> decimal
    pressure = sounding['pressure']*100. # hPa -> Pa
    
    q = epsilon*es(temperature)*rh/pressure
    tv = (temperature+273.15)*(1.+0.61*q)
    return tv
    
def parcel_sounding(sounding):    
    trunc_sounding = sounding[sounding['height'] < 10000]
    
    ## Use Espy formula to estimate LCL height in meters
    first_ind = trunc_sounding.index[0]
    t = trunc_sounding['temperature'].ix[first_ind]
    td = trunc_sounding['dewpoint'].ix[first_ind]
    rh = trunc_sounding['relative humidity'].ix[first_ind]/100.
    
    lcl = 125.*(t - td)

    parcel_temp_interp = interp1d(trunc_sounding['height'], trunc_sounding['temperature'], 'slinear')
    parcel_rh_interp = interp1d(trunc_sounding['height'], trunc_sounding['relative humidity']/100., 'slinear')
    
    dz = 10. # meters
    z0, p0 = trunc_sounding['height'].ix[first_ind], trunc_sounding['pressure'].ix[first_ind]
    
    q = epsilon*es(t)*rh/(p0*100.)
    tv = (t+273.15)*(1.+0.61*q)
    
    zs, ps, ts, tvs = [z0], [p0], [t], [tv]
    ## Integrate the dry adiabatic lapse rate up to the LCL
    while True:
        z, p, t = zs[-1], ps[-1], ts[-1]
        #print z, p, t
        if z > lcl: break
        zs.append(z+dz)
        ts.append(t-(gamma_d*dz*1e-3))
        
        q = epsilon*es(ts[-1])*parcel_rh_interp(z+dz)/(p0*100.)        
        tvs.append((ts[-1]+273.15)*(1.+0.61*q))
        
        rho = (p*100.)/(Rd*(t+273.15)) # p: hPa -> Pa; t: C -> K
        dp = -rho*g*dz
        ps.append(p+dp/100.)
        
    ## Integrate the moist adiabatic lapse rate up to 10 km. Note that gamma_s is now in dT/dp
    dp = 500. # Pa
    while True:
        z, p, t = zs[-1], ps[-1], ts[-1]
        #print z, p, t
        
        q = epsilon*es(t)/(p*100.)
        tv = (t+273.15)*(1.+0.61*q)
        rho = (p*100.)/(Rd*tv)
        dz = dp/(rho*g)
        
        if z+dz > 1e4: break
        zs.append(z+dz)
        ps.append(p-dp/100.)
        ts.append(t-(gamma_s(t, p*100.)*dp))
        tvs.append((ts[-1]+273.15)*(1.+0.61*(epsilon*es(ts[-1])/(p*100.))))        
                
    return DataFrame({'height': zs, 'pressure': ps, 'temperature': ts, 'virtual temperature': tvs})