'''
Created on Jul 31, 2012

@author: daniel
'''

from pylab import *
import pylab
from microphysics import *
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


press_levels = [1000.,987.5,975.,962.5,950.,937.5,925., 
                 900.,875.,850.,825.,800.,750.,700.,650.,  
                 600.,550.,500.,450.,400.,350.,300.,250., 
                 225.,200.,175.,150.,137.5,125.,112.5,100., 
                 87.5,75.,62.5]

## Add microphysics constants
cpmd = 0.887 # for Cp_moist
rgasmd = 0.608 # for R_moist
skewness = 75

def draw_isotherms():
    itherms = [[]]
    itherms.append([])
    for n in [-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40]:
        itherms[0].append(n)
        itherms[1].append(n + skewness)
    endvals = [1000.,100.]
    for n in range(len(itherms[1])):
        endpoints = []
        endpoints.append(itherms[0][n])
        endpoints.append(itherms[1][n])
        if itherms[0][n] < 0:
            pylab.plot(endpoints,endvals,color = 'blue', \
                linestyle='dashed', linewidth = .5)
        elif itherms[0][n] > 0:
            pylab.plot(endpoints,endvals,color = 'red', \
                linestyle='dashed', linewidth = .5)
        else:
            pylab.plot(endpoints,endvals,color = 'blue', \
                linestyle='solid', linewidth = .5)
            
def draw_dry_adiabats():
    dadiabats = []
    for n in [-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150]:
        theta = n + 273
        curadiabat = []
        for m in range(len(press_levels)):        
            curadiabat.append((theta * np.power(press_levels[m]/1000.,(287.04/1004.)) - 273) + (skewness * math.log(1000./press_levels[m],10)))
        dadiabats.append(curadiabat)
    endvals = press_levels
    for n in range(len(dadiabats)):
        pylab.plot(dadiabats[n],endvals,color = 'brown', \
            linestyle='dashed', linewidth = .5)
        
def draw_isobars():
    endxs = [-40,50]
    for n in range(100,1000,100):
        endys = [n,n]
        pylab.plot(endxs, endys, color = 'k',\
            linewidth = .5) 
        
def calc_es(T_K):
    ## Compute saturation vapor pressure from temperature in Kelvin?
    es = 6.112*np.exp(17.67*(T_K-273)/((T_K-273)+243.5))
    return es 

def draw_moist_adiabats(skewness,parcel,parcel_T,parcel_P):
    def calc_gamma_s(T,p):
        p = p/1000.
        es = .611 * math.exp(5423*((1/273.)-(1./(T))))
        r = epsilon * es / (p - es)
        gamma_d = .0098
        a = 0.28571
        b = 1.35e7
        c = 2488.4
        numer = (a * T) + (c * r)
        denom = p * (1 + (b * r / (T*T))) 
        gamma_s = numer/denom
        # print "P: %f T: %f es: %f Gs: %f" % (p, T, es, gamma_s)
        return gamma_s    

    def RK_4(p,T,dp):
        k1 = calc_gamma_s(T,p)
        k2 = calc_gamma_s((T + .5 * dp * k1),(p + .5 * dp))
        k3 = calc_gamma_s((T + .5 * dp * k2),(p + .5 * dp))
        k4 = calc_gamma_s((T + dp * k3),(p + dp))
    
        Tp = T + dp * (1/6.) * (k1 + 2. * k2 + 2. * k3 + k4)
        return Tp

    def RK_45(p,T,dp):
        k1 = calc_gamma_s(T,p)
        k2 = calc_gamma_s((T + .25 * dp * k1),(p + .25 * dp))
        k3 = calc_gamma_s((T + 3/32. * dp * k1 + 9/32. * dp * k2),(p + 3/8. * dp))
        k4 = calc_gamma_s((T + dp * 1932/2197. * k1 - 7200/2197. * dp * k2 + 7296/2197. * dp * k3),(p + 12/13. * dp))
        k5 = calc_gamma_s((T + dp * 439/216. * k1 - 8. * dp * k2 + 3680/513. * dp * k3 - 845/4104 * dp * k4),(p + dp))
        k6 = calc_gamma_s((T - dp * 8/27. * k1 + 2. * dp * k2 - 3544/2565. * dp * k3 + 1859/4104 * dp * k4 - 11/40. * dp * k5),(p + .5 * dp))
    
        Tp = T + dp * (16/135. * k1 + 6656/12825. * k3 + 28561/56430. * k4 - 9/50. * k5 + 2/55. * k6)
        return Tp

    if parcel == 0:
        for base_T in [-10,0,10,14,18,22,26,31,35]:
            madiabat = []
            madiabat_trans = []
            Tp = base_T + 273.
            
            for z in range(len(press_levels)):    
                p = press_levels[z]
                if z == 0:
                    dp = 0
                else:
                    dp = p - press_levels[z-1]        
                dt = calc_gamma_s(Tp,p) * (dp/1000.)

                new_Tp = RK_4(p*100.,Tp,(dp/10.)) 
                Tp = new_Tp
                madiabat.append(Tp)
                Tp_C = Tp - 273.
                Tp_trans = Tp_C + skewness * math.log(1000./p,10)

                madiabat_trans.append(Tp_trans)    
            pylab.semilogy(madiabat_trans, press_levels, color = 'green', basey = 10, linestyle = 'dotted', linewidth = .5)

    if parcel == 1:
        madiabat = []
        madiabat_trans = []
        Tp = parcel_T + 273.
        madiabat_press_levels = []
        for z in range(len(press_levels[parcel_P:])):    
            p = press_levels[z+parcel_P]
            if z == 0:
                dp = 0
            else:
                dp = p - press_levels[z-1 + parcel_P]        
            dt = calc_gamma_s(Tp,p) * (dp/1000.)

            new_Tp = RK_45(p*100.,Tp,(dp/10.)) 
            Tp = new_Tp
            madiabat.append(Tp)
            Tp_C = Tp - 273.
            Tp_trans = Tp_C + skewness * math.log(1000./p,10)

            madiabat_trans.append(Tp_trans)    
            madiabat_press_levels.append(p)
        return np.subtract(madiabat, 273)

def draw_parcel_trace(Tparc, Press):

    # Convert Pressures to log scale
    Pfact = np.multiply(skewness,np.log10(np.divide(1000., Press)))

    minlen = min(len(Tparc), len(Pfact))
    
    dry_parcel_trace_trans = np.add(Tparc[:minlen],Pfact[2:minlen+2])    

    pylab.semilogy(dry_parcel_trace_trans,Press[2:minlen+2],\
        basey=10, color = 'brown', linestyle = 'dashed',\
        linewidth = 1.5)

def new_draw_parcel_trace(Tb, PLCL, Press):

    # Convert Pressures to log scale
    Pfact = np.multiply(skewness,np.log10(np.divide(1000., Press)))

    parcelT = []
    flag = 1

    for p in range(len(Press)):
        if Press[p] >= PLCL:
            newTB = ((Tb + 273.) * (Press[p]/Press[0]) ** (287.04/1004.)) - 273.
            parcelT.append(newTB)
        else:
            if flag:
                if p == 0:
                    moists = draw_moist_adiabats(0, 1, Tb, 0)
                else:
                    moists = draw_moist_adiabats(0,1,parcelT[p-1], (p - 1 + len(press_levels) - len(Press)))
                for m in moists:
                    parcelT.append(m)
                flag = 0


    minlen = min(len(parcelT), len(Pfact))
    
    dry_parcel_trace = np.add(parcelT[:minlen], Pfact[:minlen])

    pylab.semilogy(dry_parcel_trace,Press[:minlen],\
        basey=10, color = 'brown', linestyle = 'dotted',\
        linewidth = 1.5)
    
    