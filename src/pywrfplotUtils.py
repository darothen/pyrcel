# -*- coding:utf-8 -*-
"""
@author Geir Arne Waagbø
@see http://code.google.com/p/pywrfplot/
"""
import math
import numpy as np
from netCDF4 import Dataset

from pywrfplotParams import *

# constants used to calculate moist adiabatic lapse rate
# See formula 3.16 in Rogers&Yau
a = 2./7.
b = eps*L*L/(R*cp)
c = a*L/R

def gamma_s(T,p):
    """Calculates moist adiabatic lapse rate for T (Celsius) and p (Pa)
    Note: We calculate dT/dp, not dT/dz
    See formula 3.16 in Rogers&Yau for dT/dz, but this must be combined with
    the dry adiabatic lapse rate (gamma = g/cp) and the 
    inverse of the hydrostatic equation (dz/dp = -RT/pg)"""
    esat = es(T)
    wsat = eps*esat/(p-esat) # Rogers&Yau 2.18
    numer = a*(T+T_zero) + c*wsat
    denom = p * (1 + b*wsat/((T+T_zero)**2)) 
    return numer/denom # Rogers&Yau 3.16

def es(T):
    """Returns saturation vapor pressure (Pascal) at temperature T (Celsius)
    Formula 2.17 in Rogers&Yau"""
    return 611.2*np.exp(17.67*T/(T+243.5))

def e(w,p):
    """Returns vapor pressure (Pa) at mixing ratio w (kg/kg) and pressure p (Pa)
    Formula 2.18 in Rogers&Yau"""
    return w*p/(w+eps)

def td(e):
    """Returns dew point temperature (C) at vapor pressure e (Pa)
    Insert Td in 2.17 in Rogers&Yau and solve for Td"""
    return 243.5 * np.log(e/611.2)/(17.67-np.log(e/611.2))

def interp(geopot, pres, p):
    """ Returns the interpolated geopotential at p using the values in pres. 
    The geopotential for an element in pres must be given by the corresponding
    element in geopot. The length of the geopot and pres arrays must be the same. 
    """
    if (len(geopot) != len(pres)):
        raise Exception, "Arrays geopot and pres must have same length"
    
    k = len(pres)-1
    while (k > 1 and pres[k-1] <= p):
        k = k-1

    if (pres[k] > p):
        w = 0.0
    elif (pres[k-1] < p):
        w = 1.0
    else: 
        w = (p-pres[k])/(pres[k-1]-pres[k])

    return (1.0-w)*geopot[k] + w*geopot[k-1]

def openWRF(nest):
    file_wrf = directory + 'wrfout_d0' + str(nest) + '_' + date + '.nc'
    try:
        return Dataset(file_wrf, 'r+')
    except Exception:
        print 'Found no WRF-file for nest: ' + str(nest) + '. Looked for: '+ file_wrf
        return None
    
def openWPS(nest):
    file_met = directory + 'met_em.d0' + str(nest) + '.' + date + '.nc'
    try:
        return Dataset(file_met, 'r+')
    except Exception:
        print 'Found no WPS-file for nest: ' + str(nest) + '. Looked for: '+ file_met
        return None

def getDimensions(nc):
    Nx = nc.getncattr('WEST-EAST_GRID_DIMENSION')-1
    Ny = nc.getncattr('SOUTH-NORTH_GRID_DIMENSION')-1
    Nz = nc.getncattr('BOTTOM-TOP_GRID_DIMENSION')-1    
    dx = nc.getncattr('DX')
    dy = nc.getncattr('DY')
    lons = nc.variables['XLONG'][0]
    lats = nc.variables['XLAT'][0]
    # find coordinates for focus point
    x,y = getXY(lons[Ny/2,:],lats[:,Nx/2])
    return Nx,Ny,Nz,lons,lats,dx,dy,x,y

def getXY(lon,lat):
    x_nr = 0
    y_nr = 0
    while (lon[x_nr] < lon_focuspoint):
        x_nr += 1
    while (lat[y_nr] < lat_focuspoint):
        y_nr += 1    
    
    print "x_nr: " + str(x_nr), " Lon: ", str(lon[x_nr])
    print "y_nr: " + str(y_nr), " Lat: ", str(lat[y_nr])
    return x_nr,y_nr


