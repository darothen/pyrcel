# -*- coding:utf-8 -*-
"""
@author Geir Arne Waagbø
@see http://code.google.com/p/pywrfplot/
"""
import matplotlib.pyplot as plt
import numpy as np

# Directory where WRF nc-files are located
directory = 'C:/Dokumenter/MET/master/gmoen/big1km/'

# The date of the WRF simulation file
date = '2010-01-15'

# The number of time steps in the WRF simulation
Nt = 73

# use 'l' for faster response, use 'h' for high resolution
mapResolution = 'h' 

# This defines the grid box for which skewT-plots etc are based
lat_focuspoint = 60.2
lon_focuspoint = 11.08

# Additional grid point that will be marked by a red dot on some of the maps
# (set to -1 to avoid this red dot!)
lat_rg = 59.37
lon_rg = 10.78

# P_top must be the same as what is used in the WRF simulation
P_top = 10**4

P_bot = 10**5

# Max height used on z-axis for xz-plots
z_max = 4000.0

# Tick increment used for xz-plot
dz = 200

# Max height used when plotting terrain contours
max_h = 800

# Pressure interval used when plotting contour lines on top level map
pmsl_int = np.arange(960.,1040.,4)

# Temperature interval used when plotting contour lines
temp_int = np.arange(-80.0,50.0,2.0)

# levels used for xz-cloud plots
xz_cloudwater_levels = np.arange(0.08, 0.7, 0.08)
xz_rain_levels = np.arange(0.003, 0.0110, 0.0015)
xz_snow_levels = np.arange(0.06, 0.17, 0.02)

# levels used for tz-cloud plots
tz_cloudwater_levels = np.arange(0.08, 0.7,0.08)
tz_rain_levels = np.arange(0.003, 0.0110, 0.0015)
tz_snow_levels = np.arange(0.02, 0.10, 0.02)

WaterColor = "#B2FFFF"
LandColor = "#FFCCB9"
barb_increments = {'half': 2.5,'full':5.0,'flag':25.0}

cmap_red = plt.get_cmap('Reds')
cmap_green = plt.get_cmap('Greens')
cmap_blue = plt.get_cmap('Blues')
cmap_grey = plt.get_cmap('Greys')
cmap_jet = plt.get_cmap('Jet')



T_base = 300.0
T_zero = 273.15
L = 2.501e6 # latent heat of vaporization
R = 287.04  # gas constant air
Rv = 461.5  # gas constant vapor
eps = R/Rv
cp = 1005.
cv = 718.
kappa = (cp-cv)/cp
g = 9.81



