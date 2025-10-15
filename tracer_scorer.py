#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 11:42:39 2024

@author: kem6245
"""

## Scorer parameter for TRACER radiosondes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits import mplot3d

import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator
from metpy.units import units
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection
from metpy.calc import wind_speed
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D
from netCDF4 import Dataset
from metpy.calc import relative_humidity_from_dewpoint
import math
from scipy import signal
import xarray as xr

WSR1 = xr.open_dataset('file')

print(WSR1.data_vars)
sfc_ht_lp = 6
# P_hPa1 = WSR1['pres'].values
# P_Pa1 = P_hPa1 * 100
# T_C1 = WSR1['tdry'].values
# T_C1_units = T_C1 * units.degC 
# td = WSR1['dp'].values * units.degC
# RH1 = relative_humidity_from_dewpoint(T_C1_units, td)
# ALT1 = WSR1['alt'].values
# HGHT1 = WSR1['alt'].values 
# U1 = WSR1['u_wind'].values
# V1 = WSR1['v_wind'].values
# U1_units = WSR1['u_wind'].values * units('m/s')
# V1_units = WSR1['v_wind'].values * units('m/s')
# lon_in = WSR1['lon'].values
# lat_in = WSR1['lat'].values
# WNDSPD1 = wind_speed(U1_units, V1_units)

P_hPa1 = WSR1['pres'].values
P_Pa1 = P_hPa1 * 100
T_C1 = WSR1['tdry'].values
RH1 = WSR1['rh'].values
ALT1 = WSR1['alt'].values
HGHT1 = ALT1 - sfc_ht_lp 
##HGHT1 = WSR1[]
U1 = WSR1['u_wind'].values
V1 = WSR1['v_wind'].values
lon_in = WSR1['lon'].values
lat_in = WSR1['lat'].values
WNDSPD1 = WSR1['wspd'].values



#CONSTANTS
Rd = 287.04 #J/kg/K - Dry air gas constant
Rv = 461.5  #J/kg/K - vapor gas constant
cp = 1005   #specific heat capacity J/kg
cpv = 1875  #specific heat capacity of vapor ar at constant pressure J/kg/K
cvv = 1410  #specific heat capacity of vapor at constant volume J/kg/K
cliq = 4200 #specific heat capacity of liquid J/kg/K
eps = 0.622 #constant in eqns
es0 = 611.2 #e_s at T = 0 degC
T0 = 273.15 #K
p0 = 100000 #sfc P in Pa
lv = 2.5e6  #J/kg enthalpy of vaporization
k = Rd/cp

#CALCULATE TEMPERAURE IN K
T_K1 = np.zeros(len(P_Pa1))
T_K1[:] = T_C1[:]+T0

#CALCUALTE SATURATION VAPOR PRESSURE (Pa)
esT_1 = np.zeros(len(P_Pa1))
esT_1[:] = es0 * np.exp((17.67*T_C1[:])/(T_C1[:] + 243.5))

#CALCULATE ACTUAL VAPOR PRESSURE USING Td (Pa)
e_1 = np.zeros(len(P_Pa1))
e_1[:] = esT_1[:] * (RH1[:] / 100)

#CALCULATE DEW POINT (C) (MARKOWSKI AND RICHARDSON 2011) - NEED TO CONVERT VAPOR PRESSURE FROM PA TO HPA   
Td_C1 = np.zeros(len(P_Pa1))
Td_C1[:] = 243.5 / ((17.67 / np.log((e_1[:]*0.01)/6.112)) - 1)


#CALCULATE VAPOR MIXING RATIO (kg/kg)
qv_1 = np.zeros(len(P_Pa1))
qv_1[:] = eps * e_1[:] / (P_Pa1[:] - e_1[:])
rv_1 = np.zeros(len(P_Pa1))
rv_1[:] = 1000 * qv_1[:]  #BOLTON EXPRESSES mixing ratio in g/kg (page 1047)


#CALCULATE POTENTIAL TEMPERATURE (BOLTON EQUATION 7)
theta_1 = np.zeros(len(P_Pa1))
theta_1[:] = T_K1[:] * ((p0/P_Pa1[:])**(k*(1-0.28e-3 * rv_1[:])))


# SMOOTH USING A SAVITZKY-GOLAY FILTER
smoothwindow = 100
Tsm_C1 = signal.savgol_filter(T_C1, window_length=smoothwindow, polyorder=3, mode="nearest")
Tdsm_C1 = signal.savgol_filter(Td_C1, window_length=smoothwindow, polyorder=3, mode="nearest")
U1_smooth =  signal.savgol_filter(U1, window_length=smoothwindow, polyorder=3, mode="nearest")
V1_smooth = signal.savgol_filter(V1, window_length=smoothwindow, polyorder=3, mode="nearest")
WNDSPD1_smooth = signal.savgol_filter(WNDSPD1, window_length=smoothwindow, polyorder=3, mode="nearest")
theta1_smooth = signal.savgol_filter(theta_1, window_length=smoothwindow, polyorder=3, mode="nearest")
ALT1_smooth = signal.savgol_filter(ALT1, window_length=smoothwindow, polyorder=3, mode="nearest")

k = 0
boreangle = 116.25 #mosaic radar 116.25; san antonio radar ranges 112.5 - 95 degrees; the angle off of mathematical zero
levels = len(P_Pa1)
degrees = np.zeros(len(P_Pa1))
diffangle = np.zeros(len(P_Pa1))
parallelwind = np.zeros(len(P_Pa1))
k2 = np.zeros(len(P_Pa1))
wavelength = 13000 #san antonio radar estimate 13 km or 13000 m is the horizontal wavelength of the waves
while k < levels:
    degrees[k] = (np.arctan2((-1*V1_smooth[k]),(-1*U1_smooth[k])))*(180/np.pi) #radians answer converted into mathematical degrees
    if degrees[k] < 0:
        degrees[k] = 360 + degrees[k] #converted to angles off of mathematical zero
    diffangle[k] = boreangle - degrees[k]
    parallelwind[k] = WNDSPD1_smooth[k] * np.cos(diffangle[k] * (np.pi / 180))
    k2[k] = ((2 * np.pi) / wavelength)**2 #Equals 3.95 x 10-7 m-2
    k = k + 1


parallelwind_smooth = signal.savgol_filter(parallelwind, window_length=smoothwindow, polyorder=3, mode="nearest") 
"""
#degreessm=savitzky_golay(degrees,101,3) #totalwind
avgparallelwindlower = np.mean(parallelwind[20:94])
avgparallelwindupper = np.mean(parallelwind[135:218])
"""

dthetadz =  np.zeros(levels-1)
dthetadz_smooth =  np.zeros(levels-1)
N2 =  np.zeros(levels-1)
N =  np.zeros(levels-1)
N2_no_nan =  np.zeros(levels-1)
N_no_nan =  np.zeros(levels-1)
N2_smooth =  np.zeros(levels-1)
N_smooth =  np.zeros(levels-1)
N2_no_nan_smooth =  np.zeros(levels-1)
N_no_nan_smooth =  np.zeros(levels-1)
avgtheta = np.zeros(levels-1)
avgtheta_smooth = np.zeros(levels-1)
P_Pa_N = np.zeros(levels-1)
alt_m_N = np.zeros(levels-1)
uavg = np.zeros(levels-1)
u2 = np.zeros(levels-1)
udiff = np.zeros(levels-1)
zdiff = np.zeros(levels-1)
dudz = np.zeros(levels-1)
zeroline = np.zeros(levels-1)
term1 = np.zeros(levels-1)
Uminusc = np.zeros(levels-1)

g = 9.8
z = 0
#c = 17.5 #NEED TO CALCULTE THIS FOR ESCAPE WAVE - theoretical wave speed
#c = 10 #san antonio radar
c = 14 #mosaic radar
while z < levels-1:
    #Calculate N
    dthetadz[z] = (theta_1[z+1] - theta_1[z]) / (ALT1_smooth[z+1] - ALT1_smooth[z])
    avgtheta[z] = (theta_1[z+1]+theta_1[z])/2
    N2[z] = ((g/avgtheta[z])*dthetadz[z])
    N[z] = (N2[z])**0.5
    dthetadz_smooth[z] = (theta1_smooth[z+1] - theta1_smooth[z]) / (ALT1_smooth[z+1] - ALT1_smooth[z])
    avgtheta_smooth[z] = (theta1_smooth[z+1]+theta1_smooth[z])/2
    N2_smooth[z] = ((g/avgtheta_smooth[z])*dthetadz_smooth[z])
    N_smooth[z] = (N2[z])**0.5
    P_Pa_N[z] = P_Pa1[z]
    alt_m_N[z] = ALT1_smooth[z]
    #Calculate avgerage parallel wind in the layer
    uavg[z] = (parallelwind_smooth[z]+parallelwind_smooth[z+1]) / 2
    #Calculate parallel wind squared
    u2[z] = (uavg[z]-c)**2
    #Calculate U - c
    Uminusc[z] = (parallelwind_smooth[z]-c)
    #Term 1
    #term1[z] = N2[z] / u2[z]
    term1[z] = N2_smooth[z] / u2[z]
    #Calculate dU/dz
    udiff[z] = parallelwind_smooth[z+1] - parallelwind_smooth[z]
    zdiff[z] = ALT1_smooth[z+1] - ALT1_smooth[z]
    dudz[z] = udiff[z] / zdiff[z]
    if term1[z] < 0:
        print(ALT1_smooth[z])
    z = z+1

a=pd.Series(N)
N_no_nan = a.interpolate(method='linear')
a=pd.Series(N2)
N2_no_nan = a.interpolate(method='linear')
a=pd.Series(N_smooth)
N_no_nan_smooth = a.interpolate(method='linear')
a=pd.Series(N2_smooth)
N2_no_nan_smooth = a.interpolate(method='linear')


dudz2 = np.zeros(levels-2)
l2 = np.zeros(levels-2)
l2_no_nan = np.zeros(levels-2)
alt_m_Sc = np.zeros(levels-2)
term2 = np.zeros(levels-2)
z = 0
while z < levels-2:
    #Calculate d2u/dz2
    dudz2[z] = (dudz[z+1]-dudz[z]) / (ALT1[z+1] - ALT1[z])
    #term 2
    term2[z] = (dudz2[z] / (uavg[z]-c)) #-1 * (dudz2[z] / (uavg[z]-c))
    #Calculate l2
    l2[z] = term1[z] - (dudz2[z] / (uavg[z]-c))
    alt_m_Sc[z] = ALT1[z]
    z = z+1
#a=pd.Series(l2)
#l2_no_nan = a.interpolate(method='linear')

###PLOT
#fig, ax1 = plt.subplots(figsize=(9.5,5))

# color palette generator https://waldyrious.net/viridis-palette-generator/
"""
VAR = 'SMOOTH_T_TD'
ax1.plot(T_C1, ALT1, color = 'black', lw = 2, alpha=1, label = hour1)
ax1.plot(Td_C1, ALT1, color ='blue', lw = 2, alpha=1, label = hour1)
ax1.plot(Tsm_C1, ALT1, color = 'red', lw = 2, alpha=1, label = hour1)
ax1.plot(Tdsm_C1, ALT1, color ='orange', lw = 2, alpha=1, label = hour1)
ax1.plot(Td1_smooth, ALT1_smooth, color = 'green', lw = 2, alpha=1, label = hour1)
ax1.plot(T1_smooth, ALT1_smooth, color ='green', lw = 2, alpha=1, label = hour1)
"""
"""
VAR = 'SMOOTH_PARALLELWIND'
ax1.plot(parallelwind, ALT1, color = 'black', lw = 2, alpha=1, label = hour1)
ax1.plot(WNDSPD1_smooth, ALT1, color ='blue', lw = 2, alpha=1, label = hour1)
ax1.plot(WNDSPD1, ALT1, color = 'red', lw = 2, alpha=1, label = hour1)
"""

"""
####OTHER
ax1.set_xlabel('m s-1', fontsize=12)
ax1.set_ylabel('height (m)', fontsize=12)
plt.title('2 June 2022 '+TRUCK+' Windsonds '+VAR+' '+UPDOWN, fontsize = 16)
ax1.set_xlim([-5,5])
#ax1.set_xticks([0,30,60,90,120,150,180,210,240,270,300,330,360])
#x_ticks_labels = (['180','210','240','270','300','330','0','30','60','90','120','150','180'])
ax1.set_ylim([0,4500])
ax1.tick_params(axis='y', which='major', labelsize=12, direction='out')
#ax1.set_xticklabels(x_ticks_labels, rotation=0, fontsize=12)
#ax1.xaxis.set_tick_params(which='major', labelsize=12, direction='out') #x-axis characteristics
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
"""


###FOR PLOTTING N2 OVERLAID WITH UWIND
# VAR = 'N2_UWIND'
# fig, ax1 = plt.subplots(figsize=(8, 6))
# ax1.set_xlabel('N$^{2}$ (s$^{-2})$', fontsize=14)
# ax1.set_ylabel('Height (m)', fontsize=14)
# xmin=-0.001
# xmax=0.001
# ymin=ALT1[0]
# ymax=5000
# ax1.axis([xmin,xmax,ymin,ymax])
# #ax1.set_xbound(-0.001, 0.001)
# #ax1.set_ybound(alt_m[0], 6000)
# #ax1.plot(N2_no_nan, alt_m_N, color ='darkred', lw = 2, alpha=0.5)
# ax1.plot(N2_no_nan_smooth, alt_m_N, color ='darkred', lw = 2, alpha=1)
# ax1.plot(zeroline, alt_m_N, color ='gray', lw = 2, linestyle = 'dashed', alpha=1)
# ax1.tick_params(axis='x', which='major', labelsize=12, direction='out', colors='darkred')
# ax1.tick_params(axis='y', which='major', labelsize=12, direction='out')
# ax1.ticklabel_format(style='sci', axis = 'x', scilimits=(-4,-4))
# #ax1.spines['bottom'].set_color('darkred') 

# ax2 = ax1.twiny()  # instantiate a second axes that shares the same x-axis
# ax2.plot(parallelwind, ALT1, color ='darkorange', lw = 2, alpha=1)
# #ax2.set_xlabel('m/2', fontsize=12)
# #ax2.set_xticks(-2.5, 15)
# #ax2.set_xbound(-2.5, 15)
# #ax2.axis([-2.5,15,ymin,ymax])
# ax2.axis([-10,10,ymin,ymax])
# #ax2.set_xticklabels([-2.5, 0, 2.5, 5, 7.5, 10, 12.5, 15])
# ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
# ax2.tick_params(axis='x', which='major', labelsize=12, direction='out', colors='darkorange')
# #title = ax1.set_title("v wind (m s$^{-1})$")
# #title = ax1.set_title("v wind (m s$^{-1})$")
# title = ax1.set_title("(a) Wind Along the Boundary-Motion Axis (m s$^{-1})$", fontsize=18)
# title.set_y(1.1)
# fig.subplots_adjust(top=0.85)
# ax1.grid(True)



# #FOR PLOTTING TERM1 AND TERM 2
# VAR = 'TERM1TERM2'
# fig, ax1 = plt.subplots(figsize=(8, 6))
# #ax1.set_xlabel('l$^{2}$ (m$^{-2}$)', fontsize=12)
# ax1.set_xlabel('d$^{2}$u/dz$^{2}$ / (u-c) [m$^{-2}$]', fontsize=14)
# ax1.set_ylabel('Height (m)', fontsize=14)
# xmin=-0.00001
# xmax=0.00001
# ymin=ALT1[0]
# ymax=5000
# ax1.axis([xmin,xmax,ymin,ymax])
# #ax1.plot(l2, alt_m_Sc, color ='cornflowerblue', lw = 2, alpha=1)
# ax1.plot(term2, alt_m_Sc, color ='cornflowerblue', lw = 2, alpha=1)
# ax1.ticklabel_format(style='sci', axis = 'x', scilimits=(-5,-5))
# ax1.tick_params(axis='x', which='major', labelsize=12, direction='out', colors='cornflowerblue')
# ax1.tick_params(axis='y', which='major', labelsize=12, direction='out')

# ax2 = ax1.twiny()  # instantiate a second axes that shares the same x-axis
# #ax2.plot(term1, alt_m_N, color ='darkblue', lw = 2, alpha=1)
# ax2.plot(term1, alt_m_N, color ='purple', lw = 2, alpha=0.7)
# ax2.axis([-0.00001,0.00001,ymin,ymax])
# ax2.ticklabel_format(style='sci', axis = 'x', scilimits=(-5,-5))
# ax2.tick_params(axis='x', which='major', labelsize=12, direction='out', colors='purple')
# title = ax1.set_title("(b) N$^{2}$ / (u-c)$^{2}$ [m$^{-2}$]", fontsize = 18)
# title.set_y(1.1)
# fig.subplots_adjust(top=0.85)
# ax1.plot(zeroline, alt_m_N, color ='black', lw = 2, linestyle = 'dashed', alpha=1)
# ax1.grid(True)


# VAR = 'SCORER'
# fig, ax1 = plt.subplots(figsize=(8, 6))
# ax1.set_xlabel('l$^{2}$ (m$^{-2}$)', fontsize=14)
# ax1.set_ylabel('Height (m)', fontsize=14)
# ax1.set_title('(c) Scorer Parameter', fontsize=18)
# xmin=-0.0000075
# xmax=0.0000075
# ymin=ALT1[0]
# ymax=5000
# ax1.axis([xmin,xmax,ymin,ymax])
# ax1.plot(l2, alt_m_Sc, color ='darkblue', lw = 2, alpha=1)
# ax1.yaxis.set_major_locator(MaxNLocator(5)) 
# ax1.ticklabel_format(style='sci', axis = 'x', scilimits=(-5,-5))
# ax1.tick_params(axis='x', which='major', labelsize=12, direction='out', colors='darkblue')
# ax1.tick_params(axis='y', which='major', labelsize=12, direction='out')
# ax1.plot(zeroline, alt_m_N, color ='magenta', lw = 2, linestyle = 'dashed', alpha=1)
# ax1.grid(True)



###FOR PLOTTING k^2 and l^2
# VAR = 'K2L2'
# fig, ax1 = plt.subplots()
# ax1.set_xlabel('l$^{2}$ (m$^{-2}$)', fontsize=14)
# ax1.set_ylabel('Height (m)', fontsize=14)
# ax1.set_title('Scorer', fontsize=18)
# xmin=-0.00003
# xmax=0.00003
# ymin=ALT1[0]
# ymax=5000
# ax1.axis([xmin,xmax,ymin,ymax])
# ax1.plot(l2, alt_m_Sc, color ='darkblue', lw = 2, alpha=1)
# ax1.ticklabel_format(style='sci', axis = 'x', scilimits=(-5,-5))
# ax1.tick_params(axis='x', which='major', labelsize=12, direction='out', colors='darkblue')
# ax1.tick_params(axis='y', which='major', labelsize=12, direction='out')

# ax2 = ax1.twiny()  # instantiate a second axes that shares the same x-axis
# ax2.plot(k2, ALT1, color ='purple', lw = 2, alpha=0.7)
# ax2.axis([-0.00003,0.00003,ymin,ymax])
# ax2.ticklabel_format(style='sci', axis = 'x', scilimits=(-5,-5))
# ax2.tick_params(axis='x', which='major', labelsize=12, direction='out', colors='purple')
# title = ax1.set_title("k$^{2}$ (m$^{-2}$)")
# title.set_y(1.1)
# fig.subplots_adjust(top=0.85)
# ax1.plot(zeroline, alt_m_N, color ='lightblue', lw = 2, linestyle = 'dashed', alpha=1)
# ax1.grid(True)



#FOR PLOTTING N OVERLAID WITH U-c
VAR = 'N_UminusC'
fig, ax1 = plt.subplots(figsize=(8, 6))
ax1.set_xlabel('N (s$^{-1})$', fontsize=14)
ax1.set_ylabel('Height (m)', fontsize=14)
xmin=-0.1
xmax=0.1
ymin=ALT1[0]
ymax=5000
ax1.axis([xmin,xmax,ymin,ymax])
ax1.plot(N_no_nan_smooth, alt_m_N, color ='darkred', lw = 2, alpha=1)
ax1.plot(zeroline, alt_m_N, color ='gray', lw = 2, linestyle = 'dashed', alpha=1)
ax1.tick_params(axis='x', which='major', labelsize=12, direction='out', colors='darkred')
ax1.tick_params(axis='y', which='major', labelsize=12, direction='out')
#ax1.ticklabel_format(style='sci', axis = 'x', scilimits=(-2,-2))


ax2 = ax1.twiny()  # instantiate a second axes that shares the same x-axis
ax2.plot(Uminusc, alt_m_N, color ='darkorange', lw = 2, alpha=1)
ax2.axis([-25,25,ymin,ymax])
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.tick_params(axis='x', which='major', labelsize=12, direction='out', colors='darkorange')
title = ax1.set_title("(d) U - c (m s$^{-1})$", fontsize=18)
title.set_y(1.1)
fig.subplots_adjust(top=0.85)
ax1.grid(True)


plt.savefig('file_figure', dpi = 300)



