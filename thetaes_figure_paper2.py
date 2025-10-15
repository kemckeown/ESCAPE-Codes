#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 17:52:11 2025

@author: kem6245
"""

## Panel (e) for synoptic figure looking at theta-es for cmas and skyler


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from datetime import timedelta
import matplotlib.dates as mdates
from matplotlib.dates import date2num, DateFormatter
import datetime as dt
from metpy.interpolate import interpolate_1d
from metpy.units import units
from metpy.calc import dewpoint_from_relative_humidity
from metpy.calc import mixing_ratio_from_relative_humidity
from metpy.calc import virtual_potential_temperature
from metpy.calc import equivalent_potential_temperature
from metpy.calc import dewpoint_from_relative_humidity
from metpy.calc import saturation_equivalent_potential_temperature
from metpy.calc import potential_temperature
import os
import metpy.calc as mpcalc
from scipy.signal import savgol_filter

def str2dt(time):
    return datetime.strptime(time, '%H:%M:%S')


ws2 = pd.read_csv('/Users/kem6245/Documents/Python Copy/ESCAPE/2 June 2022 Windsonde Data/SKYLER_Individual/2022-06-02_1417.raw_history.csv')        
ws1 = pd.read_csv('/Users/kem6245/Documents/Python Copy/ESCAPE/2 June 2022 Windsonde Data/Updated_Windsonds/CMAS/a_cmas_1428.csv')

h1 = ws1[' Altitude (m AGL)']
p1 = ws1[' Pressure (Pascal)']
p1 = pd.to_numeric(p1, errors='coerce')
t1 = ws1[' Temperature (C)']
t1 = pd.to_numeric(t1, errors='coerce')
rh1 = ws1[' Relative humidity (%)']
rh1 = pd.to_numeric(rh1, errors='coerce')
rs1 = ws1[' Rise speed (m/s)']
wdir1 = ws1[' Heading (degrees)']

## Time and Location 
time_raw1 = ws1['UTC time'] 
time1 = time_raw1.apply(str2dt)
## Now we need to get a few of our variables into unit form for the metpy calculations 
p1_hpa = (p1/100)
t_1 = t1.values * units.degC
p1_hpa1 = p1_hpa.values * units.hPa
h_1 = h1.values * units.m
rh_1 = rh1.values * units.percent

## We also must calculate theta-v 
mr1 = mixing_ratio_from_relative_humidity(p1_hpa1, t_1, rh_1)
td1 = dewpoint_from_relative_humidity(t_1, rh_1)
theta1 = potential_temperature(p1_hpa1, t_1)
thetav1 = virtual_potential_temperature(p1_hpa1, t_1, mr1)
thetae1 = equivalent_potential_temperature(p1_hpa1, t_1, td1)
thetaes1 = saturation_equivalent_potential_temperature(p1_hpa1, t_1)


h2 = ws2[' Altitude (m AGL)']
p2 = ws2[' Pressure (Pascal)']
t2 = ws2[' Temperature (C)']
rh2 = ws2[' Relative humidity (%)']
rs2 = ws2[' Rise speed (m/s)']
wdir2 = ws2[' Heading (degrees)']
#wd2_shift = wind_shift(wdir2)
wsp2 = ws2[' Speed (m/s)']
## Time and Location 
time_raw2 = ws2['UTC time'] 
time2 = time_raw2.apply(str2dt)
## Now we need to get a few of our variables into unit form for the metpy calculations 
p2_hpa = (p2/100)
t_2 = t2.values * units.degC
p2_hpa2 = p2_hpa.values * units.hPa
h_2 = h2.values * units.m
rh_2 = rh2.values * units.percent
## We also must calculate theta-v 
mr2 = mixing_ratio_from_relative_humidity(p2_hpa2, t_2, rh_2)
td2 = dewpoint_from_relative_humidity(t_2, rh_2)
theta2 = potential_temperature(p2_hpa2, t_2)
thetav2 = virtual_potential_temperature(p2_hpa2, t_2, mr2)
thetae2 = equivalent_potential_temperature(p2_hpa2, t_2, td2)
thetaes2 = saturation_equivalent_potential_temperature(p2_hpa2, t_2)

thetaes1 = savgol_filter(thetaes1.magnitude, 100, 3)
thetaes2 = savgol_filter(thetaes2.magnitude, 100, 3)



rs1s = pd.to_numeric(rs1, errors='coerce')
rs2s = pd.to_numeric(rs2, errors='coerce')

plt.figure(figsize=(7,7))
plt.plot(thetaes1[rs1s >= 1], h1[rs1s >= 1], color='mediumpurple', alpha=0.7, linewidth=4, label='CMAS')
plt.plot(thetaes2[rs2s >= 1], h2[rs2s >= 1], color='green', linewidth=4, label='SKYLER')

# Add labels, legend, and grid
plt.xlabel(r'$\theta_{es}$ (K)', fontsize=14)
plt.ylabel('Height (m)', fontsize=14)
plt.ylim(0,4500)
plt.xlim(335,370)
plt.title('(e) Saturation Equivalent Potential Temperature', fontsize=14)
plt.legend(fontsize=10)
plt.grid(True, linestyle='--', alpha=0.7)

# Show the plot
plt.tight_layout()

plt.savefig('/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 2/thetaes.png', dpi=300)



