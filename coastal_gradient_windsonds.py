#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 12:21:23 2024

@author: kem6245
"""

## Interpolation test

## Let's just interpolate one file

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

## For up, all files usable

## For down, get rid of the 'j' files

## Directories for CMAS and SKYLER Windsonds
cmas_dir = '/Users/kem6245/Documents/Python Copy/ESCAPE/2 June 2022 Windsonde Data/coastal_gradient/cmas/'
skyler_dir = '/Users/kem6245/Documents/Python Copy/ESCAPE/2 June 2022 Windsonde Data/coastal_gradient/skyler/'

# Common Height Grid
common_heights = np.arange(0, 4500, 25) * units.m

# # Function to read and process data from a single CSV file
# def process_file(filepath):
#     df = pd.read_csv(filepath)

#     # Ensure numeric conversions
#     df[' Altitude (m AGL)'] = pd.to_numeric(df[' Altitude (m AGL)'], errors='coerce')
#     df[' Pressure (Pascal)'] = pd.to_numeric(df[' Pressure (Pascal)'], errors='coerce') / 100  # Convert to hPa
#     df[' Temperature (C)'] = pd.to_numeric(df[' Temperature (C)'], errors='coerce')
#     df[' Relative humidity (%)'] = pd.to_numeric(df[' Relative humidity (%)'], errors='coerce')
    
#     # Assign variables with units
#     h = df[' Altitude (m AGL)'].values * units.m
#     t = df[' Temperature (C)'].values * units.degC
#     rh = df[' Relative humidity (%)'].values * units.percent

#     # Calculate dew point
#     td = dewpoint_from_relative_humidity(t, rh)

#     return h, t, td

# # List all CSV files in each directory
# cmas = sorted([os.path.join(cmas_dir, f) for f in os.listdir(cmas_dir) if f.endswith('.csv')])
# skyler = sorted([os.path.join(skyler_dir, f) for f in os.listdir(skyler_dir) if f.endswith('.csv')])


# # Loop through each pair of files
# for cmas_file, skyler_file in zip(cmas, skyler):
#     print(f"Processing {cmas_file} and {skyler_file}")

#     # Process files
#     h1, t1, td1 = process_file(cmas_file)
#     h2, t2, td2 = process_file(skyler_file)

#     # Interpolate onto a common height grid
#     temp_interp1 = interpolate_1d(common_heights, h1, t1)
#     temp_interp2 = interpolate_1d(common_heights, h2, t2)
def read_data(directory):
    files = sorted([os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.csv')])
    data_list = []
    
    for file in files:
        df = pd.read_csv(file)
        print(file)
        print()
        
        # Ensure numeric conversions and units
        df[' Altitude (m AGL)'] = pd.to_numeric(df[' Altitude (m AGL)'], errors='coerce')
        df[' Pressure (Pascal)'] = pd.to_numeric(df[' Pressure (Pascal)'], errors='coerce') / 100  # Convert to hPa
        df[' Temperature (C)'] = pd.to_numeric(df[' Temperature (C)'], errors='coerce')
        df[' Relative humidity (%)'] = pd.to_numeric(df[' Relative humidity (%)'], errors='coerce')
        df[' Rise speed (m/s)'] = pd.to_numeric(df[' Rise speed (m/s)'], errors='coerce')

        # Filter based on rise speed > 0
        df_filtered = df[df[' Rise speed (m/s)'] > 0]

        # Append filtered data
        data_list.append(df_filtered)
    
    return data_list

# Read data from both directories
cmas_data = read_data(cmas_dir)
skyler_data = read_data(skyler_dir)

# Prepare to store differences
temp_diff = []
dp_diff = []
theta_diff = []
theta_v_diff = []

# Process each pair of files
for file1, file2 in zip(cmas_data, skyler_data):
    # Extract necessary variables
    height1 = file1[' Altitude (m AGL)'].values * units.m
    temp1 = file1[' Temperature (C)'].values * units.degC
    rh1 = file1[' Relative humidity (%)'].values * units.percent
    pressure1 = file1[' Pressure (Pascal)'].values * units.hPa

    height2 = file2[' Altitude (m AGL)'].values * units.m
    temp2 = file2[' Temperature (C)'].values * units.degC
    rh2 = file2[' Relative humidity (%)'].values * units.percent
    pressure2 = file2[' Pressure (Pascal)'].values * units.hPa

    # Interpolate to common height grid
    temp1_interp = interpolate_1d(common_heights, height1, temp1)
    temp2_interp = interpolate_1d(common_heights, height2, temp2)
    
    # Calculate differences
    temp_diff.append(temp1_interp - temp2_interp)

    # Dew points
    dp1 = dewpoint_from_relative_humidity(temp1, rh1)
    dp2 = dewpoint_from_relative_humidity(temp2, rh2)

    # Interpolate dew points
    dp1_interp = interpolate_1d(common_heights, height1, dp1)
    dp2_interp = interpolate_1d(common_heights, height2, dp2)
    
    dp_diff.append(dp1_interp - dp2_interp)

    # Calculate potential temperature and virtual potential temperature
    theta1 = potential_temperature(pressure1, temp1)
    theta2 = potential_temperature(pressure2, temp2)
    
    # Interpolate potential temperatures
    theta_diff.append(interpolate_1d(common_heights, height1, theta1) - interpolate_1d(common_heights, height2, theta2))

    theta_v1 = virtual_potential_temperature(pressure1, temp1, rh1)
    theta_v2 = virtual_potential_temperature(pressure2, temp2, rh2)

    # Interpolate virtual potential temperatures
    theta_v_diff.append(interpolate_1d(common_heights, height1, theta_v1) - interpolate_1d(common_heights, height2, theta_v2))


print(temp1)
print()
print(temp2)
print()
print(temp1_interp)
print()
print(temp2_interp)
print()
# # Create a 2x2 figure for plotting
# fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# # Flatten the axes array for easier iteration
# axs = axs.flatten()
# cmap = plt.get_cmap('viridis')
# colors = cmap(np.linspace(0, 1, len(temp_diff)))
# ## For up:
# labels = ['1430 UTC', '1530 UTC', '1600 UTC', '1800 UTC', '1830 UTC', '1900 UTC', '2000 UTC', '2030 UTC', '2100 UTC', '2200 UTC', '2230 UTC']
# ## For down:
# ##labels = ['1430 UTC', '1530 UTC', '1600 UTC', '1800 UTC', '1830 UTC', '1900 UTC', '2000 UTC', '2030 UTC', '2100 UTC', '2230 UTC']
    
    
# # Plot temperature differences
# for i, diff in enumerate(temp_diff):
#     axs[0].plot(diff.magnitude, common_heights.magnitude/1000, color=colors[i], label=labels[i])
# axs[0].set_title('(A)', fontsize = 14)
# axs[0].set_xlabel('Temperature Difference (°C)', fontsize = 12)
# axs[0].set_ylabel('Height (km)', fontsize = 12)
# axs[0].grid(True)


# # Plot dew point differences
# for i, diff in enumerate (dp_diff):
#     axs[1].plot(diff.magnitude, common_heights.magnitude/1000, color=colors[i])
# axs[1].set_title('(B)', fontsize = 14)
# axs[1].set_xlabel('Dew Point Difference (°C)', fontsize = 12)
# axs[1].set_ylabel('Height (km)', fontsize = 12)
# axs[1].grid(True)


# # Plot potential temperature differences
# for i, diff in enumerate (theta_diff):
#     axs[2].plot(diff.magnitude, common_heights.magnitude/1000, color=colors[i])
# axs[2].set_title('(C)', fontsize = 14)
# axs[2].set_xlabel(r'$\theta$ (K)', fontsize = 12)
# axs[2].set_ylabel('Height (km)', fontsize = 12)
# axs[2].grid(True)


# # Plot virtual potential temperature differences
# for i, diff in enumerate (theta_v_diff):
#     axs[3].plot(diff.magnitude, common_heights.magnitude/1000, color=colors[i])
# axs[3].set_title('(D)', fontsize = 14)
# axs[3].set_xlabel(r'$\theta_v$ (K)', fontsize = 12)
# axs[3].set_ylabel('Height (km)', fontsize = 12)
# axs[3].grid(True)



# fig.suptitle('Coastal - Inland Gradient Ascending Windsonds', fontsize=18)
# fig.legend(labels, loc='lower center', bbox_to_anchor=(0.5, 0.05), ncol=6, fontsize=12)
# fig.tight_layout(pad=2.0, rect=[0, 0.1, 1, 1])


##plt.savefig('/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 2/gradient_windsonds_down', dpi=300)


cmap = plt.get_cmap('viridis')
colors = cmap(np.linspace(0, 1, len(temp_diff)))
## For up:
labels = ['1430 UTC', '1530 UTC', '1600 UTC', '1800 UTC']


theta_v_sls = theta_v_diff[:4]
    

cmap = plt.get_cmap('viridis')
colors = [cmap(0.0), cmap(0.33), cmap(0.66), cmap(1.0)]

fig, ax = plt.subplots(figsize=(6, 6))

ax.plot(theta_v_sls[0], common_heights, label='1430 UTC', color = cmap(0.0))
ax.plot(theta_v_sls[1], common_heights, label='1530 UTC', color = cmap(0.33))
ax.plot(theta_v_sls[2], common_heights, label='1600 UTC', color = cmap(0.66))
ax.plot(theta_v_sls[3], common_heights, label ='1800 UTC', color = cmap(1.0)) 

## Set our plot bounds
ax.set_ylabel('Height (m)', fontsize = 12)
ax.set_xlabel(r'$\theta_v$ (K)', fontsize = 12)
ax.set_ylim(0, 4500)
##ax.set_xlim(-3,4)
ax.set_title('Theta-v Gradient Between CMAS and SKYLER2', fontsize = 16)
ax.grid(True)



fig.legend(labels, loc='lower center', bbox_to_anchor=(0.5, 0.05), ncol=4, fontsize=11)
fig.tight_layout(pad=2.0, rect=[0, 0.1, 1, 1])

plt.savefig('/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 2/coastal_grad_theta_v_agu_v1.png', dpi=300)
