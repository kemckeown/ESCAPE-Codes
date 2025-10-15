#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 16:02:00 2025

@author: katiemckeown
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
from netCDF4 import Dataset

################################################################################################################
# Read in the TRACER AMF1 data

file_list = [
    '/Volumes/McKeown/Penn State/Penn State Laptop /Python Copy final/ESCAPE/ESCAPE_Data/TRACER_Lidar/houdlfptM1.b1.20220602.100113.cdf',
    '/Volumes/McKeown/Penn State/Penn State Laptop /Python Copy final/ESCAPE/ESCAPE_Data/TRACER_Lidar/houdlfptM1.b1.20220602.110113.cdf',
    '/Volumes/McKeown/Penn State/Penn State Laptop /Python Copy final/ESCAPE/ESCAPE_Data/TRACER_Lidar/houdlfptM1.b1.20220602.120113.cdf',
    '/Volumes/McKeown/Penn State/Penn State Laptop /Python Copy final/ESCAPE/ESCAPE_Data/TRACER_Lidar/houdlfptM1.b1.20220602.130112.cdf',
    '/Volumes/McKeown/Penn State/Penn State Laptop /Python Copy final/ESCAPE/ESCAPE_Data/TRACER_Lidar/houdlfptM1.b1.20220602.140112.cdf',
    '/Volumes/McKeown/Penn State/Penn State Laptop /Python Copy final/ESCAPE/ESCAPE_Data/TRACER_Lidar/houdlfptM1.b1.20220602.150112.cdf',
    '/Volumes/McKeown/Penn State/Penn State Laptop /Python Copy final/ESCAPE/ESCAPE_Data/TRACER_Lidar/houdlfptM1.b1.20220602.160112.cdf',
    '/Volumes/McKeown/Penn State/Penn State Laptop /Python Copy final/ESCAPE/ESCAPE_Data/TRACER_Lidar/houdlfptM1.b1.20220602.170112.cdf',
    '/Volumes/McKeown/Penn State/Penn State Laptop /Python Copy final/ESCAPE/ESCAPE_Data/TRACER_Lidar/houdlfptM1.b1.20220602.180112.cdf'
]

all_times = []
all_dopp_vel = []

range_gate_count = 320
snr_threshold = 0.015  # SNR threshold for filtering

for file_path in file_list: 
    with Dataset(file_path, 'r') as nc: 
        dopp_vel = nc.variables['radial_velocity'][:]
        time_var = nc.variables['time'][:]
        height = nc.variables['range'][:]
        intensity = nc.variables['intensity'][:]
        
        # Filter by SNR
        snr = intensity - 1
        dopp_vel_filtered = np.where(snr < snr_threshold, np.nan, dopp_vel)
        
        # Mask the lowest 100 m to NaN to remove surface bias
        low_alt_mask = height <= 75
        dopp_vel_filtered[:, low_alt_mask] = np.nan

        # --- Apply lidar velocity bias correction ---
        dopp_vel_filtered = dopp_vel_filtered - 0.22  # remove +0.22 m/s bias
        
        # Convert time to pandas datetime
        time_units = nc.variables['time'].units
        time_base = pd.to_datetime(time_units.split('since')[-1].strip())
        times = time_base + pd.to_timedelta(time_var, unit='s')
        
        all_times.append(times)
        all_dopp_vel.append(dopp_vel_filtered)

# Pad arrays so they can be concatenated
max_gates = max(dv.shape[0] for dv in all_dopp_vel)
all_dopp_vel_padded = [np.pad(dv, ((0, max_gates - dv.shape[0]), (0, 0)), 
                              mode='constant', constant_values=np.nan)
                       for dv in all_dopp_vel]
all_times_padded = [np.pad(t, (0, max_gates - t.shape[0]), 
                           mode='constant', constant_values=np.nan)
                    for t in all_times]

all_times = np.concatenate(all_times_padded, axis=0)
all_dopp_vel = np.concatenate(all_dopp_vel_padded, axis=0)

# Select time range
start_time = pd.Timestamp("2022-06-02 10:00:00")
end_time = pd.Timestamp("2022-06-02 18:59:59")
time_mask = (all_times >= start_time) & (all_times <= end_time)
filtered_times = all_times[time_mask]
filtered_dopp_vel = all_dopp_vel[time_mask, :]

# --- Plotting ---
fig, ax = plt.subplots(figsize=(10,6))  # Single subplot

# Copy colormap and set NaNs to white (removes gray edges)
cmap = plt.cm.seismic.copy()
cmap.set_bad(color='white')

im = ax.imshow(
    np.flip(filtered_dopp_vel.T, axis=0),
    extent=[filtered_times.min(), filtered_times.max(),
            height.min(), height.max()],
    aspect='auto',
    cmap=cmap,
    vmin=-2,
    vmax=2
)

# Add color bar
cbar = plt.colorbar(im, ax=ax, orientation='vertical')
cbar.set_label("Doppler Velocity (m/s)", fontsize=12)

ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.set_ylim(0, 2000)
ax.set_ylabel("Height (m)", fontsize=12)
ax.set_xlabel("Time (UTC)", fontsize=12)
ax.set_title('AMF1 Lidar Doppler Velocity', fontsize=16)

plt.tight_layout()
plt.savefig('/Volumes/McKeown/Penn State/Penn State Laptop /Python Copy final/ESCAPE/ESCAPE_Figures/Paper 2/paper2_tracer_lidar_updated.png', dpi=300)
plt.show()