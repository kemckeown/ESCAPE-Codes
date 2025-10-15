#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 15:54:02 2025

@author: kem6245
"""

import pyart
import fsspec
import matplotlib.pyplot as plt
import os
from io import BytesIO
import warnings
import cartopy.crs as ccrs
from metpy.plots import USCOUNTIES
import cartopy.feature as cfeature
from PIL import Image
from datetime import datetime, timedelta
from IPython.display import display, Image as IPImage
from matplotlib.ticker import FixedLocator
import numpy as np
import matplotlib.ticker as mticker

fs = fsspec.filesystem("s3", anon=True)

start_date = datetime(2022, 6, 2, 10, 0) 
end_date = datetime(2022, 6, 3, 1, 0)  
station = "KHGX"

scan_times = [ 
    datetime(2022, 6, 2, 10, 18, 31), 
    datetime(2022, 6, 2, 12, 41, 48),
    datetime(2022, 6, 2, 14, 34, 29),
    datetime(2022, 6, 2, 14, 52, 1), 
    datetime(2022, 6, 2, 15, 50, 56),
    datetime(2022, 6, 2, 16, 56, 14),
    datetime(2022, 6, 2, 17, 54, 12), 
    datetime(2022, 6, 2, 19, 4, 44),
    datetime(2022, 6, 2, 20, 3, 8)]

latitude = [28.943, 29.336, 29.670]
longitude = [-95.301, -95.687, -95.059]
labels = ["CMAS", "SKYLER2", "TRACER M1"]
markers = ["o", "v", "P"]
colors = ["black", "black","black"]

## Generate the list of files for the specified date and hour range
files = []
current_date = start_date

while current_date <= end_date:
    date_str = current_date.strftime("%Y/%m/%d")
    date_str_compact = current_date.strftime("%Y%m%d")

    if current_date.date() == start_date.date():
        start_hour = start_date.hour
    else:
        start_hour = 0

    if current_date.date() == end_date.date():
        end_hour = end_date.hour
    else:
        end_hour = 23

    for hour in range(start_hour, end_hour + 1):
        hour_str = f"{hour:02d}"
        all_files = fs.glob(
            f"s3://noaa-nexrad-level2/{date_str}/{station}/{station}{date_str_compact}_{hour_str}*"
        )
        filtered_files = [f for f in all_files if not f.endswith("_MDM")]
        files += sorted(filtered_files)

    current_date = current_date.replace(hour=0) + timedelta(days=1)


def find_times(time, files):
    for file in files:
        file_time_str = os.path.basename(file).split("_")[1][0:6]
        file_time = datetime.strptime(file_time_str, "%H%M%S")
        file_time = time.replace(hour=file_time.hour, minute=file_time.minute, second=file_time.second)
        
        if time == file_time:
            return file
        
    return None

our_files = [find_times(time,files) for time in scan_times]


if None in our_files:
    missing_times = [scan_times[i] for i, file in enumerate(our_files) if file is None]
    print(f"Warning: No radar file found for the following times: {missing_times}")
    our_files = [file for file in our_files if file is not None]
    
    
    
def read_radar_data(file_path):
    try:
        with fs.open(file_path, "rb") as f:
            radar_data = f.read()
        radar_file = BytesIO(radar_data)
        radar = pyart.io.read_nexrad_archive(radar_file)
        return radar
    except Exception as e:
        print(f"Failed to read radar data from {file_path}: {e}")
        return None
    
def filter_radar_noise(radar):
    # Create a GateFilter to remove noise
    gatefilter = pyart.correct.GateFilter(radar)
    
    # Remove gates where the reflectivity is less than a certain threshold (e.g., -5 dBZ)
    gatefilter.exclude_below('reflectivity', 10)
    
    # Optionally, remove non-meteorological echoes (speckle filter)
    gatefilter = pyart.correct.despeckle_field(radar, 'reflectivity', gatefilter=gatefilter)
    
    return gatefilter

    
    
fig, axes = plt.subplots(3, 3, figsize=(18, 18), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()  # Flatten the axes array for easier indexing

panel_labels = [f"({chr(97 + i)})" for i in range(9)]  # Generate labels a-i

for i, (file, ax, time) in enumerate(zip(our_files, axes, scan_times)):
    radar = read_radar_data(file)
    if radar is None:
        print(f"Skipping file {file} due to read error.")
        ax.set_axis_off()  # Turn off the axis for missing plots
        continue
    
    gatefilter = filter_radar_noise(radar)

    ## Display radar reflectivity on the current subplot
    radar_display = pyart.graph.RadarMapDisplay(radar)
    try:
        radar_display.plot_ppi_map(
            "reflectivity", sweep=0, vmin=10, vmax=40, ax=ax, 
            cmap="pyart_HomeyerRainbow", colorbar_flag=False,
            gatefilter=gatefilter, min_lat=28, max_lat=31, min_lon=-96.5, max_lon=-93.5,
            lat_lines = np.arange(28, 31, 1), lon_lines = np.arange(-96.5, -93.5, 1), 
            embellish=False)
        ax.set_xlim(-96.5, -93.5)
        ax.set_ylim(28, 31)
        # Access and customize the gridlines
        ax.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="lightgrey", linewidth=0.8)
        
        ax.coastlines('50m', linewidth=1)
    except Exception as e:
        print(f"Error plotting reflectivity for file {file}: {e}")
        ax.set_axis_off()  # Turn off the axis for failed plots
        continue

    ## Add location markers
    for lat, lon, label, marker, color in zip(latitude, longitude, labels, markers, colors):
        ax.plot(lon, lat, marker, label=label, color=color, transform=ccrs.PlateCarree(), markersize=10)

    ## Set plot extent
    
    #radar_display.plot_range_ring(100, linestyle='--', color="black", alpha=0.7)

    ## Add title with panel label and larger font size
    ax.set_title(f"{panel_labels[i]} KHGX {time.strftime('%H:%M:%S')} (UTC)", fontsize=18)

    ## Add legend only to the first subplot
    #if i == 0:
       # ax.legend(loc="upper right", fontsize="medium", title="Locations")

# Turn off unused axes if there are fewer than 9 plots
#for ax in axes[len(our_files):]:
   # ax.set_axis_off()

# Add a shared colorbar for all plots
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # Position of colorbar
cbar = plt.colorbar(radar_display.plots[0], cax=cbar_ax, label="Radar Reflectivity (dBZ)", extend="both")

# Increase font size for colorbar label and ticks
cbar.ax.set_ylabel("Radar Reflectivity (dBZ)", fontsize=18)  # Increase label font size
cbar.ax.tick_params(labelsize=14)  # Increase tick label font size

# Adjust layout to prevent overlap
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Leave space for the colorbar


plt.savefig("/9_panel_radar_paper2_v2.png", bbox_inches="tight", dpi=300)
#plt.close()

#rap_130_20220602_1200_000.grb2



