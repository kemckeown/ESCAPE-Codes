#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 09:56:41 2025

@author: kem6245
"""

import pyart
import fsspec
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from io import BytesIO
import cartopy.crs as ccrs
from metpy.plots import USCOUNTIES
from datetime import datetime, timedelta
from IPython.display import display, Image as IPImage

fs = fsspec.filesystem("s3", anon=True)

start_date = datetime(2022, 6, 2, 11, 0) 
end_date = datetime(2022, 6, 3, 1, 0)  
station = "KHGX"

scan_times = [ 
    datetime(2022, 6, 2, 20, 8, 0), 
    datetime(2022, 6, 2, 22, 7, 54),
    datetime(2022, 6, 2, 23, 5, 17)]

latitude = [28.943, 29.336, 29.328, 29.670]
longitude = [-95.301, -95.687, -95.741, -95.059]
labels = ["CMAS", "SKYLER2", "TRACER S3", "TRACER M1"]
markers = ["o", "^", "D", "s"]
colors = ["lightsteelblue", "hotpink","lime", "cyan"]

files = []
current_date = start_date

while current_date <= end_date:
    date_str = current_date.strftime("%Y/%m/%d")
    date_str_compact = current_date.strftime("%Y%m%d")

    start_hour = start_date.hour if current_date.date() == start_date.date() else 0
    end_hour = end_date.hour if current_date.date() == end_date.date() else 23

    for hour in range(start_hour, end_hour + 1):
        hour_str = f"{hour:02d}"
        all_files = fs.glob(f"s3://noaa-nexrad-level2/{date_str}/{station}/{station}{date_str_compact}_{hour_str}*")
        files += sorted(f for f in all_files if not f.endswith("_MDM"))

    current_date = current_date.replace(hour=0) + timedelta(days=1)

def find_times(time, files):
    for file in files:
        file_time_str = os.path.basename(file).split("_")[1][0:6]
        file_time = datetime.strptime(file_time_str, "%H%M%S")
        file_time = time.replace(hour=file_time.hour, minute=file_time.minute, second=file_time.second)
        if time == file_time:
            return file
    return None

our_files = [find_times(time, files) for time in scan_times]
our_files = [file for file in our_files if file is not None]

def read_radar_data(file_path):
    try:
        with fs.open(file_path, "rb") as f:
            radar_data = f.read()
        return pyart.io.read_nexrad_archive(BytesIO(radar_data))
    except Exception as e:
        print(f"Failed to read radar data from {file_path}: {e}")
        return None

def filter_radar_noise(radar):
    gatefilter = pyart.correct.GateFilter(radar)
    gatefilter.exclude_below('reflectivity', 0)
    return pyart.correct.despeckle_field(radar, 'reflectivity', gatefilter=gatefilter)

fig, axes = plt.subplots(1, 3, figsize=(18, 6), subplot_kw={'projection': ccrs.PlateCarree()})
mappable = None  

for i, (file, ax, time) in enumerate(zip(our_files, axes, scan_times)):
    radar = read_radar_data(file)
    if radar is None:
        print(f"Skipping file {file} due to read error.")
        continue
    
    gatefilter = filter_radar_noise(radar)
    ax.set_facecolor("black")
    radar_display = pyart.graph.RadarMapDisplay(radar)
    
    mappable = radar_display.plot_ppi_map(
        "reflectivity", sweep=0, vmin=0, vmax=40, ax=ax, 
        title=f"KHGX {time.strftime('%H:%M:%S')} (UTC)", cmap="afmhot", fontsize=14, colorbar_flag=False,
        gatefilter=gatefilter
    )
    
    ax.coastlines('50m', linewidth=1)
    
    for lat, lon, label, marker, color in zip(latitude, longitude, labels, markers, colors):
        ax.plot(lon, lat, marker, label=label, color=color, transform=ccrs.PlateCarree(), markersize=6)

    ax.set_xlim(-96, -94.5)
    ax.set_ylim(28.5, 30)
    gridlines = ax.gridlines(draw_labels=True, xlocs=[-96, -95.5, -95, -94.5], ylocs=[28.5, 29, 29.5, 30], color='gray', linewidth=0.5)
    gridlines.top_labels = False
    gridlines.right_labels = False
    gridlines.left_labels = i % 3 == 0
    gridlines.bottom_labels = i // 3 == 0

    if i == 0:
        ax.legend(loc="lower right", fontsize="medium", title="Locations")

    # Access reflectivity field directly from radar data
    reflectivity = radar.fields['reflectivity']['data']
    im = ax.imshow(reflectivity, cmap="afmhot", vmin=0, vmax=40)

        # Add colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # Position of colorbar
    plt.colorbar(im, cax=cbar_ax, label="Radar Reflectivity (dBZ)", fontsize=14)




plt.savefig("/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 2/late_period_radar_v5_zoomed_filter_0.pdf", bbox_inches="tight", dpi=300)








