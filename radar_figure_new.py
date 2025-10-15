#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 14:29:20 2024

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
import time

start_time = time.time()

fs = fsspec.filesystem("s3", anon=True)

start_date = datetime(2022, 6, 2, 14, 8, 17)  # YYYY/MM/DD HH
end_date = datetime(2022, 6, 2, 23, 37, 40)  # YYYY/MM/DD HH
station = "KHGX"

latitude = [28.943, 29.336, 29.328, 29.670]
longitude = [-95.301, -95.687, -95.741, -95.059]
labels = ["CMAS", "SKYLER2", "TRACER S3", "TRACER M1"]
markers = ["o", "o", "o", "o",]
colors = ["black", "hotpink","cyan", "purple"]

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


## Function to read radar data
def read_radar_data(file_path):
    try:
        with fs.open(file_path, "rb") as f:
            radar_data = f.read()
        radar_file = BytesIO(radar_data)
        radar = pyart.io.read_nexrad_archive(radar_file)
        # print(f"Successfully read radar data from {file_path}")
        return radar
    except Exception as e:
        print(f"Failed to read radar data from {file_path}: {e}")
        return None

def filter_radar_noise(radar):
    # Create a GateFilter to remove noise
    gatefilter = pyart.correct.GateFilter(radar)
    
    # Remove gates where the reflectivity is less than a certain threshold (e.g., -5 dBZ)
    gatefilter.exclude_below('reflectivity', 0)
    
    # Optionally, remove non-meteorological echoes (speckle filter)
    gatefilter = pyart.correct.despeckle_field(radar, 'reflectivity', gatefilter=gatefilter)
    
    return gatefilter

## Create directories for the frames
frames_dir = "frames/PPI"
os.makedirs(frames_dir, exist_ok=True)

## Create frames for the animated GIF
frames = []

## Loop through each radar file
for i, file in enumerate(files):
    radar = read_radar_data(file)
    if radar is None:
        print(f"Skipping file {file} due to read error.")
        continue

    # Create a plot for the first sweep's reflectivity
    fig = plt.figure(figsize=[12, 8])
    ax = plt.subplot(111, projection=ccrs.PlateCarree())
    
    gatefilter = filter_radar_noise(radar)
    radar_display = pyart.graph.RadarMapDisplay(radar)
    try:
        radar_display.plot_ppi_map(
            "reflectivity",
            sweep=0,
            vmin=0,
            vmax=60,
            ax=ax,
            title=f"Z for {os.path.basename(file)}",
            cmap="HomeyerRainbow",
            colorbar_flag=False, gatefilter = gatefilter,
        )  # Disable the built-in colorbar
        mappable = radar_display.plots[0]
    except Exception as e:
        print(f"Error plotting radar data for file {file}: {e}")
        plt.close(fig)
        continue

    ## Set the extent to zoom in much closer and centered on the points
    plt.xlim(-96.5, -93.5)
    plt.ylim(28, 31)
    
    
    ## Add counties
    ##ax.add_feature(USCOUNTIES, linewidth=0.5)
    ## Add the coastline without shading the land or water
    ax.coastlines('50m', linewidth=1)
    ax.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="lightgrey", linewidth=0.8)
    radar_display.plot_range_ring(100, linestyle='--', color = "black", alpha = 0.7)
    
    for lat, lon, label, marker, color in zip(
        latitude, longitude, labels, markers, colors
    ):
        ax.plot(
            lon,
            lat,
            marker,
            label=label,
            color=color,
            transform=ccrs.PlateCarree(),
            markersize=6,
        )

    ## Create a colorbar manually
    plt.tight_layout()
    cbar = plt.colorbar(
        mappable, ax=ax, orientation="vertical", fraction=0.046, pad=0.04
    )
    cbar.set_label("equivalent reflectivity factor (dBZ)")

    ## Save the plot to a file
    filename = os.path.join(frames_dir, f"radar_frame_{i}.png")
    plt.legend(loc="upper right", fontsize="large", title="Locations")
    plt.savefig(filename, bbox_inches="tight")
    plt.close(fig)

    ## Add the file to the frames list
    frames.append(filename)


## Create an animated GIF using Pillow
if frames:
    images = [Image.open(frame) for frame in frames]
    gif_filename = "radar_animation.gif"
    images[0].save(
        gif_filename, save_all=True, append_images=images[1:], duration=300, loop=0
    )  # duration in milliseconds

    ## Clean up the saved frames
    for filename in frames:
        os.remove(filename)

    print("Animated GIF created as 'radar_animation.gif'")

    # Display the GIF in the notebook
    with open(gif_filename, "rb") as f:
        display(IPImage(data=f.read(), format="gif"))
else:
    print("No frames were generated.")
    
    
end_time = time.time()  

elapsed_time = end_time - start_time  
print(f"Execution time: {elapsed_time} seconds")


