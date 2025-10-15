#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 10:24:07 2024

@author: kem6245
"""

## This map is for the locations of all the data used for part 2 of my phd

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from metpy.plots import USCOUNTIES

fig = plt.figure(figsize=(10, 10))

# Next, we add the projection to the map and add our extent
# Make sure the projection used is consistent across the board
ax = fig.add_subplot(projection=ccrs.PlateCarree())
ax.set_extent([-96.25, -94.25, 28.5, 30.25])
#ax.coastlines()

#state_borders = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines',
                                             # scale='50m', facecolor='none')
#ax.add_feature(state_borders, edgecolor='black', linewidth=1, zorder=3)

ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale='10m', edgecolor='face',
                                      facecolor=cfeature.COLORS['water'])
land = cfeature.NaturalEarthFeature('physical', 'land', scale='10m', edgecolor='face',
                                    facecolor=cfeature.COLORS['land'])

g1 = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='grey')
g1.top_labels = False
g1.left_labels = False
g1.xlines = True
g1.ylines = True
g1.xlocator = mticker.FixedLocator([-94,  -95,  -96])
g1.ylocator = mticker.FixedLocator([28, 29, 30])
g1.yformatter = LATITUDE_FORMATTER
g1.xformatter = LONGITUDE_FORMATTER
g1.xlabel_style = {'size': 16, 'color': 'black'}
g1.ylabel_style = {'size': 16, 'color': 'black'}

ax.add_feature(ocean, zorder=-1)
ax.add_feature(land, zorder=-1)
ax.add_feature(USCOUNTIES.with_scale('500k'), alpha=0.3)

# Choose a location for the scale bar (in lat/lon)
scale_bar_lon = -95.8
scale_bar_lat = 28.6

# Calculate ~0.25 degrees in lon for 25 km at this latitude (approx.)
# 1 degree longitude ~ 111 km * cos(latitude)
import numpy as np
km_per_deg_lon = 111.32 * np.cos(np.radians(scale_bar_lat))
deg_25km = 25 / km_per_deg_lon

# Draw the scale bar
ax.plot([scale_bar_lon, scale_bar_lon + deg_25km],
        [scale_bar_lat, scale_bar_lat],
        color='black', linewidth=4, transform=ccrs.PlateCarree())

# Add text label
ax.text(scale_bar_lon + deg_25km / 2, scale_bar_lat - 0.05,
        '25 km', horizontalalignment='center', verticalalignment='top',
        fontsize=14, color='black', transform=ccrs.PlateCarree())

# ax.plot(-95.688, 29.336, 'o', color ='red', markersize=8, transform=ccrs.PlateCarree(), label='SKYLER2')
# ax.plot(-94.754, 29.827, 'o', color ='red', markersize=8, transform=ccrs.PlateCarree())
# ax.plot(-95.250, 29.123, 'o', color ='red', markersize=8, transform=ccrs.PlateCarree())
# ax.plot(-95.301, 28.943, 's', color = 'blue', markersize=8, transform=ccrs.PlateCarree(), label='CMAS')
# ax.plot(-94.389, 29.550, 's', color = 'blue', markersize=8, transform=ccrs.PlateCarree())
# #ax.plot(-95.059, 29.670, 'o', color = 'cyan', markersize=8, transform=ccrs.PlateCarree(), label='La Porte - M1')  ## Land 1
# #ax.plot(-95.741, 29.328, 'o', color = 'cyan', markersize=8, transform=ccrs.PlateCarree(), label='Guy - S3')  ## Land 2
# ax.plot(-95.422, 29.381, '^', color = 'cyan', markersize=8, transform=ccrs.PlateCarree(), label='PX-1000')
# ax.plot(-95.434, 29.288, '^', color = '', markersize=8, transform=ccrs.PlateCarree())
# ax.plot(-94.502, 29.704, '^', color = 'cyan', markersize=8, transform=ccrs.PlateCarree())
# ax.plot(-95.699, 29.510, 'D', color = 'orange', markersize=8, transform=ccrs.PlateCarree(), label='RaXPol')
# ax.plot(-94.440, 29.997, 'D', color = 'orange', markersize=8, transform=ccrs.PlateCarree())
# ax.plot(-94.792, 29.543, 'D', color = 'orange', markersize=8, transform=ccrs.PlateCarree())
# ax.plot(-95.151, 29.263, 'D', color = 'orange', markersize=8, transform=ccrs.PlateCarree())
# ax.plot(-95.079, 29.472, 'X', color = 'sienna', markersize=8, transform=ccrs.PlateCarree(), label='KHGX')
# ax.plot(-95.6565, 29.6222, 'p', color = 'black', markersize=8, transform=ccrs.PlateCarree(), label='Sugarland Airport')

# ax.plot(-95.688, 29.336, 'o', color ='orange', markersize=10, transform=ccrs.PlateCarree(), label='Southwest Radars')
# ax.plot(-94.754, 29.827, 'o', color ='blue', markersize=10, transform=ccrs.PlateCarree(), label='Southeast Radars')
# ax.plot(-95.250, 29.123, 'o', color ='red', markersize=10, transform=ccrs.PlateCarree(), label='South Radars')
# ax.plot(-95.301, 28.943, 's', color = 'purple', markersize=10, transform=ccrs.PlateCarree(), label='CMAS Sites')
# ax.plot(-94.389, 29.550, 's', color = 'purple', markersize=10, transform=ccrs.PlateCarree())
# #ax.plot(-95.059, 29.670, 'o', color = 'orange', markersize=8, transform=ccrs.PlateCarree(), label='La Porte - M1')  ## Land 1
# #ax.plot(-95.741, 29.328, 'o', color = 'red', markersize=8, transform=ccrs.PlateCarree(), label='Guy - S3')  ## Land 2
# ax.plot(-95.422, 29.381, '^', color = 'orange', markersize=10, transform=ccrs.PlateCarree())
# ax.plot(-95.434, 29.288, '^', color = 'red', markersize=10, transform=ccrs.PlateCarree())
# ax.plot(-94.502, 29.704, '^', color = 'blue', markersize=10, transform=ccrs.PlateCarree())
# ax.plot(-95.699, 29.510, 'D', color = 'orange', markersize=10, transform=ccrs.PlateCarree())
# ax.plot(-94.440, 29.997, 'D', color = 'blue', markersize=10, transform=ccrs.PlateCarree())
# ax.plot(-94.792, 29.543, 'D', color = 'blue', markersize=10, transform=ccrs.PlateCarree())
# ax.plot(-95.151, 29.263, 'D', color = 'red', markersize=10, transform=ccrs.PlateCarree())
#ax.plot(-95.079, 29.472, 'X', color = 'sienna', markersize=10, transform=ccrs.PlateCarree(), label='KHGX')
# ax.plot(-95.6565, 29.6222, 'p', color = 'black', markersize=10, transform=ccrs.PlateCarree(), label='Sugarland Airport')


#ax.plot(-96.2, 30.87, 'o', color = '#440154', markersize=8, transform=ccrs.PlateCarree(), label= 'LHB')
#ax.plot(-96.37, 30.22, 'o', color = '#443a83', markersize=8, transform=ccrs.PlateCarree(), label= '11R')

#ax.plot(-95.66, 29.62, 'o', color = '#35b779', markersize=8, transform=ccrs.PlateCarree(), label= 'SGR')
#ax.plot(-95.46, 29.11, 'o', color = '#fde725', markersize=8, transform=ccrs.PlateCarree(), label= 'LBX')


# ax.plot(-95.079, 29.472, 'o', color = 'black', markersize=8, transform=ccrs.PlateCarree(), label='KHGX')
# ax.plot(-94.860, 29.265, 'o', color = 'plum', markersize=8, transform=ccrs.PlateCarree(), label= 'GLS')
# ax.plot(-95.282, 29.637, 'o', color = 'skyblue', markersize=8, transform=ccrs.PlateCarree(), label= 'HOU')
# ax.plot(-95.462, 29.109, 'o', color = 'purple', markersize=8, transform=ccrs.PlateCarree(), label= 'LBX')
# ax.plot(-95.242, 29.519, 'o', color = 'silver', markersize=8, transform=ccrs.PlateCarree(), label= 'LVJ')
# ax.plot(-95.656, 29.622, 'o', color = 'sienna', markersize=8, transform=ccrs.PlateCarree(), label= 'SGR')
##ax.plot(-95.326, 29.901, 'o', color = 'dimgrey', markersize=8, transform=ccrs.PlateCarree(), label='CLAMPS1')


## For Paper 2

ax.plot(-95.688, 29.336, 'o', color ='cyan', markersize=12, transform=ccrs.PlateCarree(), label='SKYLER-2')
ax.plot(-95.301, 28.943, '^', color = 'blue', markersize=12, transform=ccrs.PlateCarree(), label='CMAS')
ax.plot(-95.059, 29.670, 's', color = 'magenta', markersize=12, transform=ccrs.PlateCarree(), label='AMF1')
ax.plot(-95.741, 29.328, 'd', color = 'red', markersize=12, transform=ccrs.PlateCarree(), label='S3')
ax.plot(-95.079, 29.472, 'X', color = 'black', markersize=15, transform=ccrs.PlateCarree(), label='KHGX')
ax.plot(-95.3698, 29.7604, marker="*", ls = 'None', markersize=20, transform=ccrs.PlateCarree(), label='Downtown Houston')

## Transect 3
#ax.plot(-97.79, 31.42, 's', color = '#013220', markersize=8, transform=ccrs.PlateCarree(), label= 'GOP')
#ax.plot(-97.83, 31.07, 's', color = '#228B22', markersize=8, transform=ccrs.PlateCarree(), label= 'GRK')
#ax.plot(-97.44, 30.57, 's', color = '#008000', markersize=8, transform=ccrs.PlateCarree(), label= 'T74')
#ax.plot(-96.98, 30.17, 's', color = '#32CD32', markersize=8, transform=ccrs.PlateCarree(), label= 'GYB')
#ax.plot(-96.52, 29.64, 's', color = '#90EE90', markersize=8, transform=ccrs.PlateCarree(), label= '66R')
#ax.plot(-96.15, 29.25, 's', color = '#C1E1C1', markersize=8, transform=ccrs.PlateCarree(), label= 'ARM')

## Transect 2
#ax.plot(-96.97, 30.88, '^', color = '#4B0082', markersize=8, transform=ccrs.PlateCarree(), label= 'T35')
#ax.plot(-96.33, 30.72,'^', color = '#6A0DAD', markersize=8, transform=ccrs.PlateCarree(), label= 'CFD')
#ax.plot(-96.37, 30.22, '^', color = '#800080', markersize=8, transform=ccrs.PlateCarree(), label= '11R')
#ax.plot(-95.90, 29.81, '^', color = '#9B30FF', markersize=8, transform=ccrs.PlateCarree(), label= 'TME')
#ax.plot(-95.48, 29.51, '^', color = '#D8BFD8', markersize=8, transform=ccrs.PlateCarree(), label='AXH')



plt.legend(fontsize=16)
ax.legend(loc="lower right", fontsize=12)#, bbox_to_anchor=(1.12, 0.5))
plt.tight_layout()

ax.set_title("2 June 2022 Deployment Sites", fontsize=20)


plt.savefig('updated_site_locations.png', dpi=300)









