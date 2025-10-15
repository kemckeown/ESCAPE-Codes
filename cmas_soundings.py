#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 17:04:42 2024

@author: kem6245
"""

## CMAS June 2 Soundings

## Load in our modules
import xarray as xr
import matplotlib.pyplot as plt
from metpy.plots import SkewT, Hodograph
from metpy.units import units
import metpy.calc as mpcalc
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from metpy.calc import bulk_shear, wind_components
import numpy as np
import pandas as pd
import metpy 
import netCDF4 as nc
from scipy import interpolate
import scipy.interpolate as interp


def cmas_filtering(month_str,day_str,time_str,p,t,td,u,v, wind_speed, h):
    month = int(month_str)
    day = int(day_str)
    time = int(time_str)
    
    ## June 2nd 2022 IOP 3 (Surfside Beach)
    if month == 6 and day == 2 and time == 145623:
        return p[17:], t[17:], td[17:], u[17:], v[17:], wind_speed[17:], h[17:]
    elif month == 6 and day == 2 and time == 170653:
        return p[30:2431], t[30:2431], td[30:2431], u[30:2431], v[30:2431], wind_speed[30:2431], h[30:2431]
    elif month == 6 and day == 2 and time == 184909:
        return p[16:2184], t[16:2184], td[16:2184], u[16:2184], v[16:2184], wind_speed[16:2184], h[16:2184]
    elif month == 6 and day == 2 and time == 204526:
        return p[35:], t[35:], td[35:], u[35:], v[35:], wind_speed[35:], h[35:]
    elif month == 6 and day == 2 and time == 225508:
        return p[16:1861], t[16:1861], td[16:1861], u[16:1861], v[16:1861], wind_speed[16:1861], h[16:1861]
    
    else:
        
        print('Please check your sounding info. You must have made a mistake.')
        
cmas_radiosondes = '2_June_2022'

if not os.path.exists(cmas_radiosondes):
    print("❌ Directory does not exist:", cmas_radiosondes)
else:
    print("✅ Directory exists!")

# List all files in the directory (for debugging)
print("📂 Files in directory:", os.listdir(cmas_radiosondes))        

for file in os.scandir(cmas_radiosondes):
    if file.is_file() and file.name.endswith('.nc'):  # Ensure it's a .nc file
        file_name = file.path  # Get full file path
        print("📄 Found NetCDF file:", file_name)

    
    ## Read in the data and convert all of the data into arrays
        data = xr.open_dataset(file_name)
    
        p_in = data['pressure'].to_numpy()
        p_real = data['pressure'].values * units.hPa
        t_in = data['temperature'].to_numpy()
        td_in = data['dewpoint_temperature'].to_numpy()
        h_in = data['geometric_height'].values * units.m
        wind_speed_in = data['wind_speed'].to_numpy()
        windspeed = data['wind_speed'].values * (units.meter / units.second)
        wind_dir_in = data['wind_direction'].to_numpy()
        winddir = data['wind_direction'].values * units.degrees
        u_in, v_in = mpcalc.wind_components((wind_speed_in * (units.meter/units.second)).to(units.knots), (wind_dir_in * units.degrees));

        p, t, td, u, v, wind_speed, h = cmas_filtering(file_name[-14:-12],file_name[-12:-10],file_name[-9:-3],p_in,t_in,td_in,u_in,v_in, wind_speed_in, h_in)

     ## We need to add units to everything
    ## now units need to be added to each of the above
        p = units.hPa * p

    ## Temperature
        t = units.degC * t

    ## Dew Point
        td = units.degC * td
    
        h = units.m * h

    # Calculate the parcel profile.
   ## parcel_prof = mpcalc.parcel_profile(p, t[0], td[0]).to('degC')
    
        u1, v1 = wind_components(windspeed, winddir)
        print(u1)

        f = interpolate.interp1d(p, h)
    
    # # Filter data for pressures below 500 hPa
    # mask = p >= 500 *units.hPa
    # p1 = p[mask]
    # t1 = t[mask]
    # td1 = td[mask]
    # h1 = h[mask]
    
    # pw = mpcalc.precipitable_water(p1, td1)
    
    # print(pw)
    
    # Compute bulk shear
    # ubshr1, vbshr1 = mpcalc.bulk_shear(p_real, u1, v1, depth=1 * units.km)
    # bshear1 = mpcalc.wind_speed(ubshr1, vbshr1)
    # ubshr3, vbshr3 = mpcalc.bulk_shear(p_real, u1, v1, depth=3 * units.km)
    # bshear3 = mpcalc.wind_speed(ubshr3, vbshr3)
    # ubshr6, vbshr6 = mpcalc.bulk_shear(p_real, u1, v1, depth=6 * units.km)
    # bshear6 = mpcalc.wind_speed(ubshr6, vbshr6)
    
    # print("0-1 km Shear: ", bshear1)
    # print("0-3 km Shear: ", bshear3)
    # print("0-6 km Shear: ", bshear6)
    # print()
    # print()
    
    

        ##Change default to be better for skew-T
        fig = plt.figure(figsize=(6,6))
        skew = SkewT(fig,rotation=45)
        print('hi')

    ##Plot the data using normal plotting functions, in this case using
    ##log scaling in Y, as dictated by the typical meteorological plot
        skew.plot(p, t, 'r', linewidth=2)
        skew.plot(p, td, 'g', linewidth=2)
   
        skew.plot_barbs(p[::100], u[::100], v[::100])
    




    ##Add the relevant special lines
        skew.plot_dry_adiabats(t0=np.arange(233,533,10)*units('K'),alpha=0.25)
        skew.plot_moist_adiabats(t0=np.arange(233,323,5)*units('K'),alpha=0.25)
        skew.plot_mixing_lines(alpha=0.25)
        skew.ax.set_ylim(1030, 100)
        skew.ax.set_xlim(-20,40)
        skew.ax.set_xlabel('Temperature (\N{DEGREE SIGN}C)', fontsize = 14)
        skew.ax.set_ylabel('Pressure (hPa)', fontsize = 14)

    ##Calculate full parcel profile and add to plot as black line
    ##Requires that the variables have associated units
        prof = mpcalc.parcel_profile(p, t[0], td[0]).to('degC')
        skew.plot(p, prof, 'k', linewidth=3)


    ##calcualte LCL
        lcl_p,lcl_t=mpcalc.lcl(p[0],t[0],td[0]) 
    ## out_parameter.write(str(np.around(lcl_p.m))+'\n')
        print('LCL:',np.around(lcl_p))
        lcl_height = f(lcl_p)
        print('LCL:', np.around(lcl_height))

    ##calculate LFC
        lfc_p,lfc_t=mpcalc.lfc(p,t,td,prof)
    ## out_parameter.write(str(np.around(lfc_p.m))+'\n')
        print('LFC:',np.around(lfc_p))
        lfc_height = f(lfc_p)
        print('LFC:', np.around(lfc_height))

    
        skew.ax.set_title('CMAS Sounding '  + ' ' + file_name[-9:-7] + ':' + file_name[-7:-5])

    #plt.savefig('/CMAS_Sounding'+'_'+file_name[-9:-5], dpi=300)

     ## CAPE and CIN 
        cape, cin=mpcalc.surface_based_cape_cin(p,t,td)
        print(cape, cin)
    
        print('DONE')
        print()
        print()
        print()
    




# ## Let's try difference plots 

# # Define the common height grid
# common_height = np.arange(0, 10000, 20) * units.m  # Adjust the range and step size as needed

# def interpolate_to_common_height(height, temperature, common_height):
#     # Ensure height and temperature are in numpy arrays
#     height = np.asarray(height)
#     temperature = np.asarray(temperature)

#     # Create interpolation function
#     interp_func = interp.interp1d(height, temperature, kind='linear', fill_value='extrapolate')
    
#     # Interpolate to common height grid
#     return interp_func(common_height)

# # Extract and interpolate data from each file
# interpolated_data = []
# for file_name in os.scandir(cmas_radiosondes):
#     data = xr.open_dataset(file_name.path)
#     p_in = data['pressure'].to_numpy()
#     t_in = data['temperature'].to_numpy()
#     h_in = data['geometric_height'].values * units.m

#     p, t, td, u, v, wind_speed, h = cmas_filtering(
#         file_name.name[-14:-12], file_name.name[-12:-10], file_name.name[-9:-3],
#         p_in, t_in, td_in, u_in, v_in, wind_speed_in, h_in
#     )
    
#     # Interpolate temperature to common height
#     interpolated_temperature = interpolate_to_common_height(h, t, common_height)
#     interpolated_data.append(interpolated_temperature)
    
    
# # Compute differences between consecutive soundings
# temperature_differences = [interpolated_data[i + 1] - interpolated_data[i] for i in range(len(interpolated_data) - 1)]



# fig, ax = plt.subplots(figsize=(10, 6))

# # Plot each difference
# for i, temp_diff in enumerate(temperature_differences):
#     ax.plot(temp_diff.magnitude, common_height.magnitude, label=f'Difference {i+1}-{i+2}')

# # Add labels and title
# ax.set_xlabel('Temperature Difference (°C)')
# ax.set_ylabel('Height (m)')
# ax.set_title('Temperature Differences with Height')
# ax.legend()
# ax.grid(True)

# Save the plot
##plt.savefig('temperature_differences.png', dpi=300)




    
