#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 10:07:10 2024

@author: kem6245
"""

## The purpose of this code is to calculate the gradients of theta-v, theta-e
## and theta-es for the Land 1 - Gulf 1 and Land 2 - Gulf 2 stationary sites.
## This is part of the ESCAPE/TRACER datasets. 

## Load in our modules
import xarray as xr
import numpy as np
import os 
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
from metpy.plots import SkewT,Hodograph
import pandas as pd
from metpy.calc import mixing_ratio_from_relative_humidity
from metpy.calc import relative_humidity_from_dewpoint
from metpy.calc import virtual_potential_temperature
import csv
from scipy.stats import bootstrap
from array import *

#############################################################################################
## Next, we need our directories 

## Land Regime Directories
ll1 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Land_Data/CSV Files/Land 1"
ll2 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Land_Data/CSV Files/Land 2"
lg1 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Land_Data/CSV Files/Gulf 1"
lg2 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Land_Data/CSV Files/Gulf 2"

## Gulf Regime Directories
gl1 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Gulf_Data/CSV Files/Land 1"
gl2 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Gulf_Data/CSV Files/Land 2"
gg1 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Gulf_Data/CSV Files/Gulf 1"
gg2 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Gulf_Data/CSV Files/Gulf 2"

## DWL Regime Directories
b1l1 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 1 RAP Comp/CSV/Land 1"
b1l2 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 1 RAP Comp/CSV/Land 2"
b1g1 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 1 RAP Comp/CSV/Gulf 1"
b1g2 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 1 RAP Comp/CSV/Gulf 2"

## DWL Regime Directories
b2l1 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 2 RAP Comp/CSV/Land 1"
b2l2 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 2 RAP Comp/CSV/Land 2"
b2g1 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 2 RAP Comp/CSV/Gulf 1"
b2g2 = "/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 2 RAP Comp/CSV/Gulf 2"

###############################################################################################
## Next, we need to preallocate all of our arrays. We only need to do this for the variables
## we want, so temp, dew point temp, pressure, height, theta e, theta v, and theta es. 

## Start with Land 1 for our Land Regime
ll1_p = np.ones((37,11))*np.nan
ll1_h = np.ones((37,11))*np.nan
ll1_te = np.ones((37,11))*np.nan
ll1_tv = np.ones((37,11))*np.nan
ll1_tes = np.ones((37,11))*np.nan

##Land 2 for our Land Regime 
ll2_p = np.ones((37,11))*np.nan
ll2_h = np.ones((37,11))*np.nan
ll2_te = np.ones((37,11))*np.nan
ll2_tv = np.ones((37,11))*np.nan
ll2_tes = np.ones((37,11))*np.nan

##Gulf 1 for our Land Regime
lg1_p = np.ones((37,11))*np.nan
lg1_h = np.ones((37,11))*np.nan
lg1_te = np.ones((37,11))*np.nan
lg1_tv = np.ones((37,11))*np.nan
lg1_tes = np.ones((37,11))*np.nan

##Gulf 2 for our Land Regime
lg2_p = np.ones((37,11))*np.nan
lg2_h = np.ones((37,11))*np.nan
lg2_te = np.ones((37,11))*np.nan
lg2_tv = np.ones((37,11))*np.nan
lg2_tes = np.ones((37,11))*np.nan


## Start with Land 1 for our Gulf Regime
gl1_p = np.ones((37,12))*np.nan
gl1_h = np.ones((37,12))*np.nan
gl1_te = np.ones((37,12))*np.nan
gl1_tv = np.ones((37,12))*np.nan
gl1_tes = np.ones((37,12))*np.nan

##Land 2 for our Gulf Regime
gl2_p = np.ones((37,12))*np.nan
gl2_h = np.ones((37,12))*np.nan
gl2_te = np.ones((37,12))*np.nan
gl2_tv = np.ones((37,12))*np.nan
gl2_tes = np.ones((37,12))*np.nan

##Gulf 1 for our Gulf Regime
gg1_p = np.ones((37,12))*np.nan
gg1_h = np.ones((37,12))*np.nan
gg1_te = np.ones((37,12))*np.nan
gg1_tv = np.ones((37,12))*np.nan
gg1_tes = np.ones((37,12))*np.nan

##Gulf 2 for our Gulf Regime
gg2_p = np.ones((37,12))*np.nan
gg2_h = np.ones((37,12))*np.nan
gg2_te = np.ones((37,12))*np.nan
gg2_tv = np.ones((37,12))*np.nan
gg2_tes = np.ones((37,12))*np.nan


## Start with Land 1 for our DWL Regime
b1l1_p = np.ones((37,15))*np.nan
b1l1_h = np.ones((37,15))*np.nan
b1l1_te = np.ones((37,15))*np.nan
b1l1_tv = np.ones((37,15))*np.nan
b1l1_tes = np.ones((37,15))*np.nan

##Land 2 for our DWL Regime
b1l2_p = np.ones((37,15))*np.nan
b1l2_h = np.ones((37,15))*np.nan
b1l2_te = np.ones((37,15))*np.nan
b1l2_tv = np.ones((37,15))*np.nan
b1l2_tes = np.ones((37,15))*np.nan

##Gulf 1 for our DWL Regime
b1g1_p = np.ones((37,15))*np.nan
b1g1_h = np.ones((37,15))*np.nan
b1g1_te = np.ones((37,15))*np.nan
b1g1_tv = np.ones((37,15))*np.nan
b1g1_tes = np.ones((37,15))*np.nan

##Gulf 2 for our DWL Regime
b1g2_p = np.ones((37,15))*np.nan
b1g2_h = np.ones((37,15))*np.nan
b1g2_te = np.ones((37,15))*np.nan
b1g2_tv = np.ones((37,15))*np.nan
b1g2_tes = np.ones((37,15))*np.nan


## Start with Land 1 for our DWG Regime
b2l1_p = np.ones((37,15))*np.nan
b2l1_h = np.ones((37,15))*np.nan
b2l1_te = np.ones((37,15))*np.nan
b2l1_tv = np.ones((37,15))*np.nan
b2l1_tes = np.ones((37,15))*np.nan

##Land 2 for our DWG Regime
b2l2_p = np.ones((37,15))*np.nan
b2l2_h = np.ones((37,15))*np.nan
b2l2_te = np.ones((37,15))*np.nan
b2l2_tv = np.ones((37,15))*np.nan
b2l2_tes = np.ones((37,15))*np.nan

##Gulf 1 for our DWG Regime
b2g1_p = np.ones((37,15))*np.nan
b2g1_h = np.ones((37,15))*np.nan
b2g1_te = np.ones((37,15))*np.nan
b2g1_tv = np.ones((37,15))*np.nan
b2g1_tes = np.ones((37,15))*np.nan

##Gulf 2 for our DWG Regime
b2g2_p = np.ones((37,15))*np.nan
b2g2_h = np.ones((37,15))*np.nan
b2g2_te = np.ones((37,15))*np.nan
b2g2_tv = np.ones((37,15))*np.nan
b2g2_tes = np.ones((37,15))*np.nan

#############################################################################

## Now that our arrays have been preallocated, we can read in the data and add 
## everything to the variables

## First we will do this for our Land 1 site for our land regime.

## We create a file coult to index
file_count = 0

## We want to go through the entire directory file by file
for file in os.scandir(ll1):
    file_name = file.path
    print(file_name)
    
    
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Land_Data/CSV Files/Land 1/.DS_Store':
        continue
    
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
    
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
    
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
    
    
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
    
    ll1_p[:, file_count] = p
    ll1_h[:, file_count] = h
    ll1_tv[:, file_count] = tv
    ll1_te[:, file_count] = te
    ll1_tes[:, file_count] = tes
    
    ## We are going to create a function to do the math for theta e and theta es 

    file_count = file_count + 1

################################################################################
## We will repeat the above for each thing we need.
## Land 2 Land Regime
file_count = 0
    
for file in os.scandir(ll2):
    file_name = file.path
    print(file_name)
    
    
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Land_Data/CSV Files/Land 2/.DS_Store':
        continue
    
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
    
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
    
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
    
    
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
    
    ll2_p[:, file_count] = p
    ll2_h[:, file_count] = h
    ll2_tv[:, file_count] = tv
    ll2_te[:, file_count] = te
    ll2_tes[:, file_count] = tes
    
    file_count = file_count + 1

###################################################################################
## Gulf 1 Land Regime
file_count = 0
for file in os.scandir(lg1):
    file_name = file.path
    print(file_name)
    
    
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Land_Data/CSV Files/Gulf 1/.DS_Store':
        continue
    
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
    
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
    
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
    
    
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
    
    lg1_p[:, file_count] = p
    lg1_h[:, file_count] = h
    lg1_tv[:, file_count] = tv
    lg1_te[:, file_count] = te
    lg1_tes[:, file_count] = tes
    
    file_count = file_count + 1
    
#######################################################################################
## Gulf 2 Land Regime
file_count = 0
for file in os.scandir(lg2):
    file_name = file.path
    print(file_name)
    
    
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Land_Data/CSV Files/Gulf 2/.DS_Store':
        continue
    
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
    
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
    
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
    
    
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
    
    lg2_p[:, file_count] = p
    lg2_h[:, file_count] = h
    lg2_tv[:, file_count] = tv
    lg2_te[:, file_count] = te
    lg2_tes[:, file_count] = tes
    
    file_count = file_count + 1
    
  ##################################################################################  
## First we will do this for our Land 1 site for our Gulf regime.

## We create a file coult to index
file_count = 0

## We want to go through the entire directory file by file
for file in os.scandir(gl1):
    file_name = file.path
    print(file_name)
        
        
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Gulf_Data/CSV Files/Land 1/.DS_Store':
        continue
        
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
        
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
        
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
        
        
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
        
    gl1_p[:, file_count] = p
    gl1_h[:, file_count] = h
    gl1_tv[:, file_count] = tv
    gl1_te[:, file_count] = te
    gl1_tes[:, file_count] = tes
        
    file_count = file_count + 1

################################################################################
## We will repeat the above for each thing we need.
## Land 2 Gulf Regime
file_count = 0
        
for file in os.scandir(gl2):
    file_name = file.path
    print(file_name)
        
        
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Gulf_Data/CSV Files/Land 2/.DS_Store':
        continue
    
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
        
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
        
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
        
        
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
        
    gl2_p[:, file_count] = p
    gl2_h[:, file_count] = h
    gl2_tv[:, file_count] = tv
    gl2_te[:, file_count] = te
    gl2_tes[:, file_count] = tes
        
    file_count = file_count + 1

###################################################################################
 ## Gulf 1 Gulf Regime
file_count = 0
for file in os.scandir(gg1):
    file_name = file.path
    print(file_name)
        
        
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Gulf_Data/CSV Files/Gulf 1/.DS_Store':
        continue
        
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
        
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
        
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
        
        
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
        
    gg1_p[:, file_count] = p
    gg1_h[:, file_count] = h
    gg1_tv[:, file_count] = tv
    gg1_te[:, file_count] = te
    gg1_tes[:, file_count] = tes
        
    file_count = file_count + 1
        
#######################################################################################
## Gulf 2 Gulf Regime
file_count = 0
for file in os.scandir(gg2):
    file_name = file.path
    print(file_name)
        
        
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Gulf_Data/CSV Files/Gulf 2/.DS_Store':
        continue
        
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
        
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
        
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
        
        
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
        
    gg2_p[:, file_count] = p
    gg2_h[:, file_count] = h
    gg2_tv[:, file_count] = tv
    gg2_te[:, file_count] = te
    gg2_tes[:, file_count] = tes
        
    file_count = file_count + 1
    
##########################################################################################

## First we will do this for our Land 1 site for our DWL regime.

## We create a file coult to index
file_count = 0

## We want to go through the entire directory file by file
for file in os.scandir(b1l1):
    file_name = file.path
    print(file_name)
        
        
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 1 RAP Comp/CSV/Land 1/.DS_Store':
        continue
        
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
        
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
        
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
    
        
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
    
    b1l1_p[:, file_count] = p
    b1l1_h[:, file_count] = h
    b1l1_tv[:, file_count] = tv
    b1l1_te[:, file_count] = te
    b1l1_tes[:, file_count] = tes
        
    file_count = file_count + 1

 ################################################################################
## We will repeat the above for each thing we need.
## Land 2 DWL Regime
file_count = 0
        
for file in os.scandir(b1l2):
    file_name = file.path
    print(file_name)
        
        
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 1 RAP Comp/CSV/Land 2/.DS_Store':
        continue
        
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
        
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
        
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
    
        
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points.
        ## We need to store pressure, height, tv, te, tes.
    b1l2_p[:, file_count] = p
    b1l2_h[:, file_count] = h
    b1l2_tv[:, file_count] = tv
    b1l2_te[:, file_count] = te
    b1l2_tes[:, file_count] = tes
    
    file_count = file_count + 1

###################################################################################
## Gulf 1 DWL Regime
file_count = 0
for file in os.scandir(b1g1):
    file_name = file.path
    print(file_name)
        
        
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 1 RAP Comp/CSV/Gulf 1/.DS_Store':
        continue
    
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
    
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
    
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
        
        
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
        
    b1g1_p[:, file_count] = p
    b1g1_h[:, file_count] = h
    b1g1_tv[:, file_count] = tv
    b1g1_te[:, file_count] = te
    b1g1_tes[:, file_count] = tes
    
    file_count = file_count + 1
        
#######################################################################################
## Gulf 2 DWL Regime
file_count = 0
for file in os.scandir(b1g2):
    file_name = file.path
    print(file_name)
        
        
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 1 RAP Comp/CSV/Gulf 2/.DS_Store':
        continue
        
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
    
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
        
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
        
        
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
    
    b1g2_p[:, file_count] = p
    b1g2_h[:, file_count] = h
    b1g2_tv[:, file_count] = tv
    b1g2_te[:, file_count] = te
    b1g2_tes[:, file_count] = tes
        
    file_count = file_count + 1
##########################################################################################

##########################################################################################

## First we will do this for our Land 1 site for our DWG regime.

## We create a file coult to index
file_count = 0

## We want to go through the entire directory file by file
for file in os.scandir(b2l1):
    file_name = file.path
    print(file_name)
        
        
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 2 RAP Comp/CSV/Land 1/.DS_Store':
        continue
        
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
        
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
        
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
    
        
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
    
    b2l1_p[:, file_count] = p
    b2l1_h[:, file_count] = h
    b2l1_tv[:, file_count] = tv
    b2l1_te[:, file_count] = te
    b2l1_tes[:, file_count] = tes
        
    file_count = file_count + 1

 ################################################################################
## We will repeat the above for each thing we need.
## Land 2 DWG Regime
file_count = 0
        
for file in os.scandir(b2l2):
    file_name = file.path
    print(file_name)
        
        
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 2 RAP Comp/CSV/Land 2/.DS_Store':
        continue
        
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
        
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
        
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
    
        
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points.
        ## We need to store pressure, height, tv, te, tes.
    b2l2_p[:, file_count] = p
    b2l2_h[:, file_count] = h
    b2l2_tv[:, file_count] = tv
    b2l2_te[:, file_count] = te
    b2l2_tes[:, file_count] = tes
    
    print("file done")
    
    file_count = file_count + 1

###################################################################################
## Gulf 1 DWG Regime
file_count = 0
for file in os.scandir(b2g1):
    file_name = file.path
    print(file_name)
        
        
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 2 RAP Comp/CSV/Gulf 1/.DS_Store':
        continue
    
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
    
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
    
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
        
        
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
        
    b2g1_p[:, file_count] = p
    b2g1_h[:, file_count] = h
    b2g1_tv[:, file_count] = tv
    b2g1_te[:, file_count] = te
    b2g1_tes[:, file_count] = tes
    
    file_count = file_count + 1
        
#######################################################################################
## Gulf 2 DWG Regime
file_count = 0
for file in os.scandir(b2g2):
    file_name = file.path
    print(file_name)
        
        
    ## Skip the file that's causing the issue
    if file_name == '/Users/kem6245/Documents/Python Copy/ESCAPE/RAP_Data/Both Mode 2 RAP Comp/CSV/Gulf 2/.DS_Store':
        continue
        
    ## Read in the data
    ws1 = pd.read_csv(file_name, sep='\t')
    
    ## Sort into the variables we want
    p = ws1['Pressure'].values * units.hPa
    t = ws1['Temp'].values * units.degC
    td = ws1['Dew'].values * units.degC
    h = ws1['Height'].values * units.m
        
    ## Now do our calculations based on the above variables
    rh = relative_humidity_from_dewpoint(t, td)
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    tv = virtual_potential_temperature(p, t, mr)
    te = mpcalc.equivalent_potential_temperature(p, t, td)
    tes = mpcalc.saturation_equivalent_potential_temperature(p, t)
        
        
    ## Now we need to put this into our variables to store. 
    ## We don't need to store temperature and dew points. 
    ## We need to store pressure, height, tv, te, tes.
    
    b2g2_p[:, file_count] = p
    b2g2_h[:, file_count] = h
    b2g2_tv[:, file_count] = tv
    b2g2_te[:, file_count] = te
    b2g2_tes[:, file_count] = tes
        
    file_count = file_count + 1
##########################################################################################

## Lovely. 
## Now that we have all of our variables, we can do math to get our gradients
## This will significantly decrease the number of things we will be plotting
## Thank goodness. 


## We need to do interpolations for our theta-v at each height

## Land 

## Land 1
from metpy.interpolate import interpolate_1d
inh_ll1 = np.arange(0,5000,100)
inhtv_ll1 = np.zeros((inh_ll1.size, ll1_h.shape[1]))

for i in range(ll1_h.shape[1]):
    inhtv_ll1[:,i] = interpolate_1d(inh_ll1, ll1_h[:,i], ll1_tv[:,i])
    
## Land 2
inh_ll2 = np.arange(0,5000,100)
inhtv_ll2 = np.zeros((inh_ll2.size, ll2_h.shape[1]))

for i in range(ll2_h.shape[1]):
    inhtv_ll2[:,i] = interpolate_1d(inh_ll2, ll2_h[:,i], ll2_tv[:,i])

## Gulf 1
inh_lg1 = np.arange(0,5000,100)
inhtv_lg1 = np.zeros((inh_lg1.size, lg1_h.shape[1]))

for i in range(lg1_h.shape[1]):
    inhtv_lg1[:,i] = interpolate_1d(inh_lg1, lg1_h[:,i], lg1_tv[:,i])

## Gulf 2
inh_lg2 = np.arange(0,5000,100)
inhtv_lg2 = np.zeros((inh_lg2.size, lg2_h.shape[1]))

for i in range(ll1_h.shape[1]):
    inhtv_lg2[:,i] = interpolate_1d(inh_lg2, lg2_h[:,i], lg2_tv[:,i])
    
## Gulf

## Land 1
from metpy.interpolate import interpolate_1d
inh_gl1 = np.arange(0,5000,100)
inhtv_gl1 = np.zeros((inh_gl1.size, gl1_h.shape[1]))

for i in range(gl1_h.shape[1]):
    inhtv_gl1[:,i] = interpolate_1d(inh_gl1, gl1_h[:,i], gl1_tv[:,i])
    
## Land 2
inh_gl2 = np.arange(0,5000,100)
inhtv_gl2 = np.zeros((inh_gl2.size, gl2_h.shape[1]))

for i in range(gl2_h.shape[1]):
    inhtv_gl2[:,i] = interpolate_1d(inh_gl2, gl2_h[:,i], gl2_tv[:,i])

## Gulf 1
inh_gg1 = np.arange(0,5000,100)
inhtv_gg1 = np.zeros((inh_gg1.size, gg1_h.shape[1]))

for i in range(gg1_h.shape[1]):
    inhtv_gg1[:,i] = interpolate_1d(inh_gg1, gg1_h[:,i], gg1_tv[:,i])

## Gulf 2
inh_gg2 = np.arange(0,5000,100)
inhtv_gg2 = np.zeros((inh_gg2.size, gg2_h.shape[1]))

for i in range(gl1_h.shape[1]):
    inhtv_gg2[:,i] = interpolate_1d(inh_gg2, gg2_h[:,i], gg2_tv[:,i])
    
## DWL

## Land 1
from metpy.interpolate import interpolate_1d
inh_b1l1 = np.arange(0,5000,100)
inhtv_b1l1 = np.zeros((inh_b1l1.size, b1l1_h.shape[1]))

for i in range(b1l1_h.shape[1]):
    inhtv_b1l1[:,i] = interpolate_1d(inh_b1l1, b1l1_h[:,i], b1l1_tv[:,i])
    
## Land 2
inh_b1l2 = np.arange(0,5000,100)
inhtv_b1l2 = np.zeros((inh_b1l2.size, b1l2_h.shape[1]))

for i in range(b1l2_h.shape[1]):
    inhtv_b1l2[:,i] = interpolate_1d(inh_b1l2, b1l2_h[:,i], b1l2_tv[:,i])

## Gulf 1
inh_b1g1 = np.arange(0,5000,100)
inhtv_b1g1 = np.zeros((inh_b1g1.size, b1g1_h.shape[1]))

for i in range(b1g1_h.shape[1]):
    inhtv_b1g1[:,i] = interpolate_1d(inh_b1g1, b1g1_h[:,i], b1g1_tv[:,i])

## Gulf 2
inh_b1g2 = np.arange(0,5000,100)
inhtv_b1g2 = np.zeros((inh_b1g2.size, b1g2_h.shape[1]))

for i in range(b1l1_h.shape[1]):
    inhtv_b1g2[:,i] = interpolate_1d(inh_b1g2, b1g2_h[:,i], b1g2_tv[:,i])
    
## DWG

## Land 1
from metpy.interpolate import interpolate_1d
inh_b2l1 = np.arange(0,5000,100)
inhtv_b2l1 = np.zeros((inh_b2l1.size, b2l1_h.shape[1]))

for i in range(b2l1_h.shape[1]):
    inhtv_b2l1[:,i] = interpolate_1d(inh_b2l1, b2l1_h[:,i], b2l1_tv[:,i])
    
## Land 2
inh_b2l2 = np.arange(0,5000,100)
inhtv_b2l2 = np.zeros((inh_b2l2.size, b2l2_h.shape[1]))

for i in range(b2l2_h.shape[1]):
    inhtv_b2l2[:,i] = interpolate_1d(inh_b2l2, b2l2_h[:,i], b2l2_tv[:,i])

## Gulf 1
inh_b2g1 = np.arange(0,5000,100)
inhtv_b2g1 = np.zeros((inh_b2g1.size, b2g1_h.shape[1]))

for i in range(b2g1_h.shape[1]):
    inhtv_b2g1[:,i] = interpolate_1d(inh_b2g1, b2g1_h[:,i], b2g1_tv[:,i])

## Gulf 2
inh_b2g2 = np.arange(0,5000,100)
inhtv_b2g2 = np.zeros((inh_b2g2.size, b2g2_h.shape[1]))

for i in range(b2l1_h.shape[1]):
    inhtv_b2g2[:,i] = interpolate_1d(inh_b2g2, b2g2_h[:,i], b2g2_tv[:,i])    


## First, we need to take our differences. We have two that we are taking for each regime:
    ## Land 1 - Gulf 1
    ## Land 2 - Gulf 2
    
land_grad1 = np.subtract(inhtv_ll1, inhtv_lg1)
land_grad2 = np.subtract (inhtv_ll2, inhtv_lg2)

gulf_grad1 = np.subtract(inhtv_gl1, inhtv_gg1)
gulf_grad2 = np.subtract (inhtv_gl2, inhtv_gg2)

b1_grad1 = np.subtract(inhtv_b1l1, inhtv_b1g1)
b1_grad2 = np.subtract (inhtv_b1l2, inhtv_b1g2)

b2_grad1 = np.subtract(inhtv_b2l1, inhtv_b2g1)
b2_grad2 = np.subtract (inhtv_b2l2, inhtv_b2g2)


test = 1000 
boot_dist1 = np.zeros([50,test])
boot_dist2 = np.zeros([50,test])
boot_dist3 = np.zeros([50,test])
boot_dist4 = np.zeros([50,test])
boot_dist5 = np.zeros([50,test])
boot_dist6 = np.zeros([50,test])
boot_dist7 = np.zeros([50,test])
boot_dist8 = np.zeros([50,test])

## Let's try bootstrapping 
for n in range(test):
    sample_inc = np.random.choice(np.arange(0,11), size = 11, replace = True)
    sample = np.mean(land_grad1[:,sample_inc], axis = 1)
    boot_dist1[:,n] = sample   
    
##print(boot_dist1)

upper1 = np.quantile(boot_dist1, 0.925, axis = 1)
lower1 = np.quantile(boot_dist1, 0.075, axis = 1)

##print(upperl1, lowerl1)

## Land 2
for n in range(test):
    sample_inc = np.random.choice(np.arange(0,11), size = 11, replace = True)
    sample = np.mean(land_grad2[:,sample_inc], axis = 1)
    boot_dist2[:,n] = sample   
    

upper2 = np.quantile(boot_dist2, 0.925, axis = 1)
lower2 = np.quantile(boot_dist2, 0.075, axis = 1)


## Gulf 1
for n in range(test):
    sample_inc = np.random.choice(np.arange(0,12), size = 12, replace = True)
    sample = np.mean(gulf_grad1[:,sample_inc], axis = 1)
    boot_dist3[:,n] = sample   
    

upper3 = np.quantile(boot_dist3, 0.925, axis = 1)
lower3 = np.quantile(boot_dist3, 0.075, axis = 1)



## Gulf 2
for n in range(test):
    sample_inc = np.random.choice(np.arange(0,12), size = 12, replace = True)
    sample = np.mean(gulf_grad2[:,sample_inc], axis = 1)
    boot_dist4[:,n] = sample   


upper4 = np.quantile(boot_dist4, 0.925, axis = 1)
lower4 = np.quantile(boot_dist4, 0.075, axis = 1)


## DWL 1
for n in range(test):
    sample_inc = np.random.choice(np.arange(0,12), size = 12, replace = True)
    sample = np.mean(b1_grad1[:,sample_inc], axis = 1)
    boot_dist5[:,n] = sample   
    

upper5 = np.quantile(boot_dist5, 0.925, axis = 1)
lower5 = np.quantile(boot_dist5, 0.075, axis = 1)


##DWL 2
for n in range(test):
    sample_inc = np.random.choice(np.arange(0,12), size = 12, replace = True)
    sample = np.mean(b1_grad2[:,sample_inc], axis = 1)
    boot_dist6[:,n] = sample   


upper6 = np.quantile(boot_dist6, 0.925, axis = 1)
lower6 = np.quantile(boot_dist6, 0.075, axis = 1)



## DWG 1
for n in range(test):
    sample_inc = np.random.choice(np.arange(0,12), size = 12, replace = True)
    sample = np.mean(b2_grad1[:,sample_inc], axis = 1)
    boot_dist7[:,n] = sample   
    

upper7 = np.quantile(boot_dist7, 0.925, axis = 1)
lower7 = np.quantile(boot_dist7, 0.075, axis = 1)


## DWG 2
for n in range(test):
    sample_inc = np.random.choice(np.arange(0,12), size = 12, replace = True)
    sample = np.mean(b2_grad2[:,sample_inc], axis = 1)
    boot_dist8[:,n] = sample   

upper8 = np.quantile(boot_dist8, 0.925, axis = 1)
lower8 = np.quantile(boot_dist8, 0.075, axis = 1)

print(upper8, lower8)




## Now we need to get make our averages from the above variables

lg1_avg = np.nanmean(land_grad1, axis =1)
lg2_avg = np.nanmean(land_grad2, axis =1)
gg1_avg = np.nanmean(gulf_grad1, axis =1)
gg2_avg = np.nanmean(gulf_grad2, axis =1)
b1g1_avg = np.nanmean(b1_grad1, axis =1)
b1g2_avg = np.nanmean(b1_grad2, axis =1)
b2g1_avg = np.nanmean(b2_grad1, axis =1)
b2g2_avg = np.nanmean(b2_grad2, axis =1)
    
# lg1_med = np.nanmedian(land_grad1, axis =1)
# lg2_med = np.nanmedian(land_grad2, axis =1)
# gg1_med = np.nanmedian(gulf_grad1, axis =1)
# gg2_med = np.nanmedian(gulf_grad2, axis =1)
# b1g1_med = np.nanmedian(b1_grad1, axis =1)
# b1g2_med = np.nanmedian(b1_grad2, axis =1)
# b2g1_med = np.nanmedian(b2_grad1, axis =1)
# b2g2_med = np.nanmedian(b2_grad2, axis =1)

lg1_std = np.nanstd(land_grad1, axis =1)
lg2_std = np.nanstd(land_grad2, axis =1)
gg1_std = np.nanstd(gulf_grad1, axis =1)
gg2_std = np.nanstd(gulf_grad2, axis =1)
b1g1_std = np.nanstd(b1_grad1, axis =1)
b1g2_std = np.nanstd(b1_grad2, axis =1)
b2g1_std = np.nanstd(b2_grad1, axis =1)
b2g2_std = np.nanstd(b2_grad2, axis =1)

lg1_p = lg1_avg + lg1_std
lg1_n = lg1_avg - lg1_std
    
lg2_p = lg2_avg + lg2_std
lg2_n = lg2_avg - lg2_std

gg1_p = gg1_avg + gg1_std
gg1_n = gg1_avg - gg1_std
    
gg2_p = gg2_avg + gg2_std
gg2_n = gg2_avg - gg2_std

b1g1_p = b1g1_avg + b1g1_std
b1g1_n = b1g1_avg - b1g1_std
    
b1g2_p = b1g2_avg + b1g2_std
b1g2_n = b1g2_avg - b1g2_std

b2g1_p = b2g1_avg + b2g1_std
b2g1_n = b2g1_avg - b2g1_std
    
b2g2_p = b2g2_avg + b2g2_std
b2g2_n = b2g2_avg - b2g2_std
    
## For plotting ease, we will create a variable of just pressures
##pressure = [1000,975,950,925,900,875,850,825,800,775,750,725,700,675,650,625,600,575,550,525,
       ##     500, 475,450,425,400,375,350,325,300,275,250,225,200,175,150,125,100]



## Now we need to make our figure for theta v

fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots()
ax.plot(lg1_avg, inh_ll1, label = "Land 1", color = "red")
ax.plot(upper1,  inh_ll1, color = "red")
ax.plot(lower1, inh_ll1, color = "red")
plt.fill_betweenx(inh_ll1, upper1, lower1, color = "red",alpha = 0.2)
ax.plot(lg2_avg, inh_ll2, label = "Land 2", color = "red", linestyle='dashed')
ax.plot(upper2, inh_ll2,  color = "red", linestyle='dashed')
ax.plot(lower2, inh_ll2,  color = "red", linestyle='dashed')
plt.fill_betweenx(inh_ll2, lower2, upper2, color = "red",alpha = 0.2)

ax.plot(gg1_avg, inh_ll1, label = "Gulf 1", color = "blue")
ax.plot(upper3, inh_ll1, color = "blue")
ax.plot(lower3, inh_ll1, color = "blue")
plt.fill_betweenx(inh_ll1, upper3, lower3, color = "blue",alpha = 0.2)
ax.plot(gg2_avg, inh_ll1, label = "Gulf 2", color = "blue", linestyle='dashed')
ax.plot(upper4, inh_ll1,  color = "blue", linestyle='dashed')
ax.plot(lower4, inh_ll1,  color = "blue", linestyle='dashed')
plt.fill_betweenx(inh_ll1, lower4, upper4, color = "blue",alpha = 0.2)

##ax.set_title("Theta-v Gradient for Stationary Sites", fontsize = 18)
ax.set_xlabel("Theta-v (K)", fontsize = 14)
ax.set_ylabel("Height (m)", fontsize = 14)
ax.set_ylim(0, 4000)
ax.set_xlim(-3,3)
ax.legend(fontsize = 14)
ax.set_title('(B) LAND vs. GULF', fontsize=16)
ax.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, -0.45), ncol=4)
plt.tight_layout()
plt.savefig('/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 1/gradient_thetav_noshade_lg_85', dpi = 300)



fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots()
ax.plot(lg1_avg, inh_ll1, label = "Land 1", color = "red")
ax.plot(upper1, inh_ll1, color = "red")
ax.plot(lower1, inh_ll1, color = "red")
plt.fill_betweenx(inh_ll1, upper1, lower1, color = "red",alpha = 0.2)
ax.plot(lg2_avg, inh_ll1, label = "Land 2", color = "red", linestyle='dashed')
ax.plot(upper2, inh_ll1,  color = "red", linestyle='dashed')
ax.plot(lower2, inh_ll1,  color = "red", linestyle='dashed')
plt.fill_betweenx(inh_ll1, lower2, upper2, color = "red",alpha = 0.2)


ax.plot(b1g1_avg, inh_ll1, label = "DWL 1", color = "magenta")
ax.plot(upper5, inh_ll1, color = "magenta")
ax.plot(lower5, inh_ll1, color = "magenta")
plt.fill_betweenx(inh_ll1, upper5, lower5, color = "magenta",alpha = 0.2)
ax.plot(b1g2_avg, inh_ll1, label = "DWL 2", color = "magenta", linestyle='dashed')
ax.plot(lower6, inh_ll1,  color = "magenta", linestyle='dashed')
ax.plot(upper6, inh_ll1,  color = "magenta", linestyle='dashed')
plt.fill_betweenx(inh_ll1, upper6, lower6, color = "magenta",alpha = 0.2)

ax.set_title('(C) LAND vs. DWL', fontsize=16)
ax.set_xlabel("Theta-v (K)", fontsize = 14)
ax.set_ylabel("Height (m)", fontsize = 14)
ax.set_ylim(0, 4000)
ax.set_xlim(-3,3)
ax.legend(fontsize = 14)
ax.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, -0.45), ncol=4)
plt.tight_layout()
plt.savefig('/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 1/gradient_thetav_dwl_l_85', dpi = 300)

            
            
fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots()
ax.plot(lg1_avg, inh_ll1, label = "Land 1", color = "red")
ax.plot(upper1, inh_ll1, color = "red")
ax.plot(lower1, inh_ll1, color = "red")
plt.fill_betweenx(inh_ll1, upper1, lower1, color = "red",alpha = 0.2)
ax.plot(lg2_avg, inh_ll1, label = "Land 2", color = "red", linestyle='dashed')
ax.plot(upper2, inh_ll1,  color = "red", linestyle='dashed')
ax.plot(lower2, inh_ll1,  color = "red", linestyle='dashed')
plt.fill_betweenx(inh_ll1, lower2, upper2, color = "red",alpha = 0.2)

ax.plot(b2g1_avg, inh_ll1, label = "DWG 1", color = "purple")
ax.plot(upper7, inh_ll1, color = "purple")
ax.plot(lower7, inh_ll1, color = "purple")
plt.fill_betweenx(inh_ll1, lower7, upper7, color = "purple",alpha = 0.2)
ax.plot(b2g2_avg, inh_ll1, label = "DWG 2", color = "purple", linestyle='dashed')
ax.plot(upper8,  inh_ll1,  color = "purple", linestyle='dashed')
ax.plot(lower8, inh_ll1,  color = "purple", linestyle='dashed')
plt.fill_betweenx(inh_ll1, lower8, upper8, color = "purple",alpha = 0.2)

##ax.set_title("Theta-v Gradient for Stationary Sites", fontsize = 18)
ax.set_xlabel("Theta-v (K)", fontsize = 14)
ax.set_ylabel("Height (m)", fontsize = 14)
ax.set_title('(D) LAND vs. DWG', fontsize=16)
ax.set_ylim(0, 4000)
ax.set_xlim(-3,3)
ax.legend(fontsize = 14)
ax.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, -0.45), ncol=4)
plt.tight_layout()

##plt.savefig('/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 1/gradient_thetav')
plt.savefig('/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 1/gradient_thetav_dwg_l_85', dpi = 300)
# Now we need to make our figure for theta v for median

# fig = plt.figure(figsize=(8,8))
# fig, ax = plt.subplots()
# ax.plot(lg1_med, pressure, label = "Land 1", color = "red")
# ax.plot(lg2_med, pressure, label = "Land 2", color = "red", linestyle='dashed')
# ax.plot(gg1_med, pressure, label = "Gulf 1", color = "blue")
# ax.plot(gg2_med, pressure, label = "Gulf 2", color = "blue", linestyle='dashed')
# ax.plot(b1g1_med, pressure, label = "DWL 1", color = "magenta")
# ax.plot(b1g2_med, pressure, label = "DWL 2", color = "magenta", linestyle='dashed')
# ax.plot(b2g1_med, pressure, label = "DWG 1", color = "purple")
# ax.plot(b2g2_med, pressure, label = "DWG 2", color = "purple", linestyle='dashed')

# ax.set_title("Theta-v Gradient for Stationary Sites", fontsize = 18)
# ax.set_xlabel("Theta-v (K)", fontsize = 14)
# ax.set_ylabel("Pressure (hPa)", fontsize = 14)
# ax.set_ylim(1000, 800)
# ax.set_xlim(-2,3)
# ax.legend(fontsize = 14)
# ax.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, -0.45), ncol=4)
# plt.tight_layout()

# plt.savefig('/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 1/gradient_theta_med', dpi = 300)

# Now we need to make our figure for theta v for median

fig = plt.figure(figsize=(10,10))
fig, ax = plt.subplots()
ax.plot(lg1_avg, inh_ll1, label = "Land 1", color = "red")
ax.plot(lg2_avg, inh_ll1, label = "Land 2", color = "red", linestyle='dashed')
ax.plot(gg1_avg, inh_ll1, label = "Gulf 1", color = "blue")
ax.plot(gg2_avg, inh_ll1, label = "Gulf 2", color = "blue", linestyle='dashed')
ax.plot(b1g1_avg, inh_ll1, label = "DWL 1", color = "magenta")
ax.plot(b1g2_avg, inh_ll1, label = "DWL 2", color = "magenta", linestyle='dashed')
ax.plot(b2g1_avg, inh_ll1, label = "DWG 1", color = "purple")
ax.plot(b2g2_avg, inh_ll1, label = "DWG 2", color = "purple", linestyle='dashed')

##ax.set_title("Theta-v Gradient for Stationary Sites", fontsize = 18)
ax.set_title("(A) Theta-v Gradient for Stationary Sites", fontsize=16)
ax.set_xlabel("Theta-v (K)", fontsize = 14)
ax.set_ylabel("Height (m)", fontsize = 14)
ax.set_ylim(0, 4000)
ax.set_xlim(-2,3)
ax.legend(fontsize = 14)
ax.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, -0.45), ncol=4)
plt.tight_layout()

plt.savefig('/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 1/gradient_theta_mean', dpi = 300)





## this isn't working...which is fine but frustrating
## Now we are going to try and bootstrap our arrays 
# data_1 = (land_grad1)

# bootstrap_ci = bootstrap(data_1[x,:], np.median, confidence_level=0.95,
#                           random_state = 1, method = 'percentile')

# print(bootstrap_ci.confidence_interval)





# ## Let's try and do our theta-e plots now

# ## First we need to take dtheta/dz, which is done using a loop 
# ## We will also do dthetaes/dz


# ## First, we need to take our differences. We have two that we are taking for each regime:
#     ## Land 1 - Gulf 1
#     ## Land 2 - Gulf 2
    
# ## Land  
# ll1_dthetae = np.ones((37,12))*np.nan
# ll2_dthetae = np.ones((37,12))*np.nan
# lg1_dthetae = np.ones((37,12))*np.nan
# lg2_dthetae = np.ones((37,12))*np.nan

# ll1_dthetaes = np.ones((37,12))*np.nan
# ll2_dthetaes = np.ones((37,12))*np.nan
# lg1_dthetaes = np.ones((37,12))*np.nan
# lg2_dthetaes = np.ones((37,12))*np.nan

# ll1_dh = np.ones((37,12))*np.nan
# ll2_dh = np.ones((37,12))*np.nan
# lg1_dh = np.ones((37,12))*np.nan
# lg2_dh = np.ones((37,12))*np.nan

# ## Gulf  
# gl1_dthetae = np.ones((37,12))*np.nan
# gl2_dthetae = np.ones((37,12))*np.nan
# gg1_dthetae = np.ones((37,12))*np.nan
# gg2_dthetae = np.ones((37,12))*np.nan

# gl1_dthetaes = np.ones((37,12))*np.nan
# gl2_dthetaes = np.ones((37,12))*np.nan
# gg1_dthetaes = np.ones((37,12))*np.nan
# gg2_dthetaes = np.ones((37,12))*np.nan

# gl1_dh = np.ones((37,12))*np.nan
# gl2_dh = np.ones((37,12))*np.nan
# gg1_dh = np.ones((37,12))*np.nan
# gg2_dh = np.ones((37,12))*np.nan

# ## DWL 
# b1l1_dthetae = np.ones((37,15))*np.nan
# b1l2_dthetae = np.ones((37,15))*np.nan
# b1g1_dthetae = np.ones((37,15))*np.nan
# b1g2_dthetae = np.ones((37,15))*np.nan

# b1l1_dthetaes = np.ones((37,15))*np.nan
# b1l2_dthetaes = np.ones((37,15))*np.nan
# b1g1_dthetaes = np.ones((37,15))*np.nan
# b1g2_dthetaes = np.ones((37,15))*np.nan

# b1l1_dh = np.ones((37,15))*np.nan
# b1l2_dh = np.ones((37,15))*np.nan
# b1g1_dh = np.ones((37,15))*np.nan
# b1g2_dh = np.ones((37,15))*np.nan

# ## DWG
# b2l1_dthetae = np.ones((37,15))*np.nan
# b2l2_dthetae = np.ones((37,15))*np.nan
# b2g1_dthetae = np.ones((37,15))*np.nan
# b2g2_dthetae = np.ones((37,15))*np.nan

# b2l1_dthetaes = np.ones((37,15))*np.nan
# b2l2_dthetaes = np.ones((37,15))*np.nan
# b2g1_dthetaes = np.ones((37,15))*np.nan
# b2g2_dthetaes = np.ones((37,15))*np.nan

# b2l1_dh = np.ones((37,15))*np.nan
# b2l2_dh = np.ones((37,15))*np.nan
# b2g1_dh = np.ones((37,15))*np.nan
# b2g2_dh = np.ones((37,15))*np.nan

# ## Now we can loop to get our differences that we need for this analysis

# for i in range(0, 12):
#     for j in range(0,37):
#         ## Land Regime
#         ## Land 1
#         ll1_dthetae[j] = ll1_te[j]-ll1_te[j-1]
#         ll1_dthetaes[j] = ll1_tes[j]-ll1_tes[j-1]
#         ll1_dh[j] = ll1_h[j] - ll1_h[j-1]
        
#         ## Land 2
#         ll2_dthetae[j] = ll2_te[j]-ll2_te[j-1]
#         ll2_dthetaes[j] = ll2_tes[j]-ll2_tes[j-1]
#         ll2_dh[j] = ll2_h[j] - ll2_h[j-1]
        
#         ## Gulf 1
#         lg1_dthetae[j] = lg1_te[j]-lg1_te[j-1]
#         lg1_dthetaes[j] = lg1_tes[j]-lg1_tes[j-1]
#         lg1_dh[j] = lg1_h[j] - lg1_h[j-1]
        
#         ## Gulf 2
#         lg2_dthetae[j] = lg2_te[j]-lg2_te[j-1]
#         lg2_dthetaes[j] = lg2_tes[j]-lg2_tes[j-1]
#         lg2_dh[j] = lg2_h[j] - lg2_h[j-1]
        
        
#         ## Gulf Regime
#         ## Land 1
#         gl1_dthetae[j] = gl1_te[j]-gl1_te[j-1]
#         gl1_dthetaes[j] = gl1_tes[j]-gl1_tes[j-1]
#         gl1_dh[j] = gl1_h[j] - gl1_h[j-1]
        
#         ## Land 2
#         gl2_dthetae[j] = gl2_te[j]-gl2_te[j-1]
#         gl2_dthetaes[j] = gl2_tes[j]-gl2_tes[j-1]
#         gl2_dh[j] = gl2_h[j] - gl2_h[j-1]
        
#         ## Gulf 1
#         gg1_dthetae[j] = gg1_te[j]-gg1_te[j-1]
#         gg1_dthetaes[j] = gg1_tes[j]-gg1_tes[j-1]
#         gg1_dh[j] = gg1_h[j] - gg1_h[j-1]
        
#         ## Gulf 2
#         gg2_dthetae[j] = gg2_te[j]-gg2_te[j-1]
#         gg2_dthetaes[j] = gg2_tes[j]-gg2_tes[j-1]
#         gg2_dh[j] = gg2_h[j] - gg2_h[j-1]
        

# ## Need a new loop for our DW domains
# for i in range(0, 15):
#     for j in range(0,37):
#         ## DWL Regime
#         ## Land 1
#         b1l1_dthetae[j] = b1l1_te[j]-b1l1_te[j-1]
#         b1l1_dthetaes[j] = b1l1_tes[j]-b1l1_tes[j-1]
#         b1l1_dh[j] = b1l1_h[j] - b1l1_h[j-1]
        
#         ## Land 2
#         b1l2_dthetae[j] = b1l2_te[j]-b1l2_te[j-1]
#         b1l2_dthetaes[j] = b1l2_tes[j]-b1l2_tes[j-1]
#         b1l2_dh[j] = b1l2_h[j] - b1l2_h[j-1]
        
#         ## Gulf 1
#         b1g1_dthetae[j] = b1g1_te[j]-b1g1_te[j-1]
#         b1g1_dthetaes[j] = b1g1_tes[j]-b1g1_tes[j-1]
#         b1g1_dh[j] = b1g1_h[j] - b1g1_h[j-1]
        
#         ## Gulf 2
#         b1g2_dthetae[j] = b1g2_te[j]-b1g2_te[j-1]
#         b1g2_dthetaes[j] = b1g2_tes[j]-b1g2_tes[j-1]
#         b1g2_dh[j] = b1g2_h[j] - b1g2_h[j-1]
        
        
#         ## DWG Regime
#         ## Land 1
#         b2l1_dthetae[j] = b2l1_te[j]-b2l1_te[j-1]
#         b2l1_dthetaes[j] = b2l1_tes[j]-b2l1_tes[j-1]
#         b2l1_dh[j] = b2l1_h[j] - b2l1_h[j-1]
        
#         ## Land 2
#         b2l2_dthetae[j] = b2l2_te[j]-b2l2_te[j-1]
#         b2l2_dthetaes[j] = b2l2_tes[j]-b2l2_tes[j-1]
#         b2l2_dh[j] = b2l2_h[j] - b2l2_h[j-1]
        
#         ## Gulf 1
#         b2g1_dthetae[j] = b2g1_te[j]-b2g1_te[j-1]
#         b2g1_dthetaes[j] = b2g1_tes[j]-b2g1_tes[j-1]
#         b2g1_dh[j] = b2g1_h[j] - b2g1_h[j-1]
        
#         ## Gulf 2
#         b2g2_dthetae[j] = b2g2_te[j]-b2g2_te[j-1]
#         b2g2_dthetaes[j] = b2g2_tes[j]-b2g2_tes[j-1]
#         b2g2_dh[j] = b2g2_h[j] - b2g2_h[j-1]
 
# ## Now we need to skip the first row because it is dumb

# ## Land
# ll1_dthetae = ll1_dthetae[1:,:]
# ll2_dthetae = ll2_dthetae[1:,:]
# lg1_dthetae = lg1_dthetae[1:,:]
# lg2_dthetae = lg2_dthetae[1:,:]

# ll1_dthetaes = ll1_dthetaes[1:,:]
# ll2_dthetaes = ll2_dthetaes[1:,:]
# lg1_dthetaes = lg1_dthetaes[1:,:]
# lg2_dthetaes = lg2_dthetaes[1:,:]

# ll1_dh = ll1_dh[1:,:]
# ll2_dh = ll2_dh[1:,:]
# lg1_dh = lg1_dh[1:,:]
# lg2_dh = lg2_dh[1:,:]

# ## Gulf
# gl1_dthetae = gl1_dthetae[1:,:]
# gl2_dthetae = gl2_dthetae[1:,:]
# gg1_dthetae = gg1_dthetae[1:,:]
# gg2_dthetae = gg2_dthetae[1:,:]

# gl1_dthetaes = gl1_dthetaes[1:,:]
# gl2_dthetaes = gl2_dthetaes[1:,:]
# gg1_dthetaes = gg1_dthetaes[1:,:]
# gg2_dthetaes = gg2_dthetaes[1:,:]

# gl1_dh = gl1_dh[1:,:]
# gl2_dh = gl2_dh[1:,:]
# gg1_dh = gg1_dh[1:,:]
# gg2_dh = gg2_dh[1:,:]

# ## DWL
# b1l1_dthetae = b1l1_dthetae[1:,:]
# b1l2_dthetae = b1l2_dthetae[1:,:]
# b1g1_dthetae = b1g1_dthetae[1:,:]
# b1g2_dthetae = b1g2_dthetae[1:,:]

# b1l1_dthetaes = b1l1_dthetaes[1:,:]
# b1l2_dthetaes = b1l2_dthetaes[1:,:]
# b1g1_dthetaes = b1g1_dthetaes[1:,:]
# b1g2_dthetaes = b1g2_dthetaes[1:,:]

# b1l1_dh = b1l1_dh[1:,:]
# b1l2_dh = b1l2_dh[1:,:]
# b1g1_dh = b1g1_dh[1:,:]
# b1g2_dh = b1g2_dh[1:,:]

# ##DWG
# b2l1_dthetae = b2l1_dthetae[1:,:]
# b2l2_dthetae = b2l2_dthetae[1:,:]
# b2g1_dthetae = b2g1_dthetae[1:,:]
# b2g2_dthetae = b2g2_dthetae[1:,:]

# b2l1_dthetaes = b2l1_dthetaes[1:,:]
# b2l2_dthetaes = b2l2_dthetaes[1:,:]
# b2g1_dthetaes = b2g1_dthetaes[1:,:]
# b2g2_dthetaes = b2g2_dthetaes[1:,:]

# b2l1_dh = b2l1_dh[1:,:]
# b2l2_dh = b2l2_dh[1:,:]
# b2g1_dh = b2g1_dh[1:,:]
# b2g2_dh = b2g2_dh[1:,:]



# ## Now we need to calculate what we need
        
# ## Here is the division we need 

# ## Land
# ll1_dtdz = np.ones((36,12))*np.nan
# ll2_dtdz = np.ones((36,12))*np.nan
# lg1_dtdz = np.ones((36,12))*np.nan
# lg2_dtdz = np.ones((36,12))*np.nan

# ll1_dtdz = np.divide(ll1_dthetae, ll1_dh)
# ll2_dtdz = np.divide(ll2_dthetae, ll2_dh)
# lg1_dtdz = np.divide(lg1_dthetae, lg1_dh)
# lg2_dtdz = np.divide(lg2_dthetae, lg2_dh)

# ll1_dtsdz = np.ones((36,12))*np.nan
# ll2_dtsdz = np.ones((36,12))*np.nan
# lg1_dtsdz = np.ones((36,12))*np.nan
# lg2_dtsdz = np.ones((36,12))*np.nan

# ll1_dtsdz = np.divide(ll1_dthetaes, ll1_dh)
# ll2_dtsdz = np.divide(ll2_dthetaes, ll2_dh)
# lg1_dtsdz = np.divide(lg1_dthetaes, lg1_dh)
# lg2_dtsdz = np.divide(lg2_dthetaes, lg2_dh)


# ## Gulf
# gl1_dtdz = np.ones((36,12))*np.nan
# gl2_dtdz = np.ones((36,12))*np.nan
# gg1_dtdz = np.ones((36,12))*np.nan
# gg2_dtdz = np.ones((36,12))*np.nan

# gl1_dtdz = np.divide(gl1_dthetae, gl1_dh)
# gl2_dtdz = np.divide(gl2_dthetae, gl2_dh)
# gg1_dtdz = np.divide(gg1_dthetae, gg1_dh)
# gg2_dtdz = np.divide(gg2_dthetae, gg2_dh)

# gl1_dtsdz = np.ones((36,12))*np.nan
# gl2_dtsdz = np.ones((36,12))*np.nan
# gg1_dtsdz = np.ones((36,12))*np.nan
# gg2_dtsdz = np.ones((36,12))*np.nan

# gl1_dtsdz = np.divide(gl1_dthetaes, gl1_dh)
# gl2_dtsdz = np.divide(gl2_dthetaes, gl2_dh)
# gg1_dtsdz = np.divide(gg1_dthetaes, gg1_dh)
# gg2_dtsdz = np.divide(gg2_dthetaes, gg2_dh)


# ## DWL
# b1l1_dtdz = np.ones((36,15))*np.nan
# b1l2_dtdz = np.ones((36,15))*np.nan
# b1g1_dtdz = np.ones((36,15))*np.nan
# b1g2_dtdz = np.ones((36,15))*np.nan

# b1l1_dtdz = np.divide(b1l1_dthetae, b1l1_dh)
# b1l2_dtdz = np.divide(b1l2_dthetae, b1l2_dh)
# b1g1_dtdz = np.divide(b1g1_dthetae, b1g1_dh)
# b1g2_dtdz = np.divide(b1g2_dthetae, b1g2_dh)

# b1l1_dtsdz = np.ones((36,15))*np.nan
# b1l2_dtsdz = np.ones((36,15))*np.nan
# b1g1_dtsdz = np.ones((36,15))*np.nan
# b1g2_dtsdz = np.ones((36,15))*np.nan

# b1l1_dtsdz = np.divide(b1l1_dthetaes, b1l1_dh)
# b1l2_dtsdz = np.divide(b1l2_dthetaes, b1l2_dh)
# b1g1_dtsdz = np.divide(b1g1_dthetaes, b1g1_dh)
# b1g2_dtsdz = np.divide(b1g2_dthetaes, b1g2_dh)


# ## DWG
# b2l1_dtdz = np.ones((36,15))*np.nan
# b2l2_dtdz = np.ones((36,15))*np.nan
# b2g1_dtdz = np.ones((36,15))*np.nan
# b2g2_dtdz = np.ones((36,15))*np.nan

# b2l1_dtdz = np.divide(b2l1_dthetae, b2l1_dh)
# b2l2_dtdz = np.divide(b2l2_dthetae, b2l2_dh)
# b2g1_dtdz = np.divide(b2g1_dthetae, b2g1_dh)
# b2g2_dtdz = np.divide(b2g2_dthetae, b2g2_dh)

# b2l1_dtsdz = np.ones((36,15))*np.nan
# b2l2_dtsdz = np.ones((36,15))*np.nan
# b2g1_dtsdz = np.ones((36,15))*np.nan
# b2g2_dtsdz = np.ones((36,15))*np.nan

# b2l1_dtsdz = np.divide(b2l1_dthetaes, b2l1_dh)
# b2l2_dtsdz = np.divide(b2l2_dthetaes, b2l2_dh)
# b2g1_dtsdz = np.divide(b2g1_dthetaes, b2g1_dh)
# b2g2_dtsdz = np.divide(b2g2_dthetaes, b2g2_dh)

# ## Here is the gradient

# l1_grad_te = np.subtract(ll1_dtdz, lg1_dtdz)
# l2_grad_te = np.subtract(ll2_dtdz, lg2_dtdz)

# g1_grad_te = np.subtract(gl1_dtdz, gg1_dtdz)
# g2_grad_te = np.subtract(gl2_dtdz, gg2_dtdz)

# b11_grad_te = np.subtract(b1l1_dtdz, b1g1_dtdz)
# b12_grad_te = np.subtract(b1l2_dtdz, b1g2_dtdz)

# b21_grad_te = np.subtract(b2l1_dtdz, b2g1_dtdz)
# b22_grad_te = np.subtract(b2l2_dtdz, b2g2_dtdz)

# l1_grad_tes = np.subtract(ll1_dtsdz, lg1_dtsdz)
# l2_grad_tes = np.subtract(ll2_dtsdz, lg2_dtsdz)

# g1_grad_tes = np.subtract(gl1_dtsdz, gg1_dtsdz)
# g2_grad_tes = np.subtract(gl2_dtsdz, gg2_dtsdz)

# b11_grad_tes = np.subtract(b1l1_dtsdz, b1g1_dtsdz)
# b12_grad_tes = np.subtract(b1l2_dtsdz, b1g2_dtsdz)

# b21_grad_tes = np.subtract(b2l1_dtsdz, b2g1_dtsdz)
# b22_grad_tes = np.subtract(b2l2_dtsdz, b2g2_dtsdz)


# ## Here are our means 
# te_lg1_avg = np.nanmean(l1_grad_te, axis =1)
# te_lg2_avg = np.nanmean(l2_grad_te, axis =1)

# te_gg1_avg = np.nanmean(g1_grad_te, axis =1)
# te_gg2_avg = np.nanmean(g2_grad_te, axis =1)

# te_b1g1_avg = np.nanmean(b11_grad_te, axis =1)
# te_b1g2_avg = np.nanmean(b12_grad_te, axis =1)

# te_b2g1_avg = np.nanmean(b21_grad_te, axis =1)
# te_b2g2_avg = np.nanmean(b22_grad_te, axis =1)

# tes_lg1_avg = np.nanmean(l1_grad_tes, axis =1)
# tes_lg2_avg = np.nanmean(l2_grad_tes, axis =1)

# tes_gg1_avg = np.nanmean(g1_grad_tes, axis =1)
# tes_gg2_avg = np.nanmean(g2_grad_tes, axis =1)

# tes_b1g1_avg = np.nanmean(b11_grad_tes, axis =1)
# tes_b1g2_avg = np.nanmean(b12_grad_tes, axis =1)

# tes_b2g1_avg = np.nanmean(b21_grad_tes, axis =1)
# tes_b2g2_avg = np.nanmean(b22_grad_tes, axis =1)

# ## Skip pressure = 1000 because we don't need it
    
# p2 = [975,950,925,900,875,850,825,800,775,750,725,700,675,650,625,600,575,550,525,
#             500, 475,450,425,400,375,350,325,300,275,250,225,200,175,150,125,100]

# ## Now we need to plot this: theta-e

# fig = plt.figure(figsize=(8,8))
# fig, ax = plt.subplots()
# ax.plot(te_lg1_avg, p2, label = "Land 1", color = "red")
# ax.plot(te_lg2_avg, p2, label = "Land 2", color = "red", linestyle='dashed')
# ax.plot(te_gg1_avg, p2, label = "Gulf 1", color = "blue")
# ax.plot(te_gg2_avg, p2, label = "Gulf 2", color = "blue", linestyle='dashed')
# ax.plot(te_b1g1_avg, p2, label = "DWL 1", color = "magenta")
# ax.plot(te_b1g2_avg, p2, label = "DWL 2", color = "magenta", linestyle='dashed')
# ax.plot(te_b2g1_avg, p2, label = "DWG 1", color = "purple")
# ax.plot(te_b2g2_avg, p2, label = "DWG 2", color = "purple", linestyle='dashed')


# ax.set_title("d theta-e / dz of Horizontal Gradient", fontsize = 18)
# ax.set_xlabel("d theta-e / dz (K/m)", fontsize = 14)
# ax.set_ylabel("Pressure (hPa)", fontsize = 14)
# ax.set_ylim(1000, 800)
# ax.set_xlim(-0.015,0.02)
# ax.legend(fontsize = 14)
# ax.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, -0.45), ncol=4)
# plt.tight_layout()

# ##plt.savefig('/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 1/gradient_thetav')
# plt.savefig('/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 1/gradient_thetae')



# ## Now we need to plot this: theta-es

# fig = plt.figure(figsize=(8,8))
# fig, ax = plt.subplots()
# ax.plot(tes_lg1_avg, p2, label = "Land 1", color = "red")
# ax.plot(tes_lg2_avg, p2, label = "Land 2", color = "red", linestyle='dashed')
# ax.plot(tes_gg1_avg, p2, label = "Gulf 1", color = "blue")
# ax.plot(tes_gg2_avg, p2, label = "Gulf 2", color = "blue", linestyle='dashed')
# ax.plot(tes_b1g1_avg, p2, label = "DWL 1", color = "magenta")
# ax.plot(tes_b1g2_avg, p2, label = "DWL 2", color = "magenta", linestyle='dashed')
# ax.plot(tes_b2g1_avg, p2, label = "DWG 1", color = "purple")
# ax.plot(tes_b2g2_avg, p2, label = "DWG 2", color = "purple", linestyle='dashed')


# ax.set_title("d theta-es / dz of Horizontal Gradient", fontsize = 18)
# ax.set_xlabel("d theta-es / dz (K/m)", fontsize = 14)
# ax.set_ylabel("Pressure (hPa)", fontsize = 14)
# ax.set_ylim(1000, 800)
# ax.set_xlim(-0.015,0.02)
# ax.legend(fontsize = 14)
# ax.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, -0.45), ncol=4)
# plt.tight_layout()

# ##plt.savefig('/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 1/gradient_thetav')
# plt.savefig('/Users/kem6245/Documents/Python Copy/ESCAPE/ESCAPE_Figures/Paper 1/gradient_thetaes')





