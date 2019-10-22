#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:19:23 2017
Plots Contour plots created by Contour_Plots.m
@author: tojo5760
Last Updated by : Hannah Holt, 6/6/18
"""
import matplotlib.pyplot as plt
import numpy as np
import math
import datetime
from matplotlib.dates import DateFormatter
import sys

infolder = '/home/haho3703/TIEGCM/Contour_textfiles/'
Zp_name = '3,125'
z_geometric = '120'
z_geopotential = '400'
pdrag = 1                           # 0 for no ion drag, 1 for ion drag
which_vertical = 'geometric alt' # to specify the vertical coordinate to be used. Can be geometric alt, geopotential alt, or log pressure.
name = infolder + 'Neutral_Temp_'                # can be 'He_Density_' or 'Neutral_Temp_'

LT_ticks = np.array([12., 15., 18., 21., 0., 3., 6., 9., 12.])
Zp_real = Zp_name.replace(',', '.')
if pdrag == 0:
    filetype = 'ctrSS'
    description = 'no Ion Drag'
else:
    filetype = 'pdrag'
    description = 'pdrag'

#Load magnetic equator data points to be overlayed on plot
mag_equator = np.loadtxt("Magnetic_equator_lat_lon.txt",delimiter=',')
data = np.zeros((120,240))
marklon = 0

#Load data file
if which_vertical == 'log pressure':
    data = np.loadtxt(name + "Zp_" + filetype + '_' + Zp_name + ".txt", delimiter=',')
    label = '$Z_{p} = $' + Zp_real  # for graph labeling
elif which_vertical == 'geometric alt':
    data = np.loadtxt(name + filetype + '_' + z_geometric + 'km(geometric)_' + ".txt", delimiter=',')
    label = 'z = ' + z_geometric + ' km'
elif which_vertical == 'geopotential alt':
    data = np.loadtxt(name + filetype + '_' + z_geopotential + 'km(geopotential)_' + ".txt", delimiter=',')
    label = '$Z_{gp} = $' + z_geopotential + ' km'
else:
    print('ERROR: Bad vertical coordinate.')
    sys.exit()

#Create contour plot mesh of correct size
x = np.linspace(0,143/6, 144)               # longitudes
y = np.linspace(-88.75,88.75,72)            # latitudes
X, Y = np.meshgrid(x, y)
magnetic_x = np.linspace(0,143/6,360)

#Plot contour
plt.figure()
levels = np.linspace(np.amin(data),np.amax(data),100)
# levels = np.linspace(1.50E-14,7E-14,100)
myplot = plt.contourf(X, Y, data, levels, cmap='jet')
cont = plt.contour(X,Y,data, 10, colors='k')

if name == infolder + 'He_Density_':
    cbar = plt.colorbar(myplot, format='%.2e')
    cbar.ax.set_ylabel('He Density [kg/m^3]')
    plt.title('TIEGCM He Density: ' + label + ', UT=0, ' + description + ']')

elif name == infolder + 'Neutral_Temp_':
    cbar = plt.colorbar(myplot)
    cbar.ax.set_ylabel('Neutral Temp [K]')
    plt.title('TIEGCM Neutral Temp [' + label + ', UT=0, ' + description + ']')

plt.xticks(np.arange(0.,25, 3.), LT_ticks)
#plt.ylim(-41,41)
plt.plot(magnetic_x,mag_equator[:,1],'r') #Add magnetic equator line
plt.xlabel('Local Solar Time [hr]')
plt.ylabel('Latitude [deg]')
plt.show()
ctrmin = np.amin(data)
ctrmax = np.amax(data)


# extra stuff to look at the geopotential altitude for a specific pressure level
if which_vertical == 'log pressure':

    data = np.loadtxt("GeopotentialAlt_Zp_" + filetype + "_" + Zp_name + ".txt", delimiter=',')

    #Plot geometric altitude for the specfic Zp lvl
    plt.figure()
    #levels = np.linspace(1.5153e-14,6.8968999999999999e-14,100)
    levels = np.linspace(np.amin(data),np.amax(data),100)
    myplot = plt.contourf(X, Y, data,levels,cmap='jet')
    cont = plt.contour(X,Y,data,10,colors='k')
    cbar = plt.colorbar(myplot)
    cbar.ax.set_ylabel('Geopotential Altitude [km]')
    plt.xticks(np.arange(0.,25, 3.), LT_ticks)
    plt.title('Geopotential Altitude: Zp=' + Zp_real + ', UT=0 ' + description)
    plt.xlabel('Local Solar Time [hr]')
    plt.ylabel('Latitude [deg]')
    plt.show()
    ctrmin = np.amin(data)
    ctrmax = np.amax(data)