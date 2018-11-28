#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:19:23 2017 updated by HLH August 2018S

@author: tojo5760
"""
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib

infolder = '/home/haho3703/TIEGCM/Contour_textfiles/'
LT_ticks = np.array([12., 15., 18., 21., 0., 3., 6., 9., 12.])
data = np.zeros((120,240))
data = np.loadtxt(infolder + "He_dens_RealVsDiff_pdrag.txt", delimiter=',')
# data = np.loadtxt(infolder + "He_dens_Calculated_Diff_pdrag.txt", delimiter=',')
# data = np.loadtxt(infolder + "He_dens_model_400km_pdrag.txt", delimiter=',')
# data = np.loadtxt(infolder + "He_dens_model_120km_pdrag.txt", delimiter=',')

#Load magnetic equator data points to be overlayed on plot
mag_equator = np.loadtxt("Magnetic_equator_lat_lon.txt",delimiter=',')
marklon = 0
magnetic_x = np.linspace(0,143/6,360)

x = np.linspace(0,143/6, 144)
y = np.linspace(-88.75,88.75,72)
X, Y = np.meshgrid(x, y)

plt.figure()
levels = np.linspace(np.amin(data),np.amax(data),100)
cont = plt.contour(X, Y, data, 10, colors='k')

# ----- uncomment for percent differences -----------------------
# levels = np.linspace(-65, 65, 300)        # set your own limits for colorbar
# myplot = plt.contourf(X, Y, data, levels, cmap='seismic')
# cbar = plt.colorbar(myplot, format='%.0f')
# cbar.ax.set_ylabel('Percent Difference')
# plt.title('He Deviation from Diffusive Eq. [z = 400km, UT=0, pdrag]')


# ---- uncomment for Helium Density ------------------------------
myplot = plt.contourf(X, Y, data,levels,cmap='jet')
cbar = plt.colorbar(myplot, format='%.2e')
cbar.ax.set_ylabel('Number Density [#/cm^3]')
# plt.title('Calculated Diffusive Eq. He Density [z = 400km, UT=0, pdrag]')
plt.title('TIEGCM He Density [z = 120km, UT=0, pdrag]')

plt.xticks(np.arange(0.,24.,3.), LT_ticks)
plt.plot(magnetic_x,mag_equator[:,1],'r') #Add magnetic equator line
plt.xlabel('Local Solar Time [hr]')
plt.ylabel('Latitude [deg]')
plt.show()
ctrmin = np.amin(data)
ctrmax = np.amax(data)
