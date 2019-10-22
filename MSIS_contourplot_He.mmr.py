# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 12:07:31 2017

@author: tojo5760
"""
#==============================================================================
# Plot the contour plot of MSIS outputs
#==============================================================================

from pyglow.pyglow import Point
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib
import sys

lat = 0.
lon = 0.
dn = datetime(1992, 3, 20, 0, 2)#year,month,day,hour,minute
Av = 6.022141*10**23
whatplot = 'Temp'            # either He or Temp
LT_ticks = np.array([12., 15., 18., 21., 0., 3., 6., 9., 12.])
#==============================================================================
# pt = Point(dn, lat, lon, 400)
# result = pt.run_msis()
# print pt
# print result.f107
#==============================================================================


#Load magnetic equator data points to be overlayed on plot
mag_equator = np.loadtxt("Magnetic_equator_lat_lon.txt",delimiter=',')
top = mag_equator[180:360,:]
bottom = mag_equator[0:180,:]
mag_equator = np.concatenate((top,bottom),axis=0)

data = np.zeros((120,240))
marklon = 0

#Loop through lat and lon. select what MSIS values are desired
for Lon in range(-120,120,1) :
    marklat = 0
    for Lat in range(-60,60,1):
        pt = Point(dn, Lat*1.5, Lon*1.5, 400)#1.5 degree grid size
        pt.f107a = 180
        pt.f107p = 180  # previous day's F10.7
        pt.apmsis = [0,]*7
        result = pt.run_msis()


        #data[marklat,marklon] = math.log10((result.nn['HE']*4/Av)/(result.rho-result.nn['HE']*4/Av),10) #Derived helium
        if whatplot == 'He':
            data[marklat,marklon] = result.nn['HE']*.004003/Av*10**6 # MSIS helium (Should be the same as derived)
        elif whatplot == 'Temp':
            data[marklat,marklon] = result.Tn_msis #Neutral Temperature
        else:
            print('Bad Plot input')
            sys.exit()
        marklat = marklat+1
    marklon = marklon+1

print(result.f107)
print(result.f107a)
print(result.apmsis)
#np.savetxt('mmrdata',data)


print(len(data[0,:]))

#Create grid to plot data onto
x = np.linspace(0,143/6, 240)
#x = np.linspace(-180, 180, 240)
y = np.linspace(-88.75,88.75,120)
X, Y = np.meshgrid(x, y)
magnetic_x = np.linspace(0,143/6,360)

plt.figure()
#levels = np.linspace(np.amin(data),np.amax(data),100)
levels = np.linspace(890,1325,100)
myplot = plt.contourf(X, Y, data,levels,cmap='jet')
cont = plt.contour(X,Y,data,10,colors='k')
cbar = plt.colorbar(myplot)
cbar.ax.set_ylabel('Neutral Temp [K]')
plt.xticks(np.arange(0.,25, 3.), LT_ticks)
plt.plot(magnetic_x,mag_equator[:,1],'r') #Add magnetic equator line
plt.title('MSIS Neutral Temp: $Z_{gm}$=400km, UT=0')
plt.xlabel('Local Solar Time [hr]')
plt.ylabel('Latitude [deg]')
plt.show()
msismin = np.amin(data)
msismax = np.amax(data)
