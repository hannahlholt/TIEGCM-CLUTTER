# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 12:56:54 2017

@author: Torfinn Johnsrud
"""
#==============================================================================
#This plots the solar cycle variations of He and O density at 400km
#==============================================================================

from pyglow.pyglow import Point
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

month = np.array([1,4,7,10])#specify what months to pull data from
year = range(1981, 2002, 1)
He = np.zeros(len(month)*len(year))
Oxy = np.zeros(len(He))
F107 = np.zeros(len(He))
lat = 40.
lon = -80.
n = 0

for x in year:
    for y in month:
        dn = datetime(x,y,1,12,0)
        pt = Point(dn, lat, lon, 400)
        result = pt.run_msis()
        He[n] = result.nn['HE']
        Oxy[n] = result.nn['O']
        F107[n] = result.f107
        n = n+1
        
He_avg = np.sum(He)/len(He)#Find mean value
He_delta = He/He_avg*100#Find percent change
Oxy_avg = np.sum(Oxy)/len(Oxy)
Oxy_delta = Oxy/Oxy_avg*100
cycle = np.linspace(year[0], year[-1], len(He))

fig = plt.figure()
ax = plt.subplot(111)
ax.plot(cycle, He_delta, label='Helium')
ax.plot(cycle, Oxy_delta, label='Oxygen')
ax.plot(cycle, F107, label='F10.7')
plt.title('Concentration Variation with Solar Cycle 400km Altitude')
plt.ylabel('Percent Deviation from Mean')
plt.xlabel('Year')
ax.legend()
plt.show()

fig = plt.figure()
ax = plt.subplot(111)
ax.plot(cycle, He, label='Helium')
ax.plot(cycle, Oxy, label='Oxygen')
plt.title('Concentration Variation with Solar Cycle 400km Altitude')
plt.ylabel('Concentration (cm^-3)')
plt.xlabel('Year')
ax.legend()
plt.show()
