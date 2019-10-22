#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 12:06:48 2018
Plots the ETA of MSIS outputs over given longitude from latitude of -60 to 60
@author: tojo5760
"""

from pyglow.pyglow import Point
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib

lon = 30.
dn = datetime(1992, 3, 20, 0, 2)#year,month,day,hour,minute
Av = 6.022141*10**23
#==============================================================================
# pt = Point(dn, lat, lon, 400)
# result = pt.run_msis()
# print pt
# print result.f107
#==============================================================================

data = np.zeros(121)
marklon = 0
latitudes = range(-60,61,1)

#Run MSIS for range of latitudes
for Lat in range(-60,61,1) :
    pt = Point(dn, Lat, lon, 400)
    result = pt.run_msis()
    data[marklon] = result.Tn_msis
    marklon = marklon+1

f=result.f107
fa=result.f107a
apmsis=result.apmsis

#Plot the data
fig = plt.figure()
plt.plot(latitudes,data)
plt.title('Temperature ETA MSIS lon=%.1f UT=0.03 F10.7=%.1f'%(lon,f),fontweight='bold')
plt.ylabel('Neutral Temperature (K)')
plt.xlabel('Latitude')
