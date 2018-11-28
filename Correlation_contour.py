#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 13:41:14 2018
Plots the correlation matrix created by the matlab code in Contour_Plots.m
@author: tojo5760
"""

import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib

#Load data from text file
data = np.loadtxt("He_Tn_Correlation_400km_pdrag.txt",delimiter=',')

#Create contour plot mesh of correct size
x = np.linspace(-143,143,287)
y = np.linspace(-71,71,143)
X, Y = np.meshgrid(x,y)

#Plot contour
plt.figure()
levels = np.linspace(np.amin(data),np.amax(data),100)
myplot = plt.contourf(X, Y, data,levels,cmap='jet')
cont = plt.contour(X,Y,data,10,colors='k')
cbar = plt.colorbar(myplot, format='%.2e')
cbar.ax.set_ylabel('Correleation Coefficient')
plt.xticks(np.arange(-143.,143.,26.))
plt.title('Normalized He Tn Correlation plot 400km with Ion Drag')
plt.xlabel('Longitude shift')
plt.ylabel('Latitude shift')
plt.show()
