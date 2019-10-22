# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 09:01:08 2017

@author: Torfinn Johnsrud
"""

#==============================================================================
# Plot Knudsen number for a spacecraft with characeristic length of 1m and for
# the atmosphere scale height. Also plot mass percentage of Helium as a 
# function of altitude. All from MSIS model. Definition of msis 
# is in pyglow.py line 290.
# For all plots, "low" and "high" refer to the solar cycle.
#==============================================================================
from pyglow.pyglow import Point
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib

dn = datetime(1987, 3, 23, 9, 30)#year,month,day,hour,minute
dn2 = datetime(1981, 3, 23, 9, 30)
lat = 0.
lon = -80.
m = range(100, 1000, 5) #Set altitude range and step size (km)
Dens1frac = np.zeros(len(m))
Dens2frac = np.zeros(len(m))
H1 = np.zeros(len(m))#Pressure scale heights
H2 = np.zeros(len(m))
M1 = np.zeros(len(m))#Mean molecular mass
T1 = np.zeros(len(m))#Neutral temperature
T2 = np.zeros(len(m))
M2 = np.zeros(len(m))
Kn1 = np.zeros(len(m))#Knudsen number for characteristic length
Kn2 = np.zeros(len(m))
Kn1At = np.zeros(len(m))#Atmospheric knudsen number
Kn2At = np.zeros(len(m))
k = 1.3806*10**(-23) #Boltzmann's constant
g0 = 9.81 #m/s^2
R = 6373 #Radius of earth in km
Av = 6.022141*10**23 #Avogadro's Constant
n = 0
rho = np.zeros(len(m))


#Iterate through altitudes in range m
for x in m:
    pt = Point(dn, lat, lon, x)
    pt2 = Point(dn2, lat, lon, x)
    result = pt.run_msis()
    result2 = pt2.run_msis()
    
    #Knudsen number estimate
    n_dens_tot = result.nn['HE']+result.nn['N2']+result.nn['O2'] \
        +result.nn['AR']+result.nn['H']+result.nn['O']+result.nn['N']# number density per cm^3
    
    d_avg = (result.nn['HE']/n_dens_tot*260+result.nn['N2']/n_dens_tot*364\
        +result.nn['O2']/n_dens_tot*346+result.nn['AR']/n_dens_tot*340)\
        *10**(-12)+(result.nn['H']+result.nn['O']+result.nn['N'])/n_dens_tot*5.046*10**(-10)# average kinetic diameter of 
        # molecules in atmosphere. Note H and O are ignored due to no data on
        # their kinetic diameters and the breakdown of the model with ions 
        # (hard sphere collisions of molecules)   
    lam = 1/(math.sqrt(2)*math.pi*d_avg**2*n_dens_tot*10**(6))
        # lam = 1/(pi*Nv*D^2*sqrt(2))
    L = 1 #in meters. Chosen characteristic length
    Kn1[n] = lam/L
    
    g = g0*R**2/(R+x)**2
    m_dens1 = (result.nn['HE']*.004003/Av+result.nn['N2']*.0280134/Av\
    +result.nn['O2']*.032/Av+result.nn['AR']*.039948/Av+result.nn['N']\
    *.01401/Av+result.nn['H']*.001/Av+result.nn['O']*.016/Av)#kg/cm^3
    
    M1[n] = m_dens1/(n_dens_tot)
    T1[n] = result.Tn_msis
    H1[n] = k*T1[n]/(M1[n]*g)
    Kn1At[n] = lam/(H1[n])
    
    n_dens_tot2 = result2.nn['HE']+result2.nn['N2']+result2.nn['O2'] \
        +result2.nn['AR']+result2.nn['H']+result2.nn['O']+result2.nn['N']# number density per cm^3
    
    d_avg2 = (result2.nn['HE']/n_dens_tot2*260+result2.nn['N2']/n_dens_tot2*364\
        +result2.nn['O2']/n_dens_tot2*346+result2.nn['AR']/n_dens_tot2*340)\
        *10**(-12)+(result2.nn['H']+result2.nn['O']+result2.nn['N'])/n_dens_tot2*5.046*10**(-10)# average kinetic diameter of 
        # molecules in atmosphere. Note H and O are ignored due to no data on
        # their kinetic diameters and the breakdown of the model with ions 
        # (hard sphere collisions of molecules)
    lam2 = 1/(math.sqrt(2)*math.pi*d_avg2**2*n_dens_tot2*10**(6))
        # lam = 1/(pi*Nv*D^2*sqrt(2)) in meters
    Kn2[n] = lam2/L    
   
    m_dens2 = (result2.nn['HE']*.004003/Av+result2.nn['N2']*.0280134/Av\
    +result2.nn['O2']*.032/Av+result2.nn['AR']*.039948/Av+result2.nn['N']\
    *.01401/Av+result2.nn['H']*.001/Av+result2.nn['O']*.016/Av)#kg/cm^3
    
    M2[n] = m_dens2/(n_dens_tot2)
    T2[n] = result2.Tn_msis
    H2[n] = k*T2[n]/(M2[n]*g) # in meters
    Kn2At[n] = lam2/(H2[n])
    
    #Helium Density plots
    He1 = result.nn['HE']*.004003/Av 
    He2 = result2.nn['HE']*.004003/Av
    Dens1frac[n] = He1/(result.rho/1000)*100
    Dens2frac[n] = He2/(result2.rho/1000)*100
    
    n = n+1

F107_1 = result.f107  
F107_2 = result2.f107

#Plot Knudsen number as a function of altitude with charictaristic length 
#of 1 meter
#==============================================================================
fig = plt.figure()
ax = plt.subplot(111)
line1, = ax.plot(Kn1, m, label='1987 F10.7 = %.1f' %F107_1, color='red')
line2, = ax.plot(Kn2, m, label='1981 F10.7 = %.1f' %F107_2, color='blue')
line3, = ax.plot(Kn1At, m, label='Atmospheric, F10.7 = %.1f' %F107_1, color='red', ls='dashed')
line4, = ax.plot(Kn2At, m, label='Atmospheric, F10.7 = %.1f' %F107_2, color='blue', ls='dashed')

plt.xscale('log')
plt.xlabel('Knudsen Number')
plt.ylabel('Altitude (km)')
plt.title('Knudsen Number versus Altitude')
#ax.legend(loc='upper center', bbox_to_anchor=(.5, .2 ), ncol=2)

red_patch = matplotlib.patches.Patch(color='red', label='F10.7 = %.1f' %F107_1)
blue_patch = matplotlib.patches.Patch(color='blue', label='F10.7 = %.1f' %F107_2)
black_line = matplotlib.lines.Line2D([], [], color='black', label='L = 1m')
black_dash = matplotlib.lines.Line2D([], [], ls='dashed', color='black', label='L=Hp scale height')
plt.legend(handles=[red_patch, blue_patch, black_line, black_dash], loc='upper left', bbox_to_anchor=(0, 1))
plt.grid()
plt.show()

#Helium Mass percentage plot
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(Dens1frac, m, label='F10.7 = %.1f' %F107_1)
ax.plot(Dens2frac, m, label='F10.7 = %.1f' %F107_2)
plt.xscale('log')
plt.xlabel('Percentage')
plt.ylabel('Altitude (km)')
plt.title('Helium Mass Percentage vs Altitude')
plt.grid()
z = np.zeros(len(m))+50
ax.plot(z,m, label='50% Threshold',ls='dashed')
ax.legend(loc='upper left', bbox_to_anchor=(0, 1), ncol=1)
plt.show()

#Pressure scale height plot
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(H1/1000, m, label='Low')
ax.plot(H2/1000, m, label='High')
plt.title('Pressure Scale Height')
plt.ylabel('Altitude (km)')
plt.xlabel('Scale Height (km)')
ax.legend()
plt.show()

#Mean molecular mass with altitude plot
fig =plt.figure()
ax = plt.subplot(111)
ax.plot(M1*Av*1000, m, label='Low')
ax.plot(M2*Av*1000, m, label='High')
plt.xlabel('Mean Molecular Mass (kg/kmol)')
plt.ylabel('Altitude (km)')
plt.title('Avg. Molecular Mass vs Altitude')
ax.legend()
plt.show()

#Neutral Temperature with altitude plot
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(T1, m, label='Low')
ax.plot(T2, m, label='High')
plt.xlabel('Temperature (K)')
plt.ylabel('Altitude (km)')
plt.title('Temperature vs Altitude')
ax.legend()
plt.show()

#Temperature/Mass plot
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(m,T1/M1, label='Low')
ax.plot(m,T2/M2, label='High')
plt.xlabel('Altitude (km)')
plt.ylabel('T/M')
plt.title('Temperature/Molecular mass')
ax.legend()
plt.show()


