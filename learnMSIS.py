'''
**** MSIS NEEDS PYTHON 2.7 TO COMPILE AND RUN ****

learnMSIS.py is a program designed to help the user understand basic features and uses of the MSIS model.

Before you can run this program, you must download the pyglow module. This can be found at:
https://github.com/timduly4/pyglow
Select the green "clone or download" button. Then run the 'pyglow_install.sh' script from the command line. You can also
check out the README file for more help. (Remember Python 2.7!!!)
You should be able to run this program and any other pyglow features after that.

'''

from __future__ import print_function, division     # python 3 has better features for printing and math, so import those
import os
from pyglow.pyglow import Point                     # imports the pyglow Point module and all its features
from datetime import datetime                       # MSIS needs a datetime object
import numpy as np                                  # numpy is your best friend
import matplotlib.pyplot as plt                     # for plotting.

def print_header():
    now = datetime.now()
    print('\n**********************************************')
    print(os.path.basename(__file__))
    print('Written by: Hannah Holt')
    print('Last updated: 6/21/2018')
    print('**********************************************\n')


print_header()

'''
Every "Point" object needs four things to build itself

1)  Param:dn    datetime object - i.e. a date that you want to look at
2)  Param:lat   latitude for the Point (in degrees from -90 to 90)
3)  Param:long  longitude for the Point (in degrees. you can use positive and neg or from 0 to 360!)
4)  Param:alt   altitude for the Point (km - I think geographic alt as opposed to geomagnetic but I am not sure.)

***Optional Param:
5)  Param:user_ind  boolean value for geophysical indices. If false, then it grabs the indices from the date given. If true,
you must specify them (Kp and Ap index, etc.) 


Thus you can call a Point object by...

new_point = Point(date, latitude, longitude, altitude)
'''

                                # CALLING A SIMPLE POINT OBJECT
#**********************************************************************************************************************
#**********************************************************************************************************************

# To begin, lets first create a Point object from Jan 1st, 2000.

#-------------- Variables to change ---------------

date = datetime(2000, 1, 1, 0, 0)               # (yr, month, day, hour, minute) - refer to datetime module online
latitude = 45
longitude = 90
altitude = 400

#--------------------------------------------------

# for now we will not deal with inputting our out geophys indices, so we let the Point obj do this for us.
# Now make the Point object!

new_pt = Point(date, latitude, longitude, altitude)

# from this new_pt Point object, we can run MSIS (or anything else in the Point class).

result = new_pt.run_msis()          # result now contains the neutral temp, number densities (given in the run_msis fuction)
                                    # and total number density for the specific variables given above
'''
check them out...
which neutral species do you what to look at? 
Can be: 
HE, O, N2, O2, AR, H, N, O_anomalous, or rho (total number density) 

Or:
Tn_msis (neutral temperature [K])
'''

species = 'HE'      # <---- change this to whatever species you want!

density = result.nn[species]
total_density = result.rho
neutral_temp = result.Tn_msis

print('\n------Parameters------')
print('Date:', date)
print('Altitude [km]:', altitude)
print('Lat / Long [deg]:', latitude, '/', longitude)

print('\n------MSIS values------')
print(species, 'density [#/cm^3]:', density)
print('Total Density [#/cm^3]:', total_density)
print('Neutral Temperature [K]:', neutral_temp)


                                    # PLOTTING DENSITY VS ALTITUDE
#**********************************************************************************************************************
#**********************************************************************************************************************

'''
Now we want to check out the total number density of the atmosphere for a range of altitude values. We should see an
exponential decay.
'''

# make array of altitudes [km]
start_alt = 300
stop_alt = 600
stepsize = 10
alt_array = np.arange(start_alt, stop_alt, stepsize)

# make empty array of total density values to be initialized later. Declaring this now helps efficiency
tot_dens_arr = np.zeros(len(alt_array))

# Now specify your given lat and long and date. Since this is just a demo, I will keep the same values as above

# Run through altitude range and store all the MSIS total density values
for i in range(0, len(alt_array)):
    pt = Point(date, latitude, longitude, alt_array[i])
    result = pt.run_msis()
    tot_dens_arr[i] = result.rho

# Plot altitude vs density
plt.plot(tot_dens_arr, alt_array)
plt.xlabel("Total Density [#/cm^3]")
plt.ylabel('Altitude [km]')
plt.title((r'MSIS Total Number Density: ' + str(latitude) + '$\degree$N, ' + str(longitude) + '$\degree$E'))
plt.grid()
plt.show()

# we can also plot the number densities for all the individual species.
start_alt = 50
stop_alt = 600
stepsize = 5
alt_array = np.arange(start_alt, stop_alt, stepsize)

num_species = 7
dens_matrix = np.zeros((num_species, len(alt_array)))      # (row #, col #) each row is for one species (I did not do O_anomalous)

for i in range(0, len(alt_array)):
    pt = Point(date, latitude, longitude, alt_array[i])
    result = pt.run_msis()
    for j in range(0, num_species):
        if j == 0:
            dens_matrix[j, i] = result.nn['HE']
        elif j == 1:
            dens_matrix[j, i] = result.nn['O']
        elif j == 2:
            dens_matrix[j, i] = result.nn['N2']
        elif j == 3:
            dens_matrix[j, i] = result.nn['O2']
        elif j == 4:
            dens_matrix[j, i] = result.nn['AR']
        elif j == 5:
            dens_matrix[j, i] = result.nn['H']
        else:
            dens_matrix[j, i] = result.nn['N']

# Now we can plot it all.
for i in range(0, num_species):
    plt.semilogx(dens_matrix[i, :], alt_array)
plt.xlabel('Density [#/cm^3]')
plt.ylabel('Altitude [km]')
plt.legend(['He', 'O', 'N2', 'O2', 'AR', 'H', 'N'])
plt.title(('MSIS Total Number Density: ' + str(latitude) + '$\degree$N, ' + str(longitude) + '$\degree$E'))
plt.grid()
plt.show()


                                        # CONTOUR PLOTTING
#**********************************************************************************************************************
#**********************************************************************************************************************

'''
Now that we understand the basics of using MSIS, let's create a global contour plot of an output variable for every latitude
and longitude. For example, let's plot global temperature at an altitude of 400 km.
'''

# first lets define all our variables
date = datetime(1992, 3, 20, 0, 0)                      # using this date b/c is was low geomagnetic activity
gridsize = 2                                            # 2 degree grid size.
latitudes = np.arange(-90, 90, gridsize)                # lat and lon array values to pass into MSIS
longitudes = np.arange(-180, 180, gridsize)
altitude = 400
temp_matrix = np.zeros( (len(latitudes), len(longitudes)) )         # rows = lat, columns = lon


# Loop through lat and lon. select what MSIS values are desired
for i in range(0, len(longitudes)):
    for j in range(0, len(latitudes)):
        pt = Point(date, latitudes[j], longitudes[i], altitude)     # make Point and call MSIS
        result = pt.run_msis()
        temp_matrix[j, i] = result.Tn_msis                          # save the MSIS temp value

# Now plot the temperature matrix

# Create grid to plot data onto
x = np.arange(-180, 180, gridsize)
y = np.arange(-90, 90, gridsize)
X, Y = np.meshgrid(x, y)

# make a contour plot - refer to online for details
plt.figure()
levels = np.linspace(np.amin(temp_matrix), np.amax(temp_matrix), 100)
myplot = plt.contourf(X, Y, temp_matrix, levels,cmap='jet')
cont = plt.contour(X, Y, temp_matrix, 10, colors='k')
cbar = plt.colorbar(myplot)
cbar.ax.set_ylabel('Neutral Temp [K]')
plt.title('MSIS Neutral Temp: $Z_{gm}$=400km, UT=0')
plt.xlabel('Longitude [deg]')
plt.ylabel('Latitude [deg]')
plt.show()

