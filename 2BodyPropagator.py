############################################################################################
#
# Title: Two-Body Orbital Propagator
#
# Purpose: 
#
# About:  
#
# Author: Adam Farmer
# Sources: 
#
# Date Written:  6/29/18
# Date Modified: 6/29/18
#
############################################################################################

## Import libraries.

import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.axes as ax

#### Define some constants used in this program.

mEarth = 5.972 * (10**24) # kg
rEarth = 6.371 * (10**6) # meters
G_Earth = 6.67408 * (10**(-11)) # m^3 kg^-1 s^-2

#### Set intital conditions.

t0 = 0 # Seconds
dt = .01 # Seconds
tf = 10000 # Seconds
# Make time vector.
tVec = np.linspace(t0, tf, ((tf-t0)/dt))
# Initial position of the spacecraft.
ri = [2400000+rEarth, 0 , 0] # Meters
# This is the magnitude of the needed velocity to stay in a stable LEO.
V_magLEO = m.sqrt((G_Earth*mEarth)/(m.sqrt(ri[0]**2 + ri[1]**2 +ri[2]**2)))
# Initial velocity of the spacecraft.
vi = [0, 0, V_magLEO + 5000]
# Make position, velocity and acceleration vectors for the ECI frame and put intitial conditions in them.
r_ECI = np.zeros((len(tVec),3))
r_ECI[0,:] = ri
v_ECI = np.zeros((len(tVec),3))
v_ECI[0,:] = vi
a_ECI = np.zeros((len(tVec),3))
B_ECI = np.zeros((len(tVec),3))

#### Define functions used in this program.

def cartNorm(vec):
    # Take the cartesian norm of the input three element vector.
    result = m.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
    
    return result

def getUnitVec(vec):
    # returns the unit vector of the input vector.
    unitVec = vec / cartNorm(vec)

    return unitVec

def geta_ECI(r_ECI):
    # Using Newton's Law of Gravitation, find the acceleration vector of the spacecraft at
    # the current timestep.
    
    # Get the magnitude of the acceleration vector
    a_ECImag = (G_Earth*mEarth)/(cartNorm(r_ECI)**2)
    # Get the unit vector of the acceleration by noting that it is always antiparallel to the position vector.
    a_ECIhat = -getUnitVec(r_ECI)
    # Construct the final acceleration vector.
    a_ECI = a_ECImag*a_ECIhat
    
    return a_ECI

def getBvec(r_ECI, B_ECI):
    # Given some Cartesian coordinate vector in the ECI frame, this function will return the
    # approximate magnetic field vector according to the ECI frame.

    # Get the angle phi in spherical coordinates in the ECI frame.
    phi = m.acos((r_ECI[2])/(cartNorm(r_ECI)))

    # Get the components
    B_ECI[0] = -r_ECI[0]
    B_ECI[1] = -r_ECI[1]
    B_ECI[2] = m.cos(2*phi)

    # Normalize the result to get the unit vector.
    B_ECI = B_ECI / cartNorm(B_ECI)

    return B_ECI

#### Propagation of velocity and position.

for index, ti in enumerate(tVec[:-1]):
    # Evaluate the acceleration function (d2ydt2) at each point in time so we can plot it later.
    a_ECI[index+1,:] = geta_ECI(r_ECI[index,:])
    # Euler's Method the first time.
    v_ECI[index+1,:] = v_ECI[index,:] + dt*a_ECI[index,:]
    # Euler's Method the second time.
    r_ECI[index+1] = r_ECI[index] + dt*v_ECI[index,:]
    # Get Earth's magnetic field vector for the current timestep.
    B_ECI[index,:] = getBvec(r_ECI[index,:], B_ECI[index,:])


#### Plotting the results

# Turn on the minor TICKS, which are required for the minor GRID.
#plt.minorticks_on()
# Don't allow the axis to be on top of your data. (couldn't get this to work :/)
#ax.set_axisbelow(True)


mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure(1)


ax = fig.gca(aspect='equal',projection='3d')

theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
z = np.linspace(-2, 2, 100)
r = z**2 + 1
x = r * np.sin(theta)
y = r * np.cos(theta)
ax.plot(r_ECI[:,0]/1000, r_ECI[:,1]/1000, r_ECI[:,2]/1000, label='Orbital Position', color='k')

# Sphere to represent Earth
# Make data
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x1 = rEarth/1000 * np.outer(np.cos(u), np.sin(v))
y1 = rEarth/1000 * np.outer(np.sin(u), np.sin(v))
z1 = rEarth/1000 * np.outer(np.ones(np.size(u)), np.cos(v))
# Plot the surface
ax.plot_surface(x1, y1, z1, cmap='GnBu')


# Legend and labels
ax.legend()
ax.set_xlabel('X Pos (km)')
ax.set_ylabel('Y Pos (km)')
ax.set_zlabel('Z Pos (km)')
ax.set_aspect('equal')

plt.show()




