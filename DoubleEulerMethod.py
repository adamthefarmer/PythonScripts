############################################################################################
#
# Title: Double Euler Method
#
# Purpose: To numerically integrate some ODE twice, given some initial conditions and some
#          time step, to some final time.
#
# About: This program impliments Euler's Method in order to numerically integrate a given
#        ODE twice. Two sets of intial conditions will be needed in order for the program
#        to function properly. The typical example that one might want to use this code
#        for is to integrate a given acceleration function to give the position at a later
#        time. 
#
# Author: Adam Farmer
# Sources: https://calcworkshop.com/first-order-differential-equations/eulers-method-table/
#
# Date Written:  6/22/18
# Date Modified: 6/22/18
#
############################################################################################

## Import libraries (if needed).

import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.axes as ax

#### Set user defined function to numerically integrate.

def d2ydt2(y, t):
    thing = 2*t # USER CHANGE THIS!
    return thing

#### Set initial condition, endpoint for integration.

t0 = 0.0      # Seconds
yPrime0 = 0.0 # Whatever
y0 = 0.0      # Whatever
dt = 0.00001  # Seconds
tf = 10       # Seconds
# Get all time values and save them for use later.
tVec = np.linspace(t0, tf, ((tf-t0)/dt))
yDoublePrimeStar = np.zeros(len(tVec))
yPrimeStar = np.zeros(len(tVec))
yStar = np.zeros(len(tVec))

#### Euler Method.

yPrimeStar[0] = yPrime0
yStar[0] = y0
for index, ti in enumerate(tVec[:-1]):
    # Evaluate the acceleration function (d2ydt2) at each point in time so we can plot it later.
    yDoublePrimeStar[index+1] = d2ydt2(0, ti)
    # Euler's Method the first time.
    yPrimeStar[index+1] = yPrimeStar[index] + dt*d2ydt2(0, ti)
    # Euler's Method the second time.
    yStar[index+1] = yStar[index] + dt*yPrimeStar[index]

#### Print and plot the results
    
# Print the final value of the velocity.
print(yPrimeStar[-1])
# Print the final value of the position.
print(yStar[-1])

# Turn on the minor TICKS, which are required for the minor GRID
plt.minorticks_on()
# Don't allow the axis to be on top of your data (couldn't get this to work :/)
#ax.set_axisbelow(True)

# Plotting the acceleration VS time
plt.figure(1)
plt.plot(tVec, yDoublePrimeStar)
plt.title('Acceleration VS Time')
plt.xlabel('Time (seconds)')
plt.ylabel('Acceleration (meters)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)

plt.figure(2)
# Plot of the position VS time.
plt.subplot(2, 1, 1)
plt.plot(tVec, yStar)
plt.title('Velocity and Position VS Time')
plt.ylabel('Position (meters)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)
# Plot of the velocity VS time.
plt.subplot(2, 1, 2)
plt.plot(tVec, yPrimeStar)
plt.xlabel('Time (seconds)')
plt.ylabel('Velocity (meters)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)

plt.show()
