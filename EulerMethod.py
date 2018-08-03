############################################################################################
#
# Title: Euler Method
#
# Purpose: To numerically integrate some ODE, given some initial condition and some time
#          step, to some final time.
#
# Author: Adam Farmer
# Sources: 
#
# Date Written:  6/21/18
# Date Modified: 6/21/18
#
############################################################################################

## Import libraries (if needed).

import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

## Set user defined function to numerically integrate.

def dydt(y, t):
    thing = 2*t
    return thing

## Set initial condition, endpoint for integration.

t0 = 0.0 #Seconds
y0 = 0.0 #Whatever
dt = 0.001 #Seconds
tf = 10 #Seconds
#Get all time values and save them for use later.
tVec = np.linspace(t0, tf, ((tf-t0)/dt))
ystar = np.zeros(len(tVec))

## Euler Method.
ti = t0
ystar[0] = y0
for index, ti in enumerate(tVec[:-1]):
    ystar[index+1] = ystar[index] + dt*dydt(0, ti)

## Print and plot the results
print(ystar[-1])
fig = plt.plot(tVec, ystar)
plt.show()

    
