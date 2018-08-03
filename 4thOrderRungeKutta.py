############################################################################################
#
# Title: Fourth Order Runge-Kutta
#
# Purpose: To numerically integrate some ODE, given some initial condition and some time
#          step, to some final time.
#
# Author: Adam Farmer
# Sources: http://lpsa.swarthmore.edu/NumInt/NumIntFourth.html
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

## Initialize function to be numerically integrated, initial conditions and timestep.

#Initial condition.
t0 = 0.0 #Seconds
y0 = 0.0
#Final time.
tf = 10.0 #Seconds
#Time step.
dt = 0.01 #Seconds

#User defined function to be numerically integrated.
def dydt(y, t):
    dy = 2 * t #USER CHANGE THIS!
    return dy

#This information creates the time values list, an empty array for storing integrated
#values of dydt, and evaluates the function at the intial time given.
tVec = np.linspace(t0, tf, ((tf-t0)/dt))
ystar = np.zeros(len(tVec))

#### Impliment Runge-Kutta Method.
ystar[0] = y0 #Save initial condition into numerically integrated values.
for index, ti in enumerate(tVec[:-1]):
    # Define Runge-Kutta coefficients k1-k4
    k1 = dt*dydt(ystar[index], ti)
    k2 = dt*dydt(ystar[index] + k1*(dt/2), ti + dt/2)
    k3 = dt*dydt(ystar[index] + k2*(dt/2), ti + dt/2)
    k4 = dt*dydt(ystar[index] + k3*dt, ti + dt)
    
    #Save the resulting value after this timestep.
    ystar[index+1] = ystar[index] + ((k1 + 2*k2 + 2*k3 + k4)/6)
#    print(ti)
#    print(ystar[index])

## Plot the result and print the final value of the function
print(ystar[-1])
fig = plt.plot(tVec, ystar)
plt.show()


    
