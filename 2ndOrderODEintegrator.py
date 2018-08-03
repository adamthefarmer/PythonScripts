############################################################################################
#
# Title: Numerical Integration for 2nd Order ODE w/ 1st Derivative Term
#
# Purpose: The purpose of this program is to allow the user to numerically integrate
#          a second order ODE which also has a 1st order term which is not explicitely
#          known. See the 'Background' section on line 27 for more information on how
#          to prepare the program for your use.
#
#
# Author: Adam Farmer
# Sources: https://math.oregonstate.edu/home/programs/undergrad/CalculusQuestStudyGuides/ode/second/so_num/so_num.html
#
# Date Written:  6/26/18
# Date Modified: 6/26/18
#
############################################################################################

## Import libraries.

import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.axes as ax

#### Background

# Originally we have an equation which looks like the following:
# a*theta''(t) + b*theta'(t) + f(theta(t))
# Define: omega(t) = theta'(t)
# Use this to write the second order ODE as two coupled, first order ODEs
# omega(t) = theta'(t)  and  omega'(t) = -1/a [b*omega(t) + f(theta(t))]
# Then use this special implimentation of Euler's Method to numerically integrate.

#### Set the two ODEs here

def getOmegaDot(omega_i, theta_i):
    
    # Set the coefficients 'a' and 'b'
    a = 1
    b = 0.3

    omegaDot = -(1/a) * (b*omega_i + math.sin(theta_i))
    return omegaDot

#### Give time step, initial state vector and so on...

# Time stuff
t0 = 0.0
dt = 0.0001
tf = 30.0
tVec = np.linspace(t0, tf, ((tf-t0)/dt))
# Initial state vectors.
theta_i = 1.0
omega_i = 0.0
# Initialize some lists to hold angular rate and position values for the system.
theta = np.zeros(len(tVec))
omega = np.zeros(len(tVec))
theta[0] = theta_i
omega[0] = omega_i

#### Do the thing

for i, ti in enumerate(tVec[:-1]):
    # Basically Euler's Method.
    theta[i+1] = theta[i] + dt*omega[i]
    omega[i+1] = omega[i] + dt*getOmegaDot(omega[i], theta[i])

#### Plot the results
    
# Turn on the minor TICKS, which are required for the minor GRID.
plt.minorticks_on()
# Don't allow the axis to be on top of your data. (couldn't get this to work :/)
#ax.set_axisbelow(True)

plt.figure(1)
# Plot of the position VS time.
plt.subplot(2, 1, 1)
plt.plot(tVec, theta)
plt.title('Angular Velocity and Position VS Time')
plt.ylabel('Position (degrees)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)
# Plot of the velocity VS time.
plt.subplot(2, 1, 2)
plt.plot(tVec, omega)
plt.xlabel('Time (seconds)')
plt.ylabel('Velocity (deg/sec)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)

plt.show()

