############################################################################################
#
# Title: Double Numerical Integration w/ Fourth Order Runge-Kutta Method & Euler's Method
#
# Purpose: To numerically integrate some ODE twice, given some initial and final time and an
#          intial state vector which describes the system at that intial time. The classical
#          example for which this program would be useful for is to numercally integrate an
#          objects known acceleration VS time to know that objects velocity and postion at
#          a later point in time with reasonable accuracy. This program uses the 4th order
#          Runge-Kutta Method for the 1st integration, and then uses Euler's Method for the
#          second integration.
#
# Author: Adam Farmer
# Sources: http://lpsa.swarthmore.edu/NumInt/NumIntFourth.html
#          https://calcworkshop.com/first-order-differential-equations/eulers-method-table/
#
# Date Written:  6/23/18
# Date Modified: 6/23/18
#
############################################################################################

## Import libraries (if needed).

import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#### Initialize function to be numerically integrated, initial conditions and timestep.

#Initial condition.
t0 = 0.0 #Seconds
yPrime0 = 0.0
y0 = 0.0
#Final time.
tf = 10.0 #Seconds
#Time step.
dt = 0.00001 #Seconds

#User defined function to be numerically integrated.
def d2ydt2(y, t):
    thing2 = 2*t #USER CHANGE THIS!
    return thing2

#This information creates the time values list, an empty array for storing integrated
#values of dydt, and evaluates the function at the intial time given.
tVec = np.linspace(t0, tf, ((tf-t0)/dt))
yPrimeStar = np.zeros(len(tVec))
yStar = np.zeros(len(tVec))

#### Impliment Runge-Kutta Method.

#Save initial condition into numerically integrated values.
yPrimeStar[0] = yPrime0
yStar[0] = y0
for index, ti in enumerate(tVec[:-1]):
    # Define Runge-Kutta coefficients k1-k4 for the 1st integration.
    k1 = d2ydt2(yPrimeStar[index], ti)
    k2 = d2ydt2(yPrimeStar[index] + k1*(dt/2), ti + dt/2)
    k3 = d2ydt2(yPrimeStar[index] + k2*(dt/2), ti + dt/2)
    k4 = d2ydt2(yPrimeStar[index] + k3*dt, ti + dt)
    #Save the resulting value after this timestep.
    yPrimeStar[index+1] = yPrimeStar[index] + ((k1 + 2*k2 + 2*k3 + k4)*(dt/6))

    # Use Euler's Method for the 2nd integration.
    yStar[index+1] = yStar[index] + dt*yPrimeStar[index]


#### Plot the result and print the final value of the function

# Print the final values of velocity and position.
print(yPrimeStar[-1])   
print(yStar[-1])

# Turn on the minor TICKS, which are required for the minor GRID
plt.minorticks_on()
# Don't allow the axis to be on top of your data (couldn't get this to work :/)
#ax.set_axisbelow(True)

plt.figure(1)
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

    
