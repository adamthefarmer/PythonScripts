############################################################################################
#
# Title: Passive ACS Simulation
#
# Purpose: 
#
# About: 
#
# Author: Adam Farmer
# Sources: http://lasp.colorado.edu/home/csswe/files/2012/06/Gerhardt_SSC10_PMAC.pdf
#          https://slideplayer.com/slide/5966557/
#
# Date Written:  7/2/18
# Date Modified: 7/2/18
#
############################################################################################

## Import libraries.

import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.axes as ax
import mpl_toolkits
from mpl_toolkits.mplot3d import Axes3D

#### Define some constants used in this program.

## Orbit propagation
mEarth = 5.972 * (10**24) # kg
rEarth = 6.371 * (10**6) # meters
G_Earth = 6.67408 * (10**(-11)) # m^3 kg^-1 s^-2


# Earth's magnetic field in LEO varies between .2 and .45 Gauss
# We can choose that B = approx 3.2 * 10^ -5 Tesla
Mz = 8*(10**15) #T*m^3

## Hysteresis rod stuff
mu_0 = 4*np.pi*(10**(-7))
mu_hyst = 1.5*(10**4)
Length_hyst = 0.095 # Length of the hysteresis rod, meters
Diameter_hyst = 0.001 # Diameter of hysteresis rod, meters
V_hyst = np.pi*((Diameter_hyst/2)**2)*Length_hyst # Volume of the hysteresis rod, m^3
m_hyst = np.zeros(3)
Br = 0.35 # Tesla
Bs = 0.74 # Tesla
Hc = 0.96 # A/m
numHyst = [2.0, 2.0, 0.0] # Define the number of hysteresis rods on each axis.
B_hyst = np.array([0.0, 0.0, 0.0])
m_hyst = np.array([0.0, 0.0, 0.0])

#### Definitions of functions used in this program

# Sine and cosine functions, 
def s(angle):
    result = np.sin(np.deg2rad(angle))
##    print(angle)
##    if angle > 3e6:
##        raise Exception("WOOOOAH")
    return result

def c(angle):
    result = np.cos(np.deg2rad(angle))
    return result

def normalizeAngle(angle):
    # Take some angle in degrees and normalize it to a value between 0 and 360.
    normalizedAngle = angle % 360.0
    #normalizedAngle = angle % (2*np.pi)

    return normalizedAngle

def DCMrotation(theta, vector):
    # Take a vector in the ECI frame and use the DCM to describe it in terms of the body frame.
    # Need current theta passed in to get the correct DCM

    DCM = np.array([[c(theta[1])*c(theta[0]), c(theta[1])*s(theta[0]), -s(theta[1])],
             [s(theta[2])*s(theta[1])*c(theta[0]) - c(theta[2])*s(theta[0]), s(theta[2])*s(theta[1])*s(theta[0]) + c(theta[2])*c(theta[0]), s(theta[2])*c(theta[1])],
             [c(theta[2])*s(theta[1])*c(theta[0]) + s(theta[2])*s(theta[0]), c(theta[2])*s(theta[1])*s(theta[0]) - s(theta[2])*c(theta[0]), c(theta[2])*c(theta[1])]])

    # Do the rotation.
    result = np.matmul(DCM, vector)
    
    return result

def getThetaDot(theta, omega, index, lastThetaDot):
    # Get the Euler angular rates given a current body fixed angular rate and current Euler angle.

    # From top of pg.4 we have the rotational kinematics equations:
    correctionMatrix = np.array([[0, s(theta[2]), c(theta[2])],
                     [0, c(theta[2])*c(theta[1]), -s(theta[2])*c(theta[1])],
                     [c(theta[1]), s(theta[2])*c(theta[1]), c(theta[2])*s(theta[1])]])

    ## Make sure we don't hit the singularity at theta[1] = +-90deg
    if (normalizeAngle(theta[1]) - 90) < 0.5 and (normalizeAngle(theta[1]) - 90) > 0 and lastThetaDot[1] > 0: # case for theta greater than zero
        #print("why tho?", normalizeAngle(theta[1]) - 90)
        thetaDot = 10 * (np.matmul(correctionMatrix, omega))
        #print("case for theta greater than zero: ", index, theta[1])
    elif (normalizeAngle(theta[1]) - 270) < 0.5 and (normalizeAngle(theta[1]) - 270) > 0 and lastThetaDot[1] < 0: # case for theta less than zero
        thetaDot = -10 * (np.matmul(correctionMatrix, omega))
        #print("case for theta less than zero: ", index, theta[1])
    elif theta[1] is 0: # case for theta = 0
        if index is 0:
            raise Exception("Theta[1] may have initial conditions at a sigularity value. This is not recommended.")
            thetaDot = 1/c(theta[1]) * (np.matmul(correctionMatrix, omega))
        else:
            sign = np.sign(thetaDot[index-1])
            thetaDot = sign*10*(np.matmul(correctionMatrix, omega))
            #print("case for theta equal to zero: ", index)
    else:
        thetaDot = 1/c(theta[1]) * (np.matmul(correctionMatrix, omega))
        #print("heck")
        
        
    return thetaDot

def getExternalTorque(B_current, B_last, theta_current, L_current, index):
    global B_hyst
    global m_hyst
    # Get the external torques on the spacecraft based on current attitude. Currently the external torques on the spacecraft
    # are the torque produced by the interaction of the on-board bar magnet with the Earth's magnetic field, and the damping
    # torque produced by the changing magnetic field inducing currents in the hysteresis rods.
    
    # Since the magnetic diple moment is always aligned with the z-axis this will never change.
    # The torque produced by the bar magnet in Earth's B field is then given by: tao = m cross B
    L_current[0:3] = np.cross(m_bar,B_current)

    # The damping torque produced by the hysteresis rods is described on the bottom of pg.3 in the paper cited.
    
    p = (1/Hc)*np.tan((np.pi*Br)/(2*Bs)) # Constant, same for all rods.
    
    sign = np.array([1.0, 1.0, 1.0])
    # X-Axis
    if B_last[0] < B_current[0]: # Mag field strength is increasing.
        sign[0] = -1.0
    else:
        sign[0] = 1.0
    # Y-Axis        
    if B_last[1] < B_current[1]: # Mag field strength is increasing.
        sign[1] = -1.0
    else:
        sign[1] = 1.0
    # Z-Axis        
    if B_last[2] < B_current[2]: # Mag field strength is increasing.
        sign[2] = -1.0
    else:
        sign[2] = 1.0

    B_hyst[0] = (2/np.pi)*Bs*np.arctan(p*(B_current[0] + sign[0]*Hc))
    B_hyst[1] = (2/np.pi)*Bs*np.arctan(p*(B_current[1] + sign[1]*Hc))
    B_hyst[2] = (2/np.pi)*Bs*np.arctan(p*(B_current[2] + sign[2]*Hc))

    # Hysteresis rod magnetic moment components.
    m_hyst = numHyst * ((B_hyst * V_hyst) / mu_0)
    
    L_current[3:] = np.cross(m_hyst,B_current)
    
    return L_current
    
    
def getOmegaDot(omega, L, index):
    # Get the change in the body fixed angular rates WRT time.
##    omegaDot = np.array([1/I[0]*(-(I[2]-I[1])*omega[1]*omega[2] + L[0] - L[3]),
##                        1/I[1]*(-(I[0]-I[2])*omega[2]*omega[0] + L[1] - L[4]),
##                        1/I[2]*(-(I[1]-I[0])*omega[0]*omega[1] + L[2] - L[5])])

##    omegaDot = np.array([1/I[0]*(-(I[2]-I[1])*omega[1]*omega[2]),
##                        1/I[1]*(-(I[0]-I[2])*omega[2]*omega[0]),
##                        1/I[2]*(-(I[1]-I[0])*omega[0]*omega[1])])

    omegaDot = np.array([1/I[0]*(-(I[2]-I[1])*omega[1]*omega[0]),
                        1/I[1]*(-(I[0]-I[2])*omega[0]*omega[2]),
                        1/I[2]*(-(I[1]-I[0])*omega[2]*omega[1])])

    return omegaDot

def cartNorm(vec):
    # Take the cartesian norm of the input three element vector.
    result = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
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

def getBvec(r_ECI, theta, B, index):
    # Given some Cartesian coordinate vector in the ECI frame, this function will return the
    # approximate magnetic field vector according to the ECI frame and also the body frame.

    # Good approximation to Earth's magnetic field
    B[0] = -Mz*(3*r_ECI[0]*r_ECI[2]) / (((cartNorm(r_ECI))**5))
    B[1] = -Mz*(3*r_ECI[1]*r_ECI[2]) / (((cartNorm(r_ECI))**5))
    B[2] = Mz*(3*(r_ECI[2]**2) - (((cartNorm(r_ECI))**2))) / (((cartNorm(r_ECI))**5))

    ## Also use the current DCM to get the B vector in the body frame as well
    # Define the DCM
    DCM = np.array([[c(theta[1])*c(theta[0]), c(theta[1])*s(theta[0]), -s(theta[1])],
             [s(theta[2])*s(theta[1])*c(theta[0]) - c(theta[2])*s(theta[0]), s(theta[2])*s(theta[1])*s(theta[0]) + c(theta[2])*c(theta[0]), s(theta[2])*c(theta[1])],
             [c(theta[2])*s(theta[1])*c(theta[0]) + s(theta[2])*s(theta[0]), c(theta[2])*s(theta[1])*s(theta[0]) - s(theta[2])*c(theta[0]), c(theta[2])*c(theta[1])]])

    B[3:] = np.matmul(DCM, B[:3])

    return B

#### Initial state vectors, time and whatnot

## Position propagation stuff:
# Initial position of the spacecraft.
ri = [600000+rEarth, 0.0 , 0.0] # Meters

# This is the magnitude of the needed velocity to stay in a stable LEO.
V_magLEO = np.sqrt((G_Earth*mEarth)/(np.sqrt(ri[0]**2 + ri[1]**2 +ri[2]**2)))
# Initial velocity of the spacecraft.
vi = [0.0, 0.0, V_magLEO]

## Time stuff
# Find the orbital period
orbitalPeriod = 2*np.pi*cartNorm(ri)/V_magLEO
t0 = 0.0 # Seconds
dt = 0.1 # Seconds
tf = orbitalPeriod # Seconds
tVec = np.linspace(t0, tf, ((tf-t0)/dt)) # Make time vector.

# Make position, velocity and acceleration vectors for the ECI frame and put intitial conditions in them.
r_ECI = np.zeros((len(tVec),3))
r_ECI[0,:] = ri
v_ECI = np.zeros((len(tVec),3))
v_ECI[0,:] = vi
a_ECI = np.zeros((len(tVec),3))
B = np.zeros((len(tVec),6))

## Attitude propagation stuff:
# Moments of inertia about the spacecraft's principle axes.
I = [0.1, 0.1, 0.1]
# Bipole moment of the on-board magnet. Units are A*m^2
m_bar = [0.0, 0.0, 0.29]
# Initialize external torques array.
L = np.zeros((len(tVec),6))
# Initial state vectors as 1, 2, 3.
theta_i = np.array([0.0, 0.0, 0.0])  # Initial Euler anglar postion
omega_i = np.array([0.0, 0.0, 0.0])  # Initial body fixed angular velocity in degrees/sec
# Initialize some lists to hold angular rate and position values for the system.
theta = np.zeros((len(tVec),3))
omega = np.zeros((len(tVec),3))
thetaDot = np.zeros((len(tVec),3))
omegaDot = np.zeros((len(tVec),3))
# Put initial state vectors into the arrays for angular rate and position.
theta[0,:] = theta_i
omega[0,:] = omega_i
thetaDot[0,:] = getThetaDot(theta[0,:], omega[0,:], 0, [0.0, 0.0, 0.0])
omegaDot[0,:] = getOmegaDot(omega[0,:], L[0,:], 0)

#### Propagation of velocity, position and attitude.

for index, ti in enumerate(tVec[:-1]):

    ## Position Stuff:
    # Evaluate the acceleration function (d2ydt2) at each point in time so we can plot it later.
    a_ECI[index+1,:] = geta_ECI(r_ECI[index,:])
    # Euler's Method the first time.
    v_ECI[index+1,:] = v_ECI[index,:] + dt*a_ECI[index,:]
    # Euler's Method the second time.
    r_ECI[index+1] = r_ECI[index] + dt*v_ECI[index,:]
    # Get Earth's magnetic field vector for the current timestep. Reports in both the body and ECI frame!
    B[index,:] = getBvec(r_ECI[index,:], theta[index,:], B[index,:], index)

    ## Attitude stuff: Get a new Euler angular rate for this time step
    # If this is the first timestep, dB/dt = 0
    if index == 0:
        L[index,:] = getExternalTorque(B[index,3:], B[index,3:], theta[index,:], L[index,:], index)
        thetaDot[index,:] = getThetaDot(theta[index,:], omega[index,:], index, [0.0, 0.0, 0.0])
    else:
        L[index,:] = getExternalTorque(B[index,3:], B[index-1,3:], theta[index,:], L[index,:], index)
        thetaDot[index,:] = getThetaDot(theta[index,:], omega[index,:], index, thetaDot[index-1,:])
    omegaDot[index,:] = getOmegaDot(omega[index,:], L[index,:], index)
    
    ## Euler's Method to propagate theta and omega
    #theta[index+1,:] = normalizeAngle(theta[index,:] + dt*thetaDot[index,:])
    theta[index+1,:] = theta[index,:] + dt*thetaDot[index,:]
    omega[index+1,:] = omega[index,:] + dt*omegaDot[index,:]

    if thetaDot[index,0] > 1800 or thetaDot[index,1] > 1800 or thetaDot[index,2] > 1800 or thetaDot[index,0] < -1800 or thetaDot[index,1] < -1800 or thetaDot[index,2] < -1800:
        print(index)
        print("Omega Before: {}".format(omega[index-1,:]))
        print("Omega Current: {}".format(omega[index,:]))
        print("Omega After: {}".format(omega[index+1,:]))
        print("ThetaDot Before Before: {}".format(thetaDot[index-2,:]))
        print("ThetaDot Before: {}".format(thetaDot[index-1,:]))
        print("ThetaDot Current: {}".format(thetaDot[index,:]))
        print("ThetaDot After: {}".format(thetaDot[index+1,:]))
        print("Theta Before Before: {}".format(theta[index-2,:]))
        print("Theta Before: {}".format(theta[index-1,:]))
        print("Theta Current: {}".format(theta[index,:]))
        print("Theta After: {}".format(theta[index+1,:]))
    
    #print("Ang. mom.: {}".format(omega[index,:]*I))
    #print("Mag or ang. mom.: {}".format(cartNorm(omega[index+1,:]*I)))

    ## Runge-Kutta Method for numerical integration:
##
##    # Theta Integration
##    # Define Runge-Kutta coefficients k1-k4 for the 1st integration.
##    k1 = thetaDot[index]
##    k2 = 
    
##    k1 = d2ydt2(yPrimeStar[index], ti)
##    k2 = d2ydt2(yPrimeStar[index] + k1*(dt/2), ti + dt/2)
##    k3 = d2ydt2(yPrimeStar[index] + k2*(dt/2), ti + dt/2)
##    k4 = d2ydt2(yPrimeStar[index] + k3*dt, ti + dt)
##    #Save the resulting value after this timestep.
##    yPrimeStar[index+1] = yPrimeStar[index] + ((k1 + 2*k2 + 2*k3 + k4)*(dt/6))

#quit()

##for i in range(0, 7):
##    print(i + 18176)
##    print("Omega: {}".format(omega[i+18176,:]))
##    print("ThetaDot: {}".format(thetaDot[i+18176,:]))
##    print("Theta: {}".format(theta[i+18176,:]))

#### Plotting stuff.

# Turn on the minor TICKS, which are required for the minor GRID.
mpl.rcParams['legend.fontsize'] = 10

## Position plotting:
fig = plt.figure(1)
ax = fig.gca(aspect='equal',projection='3d')

ax.plot(r_ECI[:,0]/1000, r_ECI[:,1]/1000, r_ECI[:,2]/1000, label='Orbital Position', color='k')

# Sphere to represent Earth
# Make data
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
xSphere = rEarth/1000 * np.outer(np.cos(u), np.sin(v))
ySphere = rEarth/1000 * np.outer(np.sin(u), np.sin(v))
zSphere = rEarth/1000 * np.outer(np.ones(np.size(u)), np.cos(v))
# Plot the surface
ax.plot_surface(xSphere, ySphere, zSphere, cmap='GnBu')

# Legend and labels.
ax.legend()
ax.set_xlabel('X Pos (km)')
ax.set_ylabel('Y Pos (km)')
ax.set_zlabel('Z Pos (km)')
plt.title('Orbit Plotting')
#ax.set_aspect('equal')

## Attitude Plotting:
plt.figure(2)
plt.plot(tVec, theta[:,0])
plt.plot(tVec, theta[:,1])
plt.plot(tVec, theta[:,2])
plt.title('Angular Position VS Time in the ECI Frame')
plt.legend(['x-axis','y-axis','z-axis'])
plt.ylabel('Ang. Pos. (degrees)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)

## Angular velocity in the body & global frame plotting:
plt.figure(3)
plt.minorticks_on()
# Plot of the agular velocity of the spacecraft in the body frame VS time.
plt.subplot(2, 1, 1)
plt.plot(tVec, omega[:,0])
plt.plot(tVec, omega[:,1])
plt.plot(tVec, omega[:,2])
plt.title('Angular Velocity: ECI & Body Frame VS Time')
plt.legend(['z-axis','y-axis','x-axis'])
plt.ylabel('Ang. Vel. Body (deg/sec)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)
# Plot of the angular velocity of the spacecraft in the ECI frame VS time.
plt.subplot(2, 1, 2)
plt.plot(tVec, thetaDot[:,0])
plt.plot(tVec, thetaDot[:,1])
plt.plot(tVec, thetaDot[:,2])
plt.legend(['x-axis','y-axis','z-axis'])
plt.xlabel('Time (seconds)')
plt.ylabel('Ang. Vel. ECI (deg/sec)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)

## Plotting the magnetic field
plt.figure(4)
plt.minorticks_on()
# Plot of B in the ECI frame VS time.
plt.subplot(2, 1, 1)
plt.plot(tVec, B[:,0:3])
plt.title('Magnetic Field in the ECI Frame')
plt.legend(['x-axis','y-axis','z-axis'])
plt.ylabel('B_ECI (Tesla)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)
# Plot of B in the body frame VS time.
plt.subplot(2, 1, 2)
plt.plot(tVec, B[:,3:])
plt.legend(['x-axis','y-axis','z-axis'])
plt.xlabel('Time (seconds)')
plt.ylabel('B_body (Tesla)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)

## Plotting the external torques
plt.figure(5)
plt.minorticks_on()
# Plot of external torque from the bar magnet VS time.
plt.subplot(2, 1, 1)
plt.plot(tVec, L[:,0:3])
plt.title('External Torques on Spacecraft VS Time')
plt.legend(['x-axis','y-axis','z-axis'])
plt.ylabel('Bar Magnet Torque (N*m)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)
# Plot of B in the body frame VS time.
plt.subplot(2, 1, 2)
plt.plot(tVec, L[:,3:])
plt.legend(['x-axis','y-axis','z-axis'])
plt.xlabel('Time (seconds)')
plt.ylabel('Hysteresis Torque (N*m)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)

plt.show()
   


