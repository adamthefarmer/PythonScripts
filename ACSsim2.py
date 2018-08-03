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

def cartNorm(vec):
    # Take the cartesian norm of the input three element vector.
    result = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
    return result

def normalizeAngle(angle):
    # Take some angle in degrees and normalize it to a value between 0 and 360.
    normalizedAngle = angle % 360.0
    #normalizedAngle = angle % (2*np.pi)

    return normalizedAngle

def getDCM(EP):

    DCMcurrent = np.array([[EP[0]**2 + EP[1]**2 - EP[2]**2 - EP[3]**2, 2*(EP[1]*EP[2] + EP[0]*EP[3]), 2*(EP[1]*EP[3] - EP[0]*EP[2])],
                    [2*(EP[1]*EP[2] - EP[0]*EP[3]), EP[0]**2 - EP[1]**2 + EP[2]**2 - EP[3]**2, 2*(EP[2]*EP[3] + EP[0]*EP[1])],
                    [2*(EP[1]*EP[3] + EP[0]*EP[2]), 2*(EP[2]*EP[3] - EP[0]*EP[2]), EP[0]**2 - EP[1]**2 - EP[2]**2 + EP[3]**2]])

    return DCMcurrent

######################################
# Checks if a matrix is a valid rotation matrix.
def isRotationMatrix(DCM) :
    DCMt = np.transpose(DCM)
    shouldBeIdentity = np.dot(DCMt, DCM)
    I = np.identity(3, dtype = DCM.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6
 
# Calculates rotation matrix to euler angles
def DCMtoEulerAngles(DCM):
 
    #assert(isRotationMatrix(DCM))
     
    sy = np.sqrt(DCM[0,0] * DCM[0,0] +  DCM[1,0] * DCM[1,0])
     
    singular = sy < 1e-6
 
    if  not singular :
        x = np.arctan2(DCM[2,1] , DCM[2,2])
        y = np.arctan2(-DCM[2,0], sy)
        z = np.arctan2(DCM[1,0], DCM[0,0])
    else :
        x = np.arctan2(-DCM[1,2], DCM[1,1])
        y = np.arctan2(-DCM[2,0], sy)
        z = 0
 
    return np.array([x, y, z])
####################################

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

def getEPdot(EP, omega, index):
    # Get the Euler angular rates given a current body fixed angular rate and current Euler angle.

    # From the coupled kinematic differential equations of motion for Euler parameters (quaternions) pg.100 of Analytical Mechanics of Aerospace Systems
    correctionMatrix = np.array([[EP[0], -EP[1], -EP[2], -EP[3]],
                                 [EP[1], EP[0], -EP[3], EP[2]],
                                 [EP[2], EP[3], EP[0], -EP[1]],
                                 [EP[3], -EP[2], EP[1], EP[0]]])

    # Needed to add an extra zero in front of the omega vector.
    tempOmega = np.array([0, omega[0], omega[1], omega[2]])

    # pg.100 of Analytical Mechanics of Aerospace Systems
    EPdot = 0.5*(np.matmul(correctionMatrix, tempOmega))

    return EPdot
    
def getOmegaDot(omega, L, index):
    # Get the change in the body fixed angular rates WRT time.
##    omegaDot = np.array([1/I[0]*(-(I[2]-I[1])*omega[1]*omega[2] + L[0] - L[3]),
##                         1/I[1]*(-(I[0]-I[2])*omega[2]*omega[0] + L[1] - L[4]),
##                         1/I[2]*(-(I[1]-I[0])*omega[0]*omega[1] + L[2] - L[5])])

    omegaDot = np.array([1/I[0]*(-(I[2]-I[1])*omega[1]*omega[2]),
                         1/I[1]*(-(I[0]-I[2])*omega[2]*omega[0]),
                         1/I[2]*(-(I[1]-I[0])*omega[0]*omega[1])])

    return omegaDot

def getUnitVec(vec):
    # returns the unit vector of the input vector.
    unitVec = vec / cartNorm(vec)

    return unitVec

def geta_ECI(r_ECI):
    # Using Newton's Law of Gravitation, find the acceleration vector of the spacecraft at the current timestep.
    
    # Get the magnitude of the acceleration vector
    a_ECImag = (G_Earth*mEarth)/(cartNorm(r_ECI)**2)
    # Get the unit vector of the acceleration by noting that it is always antiparallel to the position vector.
    a_ECIhat = -getUnitVec(r_ECI)
    # Construct the final acceleration vector.
    a_ECI = a_ECImag*a_ECIhat

    return a_ECI


def getBvec(r_ECI, B, index):
    # Given some Cartesian coordinate vector in the ECI frame, this function will return the
    # approximate magnetic field vector according to the ECI frame and also the body frame.

    # Good approximation to Earth's magnetic field
##    B[0] = -Mz*(3*r_ECI[0]*r_ECI[2]) / (((cartNorm(r_ECI))**5))
##    B[1] = -Mz*(3*r_ECI[1]*r_ECI[2]) / (((cartNorm(r_ECI))**5))
##    B[2] = Mz*(3*(r_ECI[2]**2) - (((cartNorm(r_ECI))**2))) / (((cartNorm(r_ECI))**5))

    # Trivial case
    B[0] = 2
    B[1] = 0
    B[2] = 0

    ## Also use the current DCM to get the B vector in the body frame as well

    # Get the DCM for the current EP.
    DCM = getDCM(EP[index,:])
    # Rotate the magnetic field vector from the ECI frame into the body frame.
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
EP_i = np.array([1.0, 0.0, 0.0, 0.0])  # Initial Euler parameters.
omega_i = np.array([.05, 0.0, 0.0])  # Initial body fixed angular velocity in degrees/sec
# Initialize some lists to hold angular rate and position values for the system.
theta = np.zeros((len(tVec),3))
EP = np.zeros((len(tVec),4))
omega = np.zeros((len(tVec),3))
EPdot = np.zeros((len(tVec),4))
omegaDot = np.zeros((len(tVec),3))
# Put initial state vectors into the arrays for angular rate and position.
EP[0,:] = EP_i
omega[0,:] = omega_i
EPdot[0,:] = getEPdot(EP[0,:], omega[0,:], 0)
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
    B[index,:] = getBvec(r_ECI[index,:], B[index,:], index)

    ## Attitude stuff:
    if index == 0: # If this is the first timestep, dB/dt = 0
        L[index,:] = getExternalTorque(B[index,3:], B[index,3:], EP[index,:], L[index,:], index)
    else:
        L[index,:] = getExternalTorque(B[index,3:], B[index-1,3:], EP[index,:], L[index,:], index)

    # Get the change in body angular velocity (omegaDot) and the change in the Euler parameters (EP).
    omegaDot[index,:] = getOmegaDot(omega[index,:], L[index,:], index)
    EPdot[index,:] = getEPdot(EP[index,:], omega[index,:], index)
    print(EPdot[index,:])
    
    ## Euler's Method to propagate EP and omega
    EP[index+1,:] = EP[index,:] + dt*EPdot[index,:]
    omega[index+1,:] = omega[index,:] + dt*omegaDot[index,:]

    # Get the 3-2-1 Euler angle sequence from the current EP for later plotting
    print(EP[index,:])
    theta = DCMtoEulerAngles(getDCM(EP[index,:]))

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
# Plot of the agular velocity and acceleration of the spacecraft in the body frame VS time.
plt.subplot(2, 1, 1)
plt.plot(tVec, omega[:,0])
plt.plot(tVec, omega[:,1])
plt.plot(tVec, omega[:,2])
plt.title('Angular Vel. & Acc. in The Body Frame VS Time')
plt.legend(['z-axis','y-axis','x-axis'])
plt.ylabel('Ang. Vel. Body (deg/sec)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)
# Plot of change in EP's VS time.
plt.subplot(2, 1, 2)
plt.plot(tVec, omegaDot[:,0])
plt.plot(tVec, omegaDot[:,1])
plt.plot(tVec, omegaDot[:,2])
plt.legend(['x-axis','y-axis','z-axis'])
plt.xlabel('Time (seconds)')
plt.ylabel('Ang. Acc. Body (deg/sec^2)')
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
   


