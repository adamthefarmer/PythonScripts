############################################################################################
#
# Title: Passive ACS Simulation
#
# Purpose: The purpose of this simulation is to test the effectiveness of PolarCube's
#          passive ACS system. Given a moment of inertia tensor, a bar-magnet strength and
#          knowledge about the hysteresis rods being used on the satellite, this simulation
#          will be able to tell the user the settling time 
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
import csv

#### Define some constants used in this program.

## Orbit propagation
mEarth = 5.972 * (10**24) # kg
rEarth = 6.371 * (10**6) # meters
G_Earth = 6.67408 * (10**(-11)) # m^3 kg^-1 s^-2

# Earth's magnetic field in LEO varies between 0.2 and 0.45 Gauss
# We can choose that B = approx 3.2 * 10^ -5 Tesla
Mz = 8*(10**15) #T*m^3

## Choose test case: 1=PolarCbue, 2=CSSWE
scenario = 1

if scenario == 1: # PolarCube parameters
    ## Hysteresis rod stuff
    mu_0 = 4*np.pi*(10**(-7))
    mu_hyst = 1.5*(10**4)
    lengthHystXY = 0.0762 # Length of the x and y axis hysteresis rod, meters
    lengthHystZ = 0.1524 # length of the z-axis hysteresis rods, meters
    Diameter_hyst = 0.000762 # Diameter of hysteresis rod, meters
    V_hyst = [np.pi*((Diameter_hyst/2)**2)*lengthHystXY, np.pi*((Diameter_hyst/2)**2)*lengthHystXY, np.pi*((Diameter_hyst/2)**2)*lengthHystZ] # Volume of the hysteresis rod, m^3
    m_hyst = np.zeros(3)
    Br = 0.004 # Tesla
    Bs = 0.025 # Tesla
    Hc = 12.0 # A/m
    numHyst = [16.0, 16.0, 6.0] # Define the number of hysteresis rods on each axis.

    ## Position propagation stuff:
    # Initial position of the spacecraft.
    ri = [0.0, 0.0, 500000+rEarth] # Meters
    # This is the magnitude of the needed velocity to stay in a stable LEO.
    V_magLEO = np.sqrt((G_Earth*mEarth)/(np.sqrt(ri[0]**2 + ri[1]**2 +ri[2]**2)))
    # Initial velocity of the spacecraft.
    vi = [V_magLEO, 0.0, 0.0]

    ## Attitude propagation stuff:
    # Moments of inertia about the spacecraft's principle axes.
    I = [0.112209, 0.124612, 0.016524]
    # Maximum bipole moment of the on-board magnet. Units are A*m^2
    mBarMax = np.array([0.0, 0.0, 0.3])
    EP_i = np.array([np.cos(np.deg2rad(45/2)), np.sin(np.deg2rad(45/2)), 0.0, 0.0])  # Initial Euler parameters.
    omega_i = np.radians(np.array([0.0, 0.0, 0.0]))  # Initial body fixed angular velocity in degrees/sec
    
elif scenario == 2: # CSSWE parameters
    ## Hysteresis rod stuff
    mu_0 = 4*np.pi*(10**(-7))
    mu_hyst = 1.5*(10**4)
    Length_hyst = 0.095 # Length of the hysteresis rod, meters
    Diameter_hyst = 0.001 # Diameter of hysteresis rod, meters
    V_hyst = np.pi*((Diameter_hyst/2)**2)*Length_hyst # Volume of the hysteresis rod, m^3
    m_hyst = np.zeros(3)
    Br = 0.004 # Tesla
    Bs = 0.025 # Tesla
    Hc = 12.0 # A/m
    numHyst = [2.0, 2.0, 0.0] # Define the number of hysteresis rods on each axis.

    ## Position propagation stuff:
    # Initial position of the spacecraft.
    ri = [600000+rEarth, 0.0, 0.0] # Meters
    # This is the magnitude of the needed velocity to stay in a stable LEO.
    V_magLEO = np.sqrt((G_Earth*mEarth)/(np.sqrt(ri[0]**2 + ri[1]**2 +ri[2]**2)))
    # Initial velocity of the spacecraft.
    vi = [0.0, np.cos(np.deg2rad(55))*V_magLEO, np.sin(np.deg2rad(55))*V_magLEO]

    ## Attitude propagation stuff:
    # Moments of inertia about the spacecraft's principle axes.
    I = [0.00551, 0.02552, 0.02565]
    # Bipole moment of the on-board magnet. Units are A*m^2
    mBar = [0.0, 0.0, 0.29]
    EP_i = np.array([1.0, 0.0, 0.0, 0.0])
    omega_i = np.radians(np.array([10.0, 5.0, 5.0]))  # Initial body fixed angular velocity in degrees/sec
    
else:
    print("That is not a valid scenario brah..")
    quit()

#### Initial state vectors, time and whatnot

## Time stuff
# Find the orbital period
orbitalPeriod = 2*np.pi*np.linalg.norm(ri)/V_magLEO
t0 = 0.0 # Seconds
dt = .1 # Seconds
tf =  5.5*86400 #1*orbitalPeriod # Seconds
tVec = np.linspace(t0, tf-dt, ((tf-t0)/dt)) # Make time vector.

# Initialize external torques, and hysteresis arrays.
L = np.zeros((len(tVec),6))
B_hyst = np.array([0.0, 0.0, 0.0])
m_hyst = np.array([0.0, 0.0, 0.0])

# Initialize some lists to hold angular rate and position values for the system.
global state
global deltaState
state = np.zeros((len(tVec),25))
deltaState = np.zeros((len(tVec),13))

#### Definitions of functions used in this program
 
def s(angle):
    result = np.sin(np.deg2rad(angle))
    return result

def c(angle):
    result = np.cos(np.deg2rad(angle))
    return result

def normalizeAngle(angle):
    # Take some angle in degrees and normalize it to a value between 0 and 360.
    normalizedAngle = angle % 360.0
    #normalizedAngle = angle % (2*np.pi)

    return normalizedAngle

def getCurrentState():
    # Gets the current state of the system. This will be a vector in the format:
    # state = [r_ECI_x, r_ECI_y, r_ECI_z,               index: 0-2
    #          v_ECI_x, v_ECI_y, v_ECI_z,               index: 3-5
    #          omega_x, omega_y, omega_z                index: 6-8
    #          EP_0, EP_1, EP_2, EP_3,                  index: 9-12 
    #          B_ECIx, B_ECIy, B_ECIz                   index: 13-15
    #          B_body_x, B_body_y, B_body_z             index: 16-18
    #          L_1, L_2, L_3, L_4, L_5, L_6]            index: 19-24

    # Orbit stuff
    if index is 0:
        state[index,0:3] = ri                                                              # Position
        state[index,3:6] = vi                                                              # Velocity
        state[index,9:13] = EP_i                                                           # Quaternion
        state[index,6:9] = omega_i                                                         # Omega
        state[index,13:19] = getBvec(state[index,0:3], state[index,9:13], index)           # B
        state[index,19:] = getExternalTorque(state[index,16:19], state[index,16:19])       # L
    else:
        state[index,9:13] =  state[index,9:13] / np.linalg.norm(state[index,9:13])
        state[index,13:19] = getBvec(state[index,0:3], state[index,9:13], index)           # B
        state[index,19:] = getExternalTorque(state[index,16:19], state[index-1,16:19])     # L

def getChangeInState(RKstate, ti):
    # Returns the change in the state of the system for this timestep. This will be a vector in the format:
    # deltaState = [v_ECI_x, v_ECI_y, v_ECI_z,              index: 0-2
    #               a_ECI_x, a_ECI_y, a_ECI_z,              index: 3-5
    #               omegaDot_x, omegaDot_y, omegaDot_z,     index: 6-8
    #               EPdot_0, EPdot_1, EPdot_2, EPdot_3]     index: 9-12 -AF

    # Orbit stuff at current timestep.
    deltaState[index,0:3] = RKstate[3:6]                                                    # Velocity
    deltaState[index,3:6] = geta_ECI(RKstate[0:3])                                          # Accelteration

    # Attitude stuff at current timestep.
    deltaState[index,6:9] = getOmegaDot(RKstate[6:9], state[index,19:], index)              # OmegaDot
    deltaState[index,9:] = getEPdot(RKstate[9:13], RKstate[6:9], index)                     # EPdot
    
    return deltaState[index,:]

def RungeKuttaIntegrator(currentState, t_current, dydt):
    # A 4th order Runge-Kutta integrator.

    # Define Runge-Kutta coefficients
    k1 = dt*dydt(currentState[0:13], t_current)
    k2 = dt*dydt(currentState[0:13] + k1*(1/2), t_current + dt/2)
    k3 = dt*dydt(currentState[0:13] + k2*(1/2), t_current + dt/2)
    k4 = dt*dydt(currentState[0:13] + k3*1, t_current + dt)

    # Calculate and return new value of y
    newState = currentState[0:13] + ((k1 + 2*k2 + 2*k3 + k4)/6)
    
    return newState

def getDCM(EP):

    DCMcurrent = np.array([[EP[0]**2 + EP[1]**2 - EP[2]**2 - EP[3]**2, 2*(EP[1]*EP[2] + EP[0]*EP[3]), 2*(EP[1]*EP[3] - EP[0]*EP[2])],
                           [2*(EP[1]*EP[2] - EP[0]*EP[3]), EP[0]**2 - EP[1]**2 + EP[2]**2 - EP[3]**2, 2*(EP[2]*EP[3] + EP[0]*EP[1])],
                           [2*(EP[1]*EP[3] + EP[0]*EP[2]), 2*(EP[2]*EP[3] - EP[0]*EP[1]), EP[0]**2 - EP[1]**2 - EP[2]**2 + EP[3]**2]])

    return DCMcurrent

def isRotationMatrix(DCM):
    # Checks if a matrix is a valid rotation matrix.
    DCMt = np.transpose(DCM)
    shouldBeIdentity = np.dot(DCMt, DCM)
    I = np.identity(3, dtype = DCM.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6
 
def DCMtoEulerAngles(DCM):
    # Calculates rotation matrix to euler angles
 
    assert(isRotationMatrix(DCM))
     
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

def getMbar(mBarMax):
    ## Control law for PolarCube's magnetic torque rod. If the angle between the body z-axis and the
    ## Earth's magnetic field vector is increasing then we turn the rod on full throttle. If the
    ## angle is decreasing but not fast enough, then we turn on the rod full throttle.
    
    # Get the angle between the z-axis and Earth's magnetic field vector in the body frame for this and the last timestep.
    angCurrent = np.rad2deg(np.arccos(np.dot(state[index,16:19],[0,0,1]) / (np.linalg.norm(state[index,16:19]))))
    angLast = np.rad2deg(np.arccos(np.dot(state[index-1,16:19],[0,0,1]) / (np.linalg.norm(state[index-1,16:19]))))

    print(angLast - angCurrent)

    if (angLast-angCurrent) < 3*dt: # If the angle is decreasing but not fast enough turn m_bar on.
        mBar = mBarMax
    else:
        mBar = [0.0, 0.0, 0.0]
        print("this happened")

    return mBar
    
    
def getExternalTorque(B_current, B_last):
    global B_hyst
    global m_hyst
    # Get the external torques on the spacecraft based on current attitude. Currently the external torques on the spacecraft
    # are the torque produced by the interaction of the on-board bar magnet with the Earth's magnetic field, and the damping
    # torque produced by the changing magnetic field inducing currents in the hysteresis rods.
    L_current = np.zeros(6)
    
    # Since the magnetic diple moment is always aligned with the z-axis this will never change.
    # The torque produced by the bar magnet in Earth's B field is then given by: tao = m cross B
    if scenario == 1:
        if index == 0:
            mBar = mBarMax
        else:
            mBar = getMbar(mBarMax)
    elif scenario == 2:
        pass
    else:
        print("How did you even get this far? That is truly troubling...")
        quit()

    L_current[0:3] = np.cross(mBar,B_current)
    
    # The damping torque produced by the hysteresis rods is described on the bottom of pg.3 in the paper cited.
    p = (1/Hc)*np.tan((np.pi*Br)/(2*Bs)) # Function of material properties, so constant for all rods.

    sign =[-1.0 if c - l > 0 else 1.0 for c, l in zip(B_current, B_last)]
    
    B_hyst[0] = (2/np.pi)*Bs*np.arctan(p*(B_current[0]/mu_0 + sign[0]*Hc))
    B_hyst[1] = (2/np.pi)*Bs*np.arctan(p*(B_current[1]/mu_0 + sign[1]*Hc))
    B_hyst[2] = (2/np.pi)*Bs*np.arctan(p*(B_current[2]/mu_0 + sign[2]*Hc))

    # Hysteresis rod magnetic moment components.
    m_hyst = numHyst * np.multiply(B_hyst,V_hyst) / mu_0
    
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

    omegaDot = np.array([1/I[0]*(-(I[2]-I[1])*omega[1]*omega[2] + L[0] + L[3]),
                         1/I[1]*(-(I[0]-I[2])*omega[2]*omega[0] + L[1] + L[4]),
                         1/I[2]*(-(I[1]-I[0])*omega[0]*omega[1] + L[2] + L[5])])

    # Test case with no torques:
##    omegaDot = np.array([1/I[0]*(-(I[2]-I[1])*omega[1]*omega[2]),
##                         1/I[1]*(-(I[0]-I[2])*omega[2]*omega[0]),
##                         1/I[2]*(-(I[1]-I[0])*omega[0]*omega[1])])

    # Test case with artificial damping:
##    omegaDot = np.array([1/I[0]*(-(I[2]-I[1])*omega[1]*omega[2] + L[0] - omega[0]*.01),
##                         1/I[1]*(-(I[0]-I[2])*omega[2]*omega[0] + L[1] - omega[0]*.01),
##                         1/I[2]*(-(I[1]-I[0])*omega[0]*omega[1] + L[2] - omega[0]*.01)])

    return omegaDot

def getUnitVec(vec):
    # returns the unit vector of the input vector.
    unitVec = vec / np.linalg.norm(vec)

    return unitVec

def geta_ECI(r_ECI):
    # Using Newton's Law of Gravitation, find the acceleration vector of the spacecraft at the current timestep.
    
    # Get the magnitude of the acceleration vector
    a_ECImag = (G_Earth*mEarth)/(np.linalg.norm(r_ECI)**2)
    # Get the unit vector of the acceleration by noting that it is always antiparallel to the position vector.
    a_ECIhat = -getUnitVec(r_ECI)
    # Construct the final acceleration vector.
    a_ECI = a_ECImag*a_ECIhat

    return a_ECI

def getBvec(r_ECI, EP, index):
    # Given some Cartesian coordinate vector in the ECI frame, this function will return the
    # approximate magnetic field vector according to the ECI frame and also the body frame.
    B = np.zeros(6)

    # Good approximation to Earth's magnetic field
    B[0] = -Mz*(3*r_ECI[0]*r_ECI[2]) / (((np.linalg.norm(r_ECI))**5))
    B[1] = -Mz*(3*r_ECI[1]*r_ECI[2]) / (((np.linalg.norm(r_ECI))**5))
    B[2] = Mz*(3*(r_ECI[2]**2) - (((np.linalg.norm(r_ECI))**2))) / (((np.linalg.norm(r_ECI))**5))

    # Trivial case of homogeneous B field.
##    B[0] = 0.0
##    B[1] = 0.2e-5
##    B[2] = 0.0

    ## Also use the current DCM to get the B vector in the body frame as well

    # Get the DCM for the current EP.
    DCM = getDCM(EP)
    # Rotate the magnetic field vector from the ECI frame into the body frame.
    B[3:] = np.matmul(DCM, B[:3])

    return B

#### Loop which calls functions to get the current state of the system and integrate to the next timestep

for index, ti in enumerate(tVec[:-1]):
    
    getCurrentState()

    state[index+1,0:13] = RungeKuttaIntegrator(state[index,:], ti, getChangeInState)

#### Writing the quaternions to a CSV file.
outBase = 'quaterionsOutput.csv'
outDir = '/home/adamthefarmer/Dropbox/PolarCube/ADCS/PassiveACSSsimulation/'

data = np.zeros((len(tVec),5))
data[:,0] = tVec
data[:,1:] = state[:,9:13]
np.savetxt(outDir+outBase, data, delimiter=",")

#### Plotting the results

# Turn on the minor TICKS, which are required for the minor GRID.
mpl.rcParams['legend.fontsize'] = 10

## Position plotting:
fig = plt.figure(1)
ax = fig.gca(aspect='equal',projection='3d')

ax.plot(state[:,0]/1000, state[:,1]/1000, state[:,2]/1000, label='Orbital Position', color='k')

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
ax.set_aspect('equal')

## Plotting the magnetic field
plt.figure(2)
plt.minorticks_on()
# Plot of B in the ECI frame VS time.
plt.subplot(2, 1, 1)
plt.plot(tVec, state[:,13:16])
plt.title('Magnetic Field in the ECI Frame')
plt.legend(['x-axis','y-axis','z-axis'])
plt.ylabel('B_ECI (Tesla)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)
# Plot of B in the body frame VS time.
plt.subplot(2, 1, 2)
plt.plot(tVec, state[:,16:19])
plt.legend(['x-axis','y-axis','z-axis'])
plt.xlabel('Time (seconds)')
plt.ylabel('B_body (Tesla)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)

## Plotting the external torques
plt.figure(3)
plt.minorticks_on()
# Plot of external torque from the bar magnet VS time.
plt.subplot(2, 1, 1)
plt.plot(tVec, state[:,19:22])
plt.title('External Torques on Spacecraft VS Time')
plt.legend(['x-axis','y-axis','z-axis'])
plt.ylabel('Bar Magnet Torque (N*m)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)
# Plot of B in the body frame VS time.
plt.subplot(2, 1, 2)
plt.plot(tVec, state[:,22:])
plt.legend(['x-axis','y-axis','z-axis'])
plt.xlabel('Time (seconds)')
plt.ylabel('Hysteresis Torque (N*m)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)

## Plotting of the quaternions and their change with respect to time.
plt.figure(4)
plt.minorticks_on()
# Plot of the quaternions VS time.
plt.subplot(2, 1, 1)
plt.plot(tVec, state[:,9])
plt.plot(tVec, state[:,10])
plt.plot(tVec, state[:,11])
plt.plot(tVec, state[:,12])
plt.title('Euler Parameters VS Time')
plt.legend(['EP_0','EP_1','EP_2','EP_3'])
plt.ylabel('Euler Parameters (dimensionless)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)
# Plot of change in EP's VS time.
plt.subplot(2, 1, 2)
plt.plot(tVec, deltaState[:,9])
plt.plot(tVec, deltaState[:,10])
plt.plot(tVec, deltaState[:,11])
plt.plot(tVec, deltaState[:,12])
plt.legend(['EPdot_0','EPdot_1','EPdot_2','EPdot_3'])
plt.xlabel('Time (seconds)')
plt.ylabel('EPdot (dimensionless)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)

## Plotting the angular velocity in the body frame.
plt.figure(5)
plt.minorticks_on()
# Plot of angular velocity VS time.
plt.subplot(2, 1, 1)
plt.plot(tVec[:], state[:,6])
plt.plot(tVec[:], state[:,7])
plt.plot(tVec[:], state[:,8])
plt.title('Angular Vel. in The Body Frame VS Time')
plt.legend(['x-axis','y-axis','z-axis'])
plt.ylabel('Ang. Vel. Body (deg/sec)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)
# Plot of change in EP's VS time.
plt.subplot(2, 1, 2)
plt.plot(tVec[:], deltaState[:,6])
plt.plot(tVec[:], deltaState[:,7])
plt.plot(tVec[:], deltaState[:,8])
plt.legend(['omegaDotx','omegaDoty','omegaDotz'])
plt.xlabel('Time (seconds)')
plt.ylabel('omegaDot (deg/s^2)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)

## Plot of rotational energy over time.
plt.figure(6)
plt.minorticks_on()
plt.plot(tVec[:], 0.5 * (I[0]*state[:,6]**2 + I[1]*state[:,7]**2 + I[2]*state[:,8]**2))
plt.plot(tVec[:], 0.5*I[0]*state[:,6]**2)
plt.plot(tVec[:], 0.5*I[1]*state[:,7]**2)
plt.plot(tVec[:], 0.5*I[2]*state[:,8]**2)
plt.title('Rotational Energy VS Time')
plt.legend(['Total','x-axis','y-axis','z-axis'])
plt.xlabel('Time (seconds)')
plt.ylabel('Rotational Energy (Joules)')
plt.grid(which='major', color='k', linestyle='-', linewidth=.1)
plt.grid(which='minor', color='k', linestyle=':', linewidth=.1)

plt.show()
