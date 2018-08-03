############################################################################################
#
# Title: 3-2-1 Euler Sequence to Quaternion
#
# Purpose: Take a 3-2-1 Euler sequence and convert it into quaternion
#
# Author: Adam Farmer
# Sources: Analytical Mechanics of Space Systems - 3rd Edition, pg.82
#          http://mathworld.wolfram.com/EulerAngles.html
#
# Date Written:  7/16/18
# Date Modified: 7/16/18
#
############################################################################################

## Import libraries.

import math
import numpy as np
import matplotlib as mpl

# Sine and cosine functions 
def s(angle):
    result = np.sin(np.deg2rad(angle))
    return result

def c(angle):
    result = np.cos(np.deg2rad(angle))
    return result

## Definition of the function which does the conversion
def Euler321toDCM(EA):

    DCM = np.array([[c(EA[1])*c(EA[0]), c(EA[1])*s(EA[0]), -s(EA[1])],
                    [s(EA[2])*s(EA[1])*c(EA[0]) - c(EA[2])*s(EA[0]), s(EA[2])*s(EA[1])*s(EA[0]) + c(EA[2])*c(EA[0]), s(EA[2])*c(EA[1])],
                    [c(EA[2])*s(EA[1])*c(EA[0]) + s(EA[2])*s(EA[0]), c(EA[2])*s(EA[1])*s(EA[0]) - s(EA[2])*c(EA[0]), c(EA[2])*c(EA[1])]])
    
    return DCM

def Euler313toDCM(EA):

    DCM = np.array([[c(EA[2])*c(EA[0]) - s(EA[2])*c(EA[1])*s(EA[0]), c(EA[2])*s(EA[0]) + s(EA[2])*c(EA[1])*c(EA[0]), s(EA[2])*s(EA[1])],
                    [-s(EA[2])*c(EA[0]) - c(EA[2])*c(EA[1])*s(EA[0]), -s(EA[2])*s(EA[0]) + c(EA[2])*c(EA[1])*c(EA[0]), c(EA[2])*s(EA[1])],
                    [s(EA[1])*s(EA[0]), -s(EA[1])*c(EA[0]), c(EA[1])]])

    return DCM

def DCMtoEP(DCM):
    EP = np.array([0.0, 0.0, 0.0, 0.0])

    EP[0] = 0.5*np.sqrt(DCM[0,0] + DCM[1,1] + DCM[2,2] + 1)
    EP[1] = (DCM[1,2] - DCM[2,1])/(4*EP[0])
    EP[2] = (DCM[2,0] - DCM[0,2])/(4*EP[0])
    EP[3] = (DCM[0,1] - DCM[1,0])/(4*EP[0])

    return EP

def DCMto321Euler(DCM):
    EA = np.array([0.0, 0.0, 0.0])

    #Returns the Euler Angles in radians
    EA[0] = np.arctan(DCM[0,1]/DCM[0,0])
    EA[1] = -np.arcsin(DCM[0,2])
    EA[2] = np.arctan(DCM[1,2]/DCM[2,2])

    # Convert the Euler Angles into degrees
    EA = EA*(180/np.pi)

    return EA

def DCMto313Euler(DCM):
    EA = np.array([0.0, 0.0, 0.0])
    
    EA[0] = np.arctan(DCM[2,0]/(-DCM[2,1]))
    EA[1] = np.arccos(DCM[2,2])
    EA[2] = np.arctan(DCM[0,2]/DCM[1,2])

    # Convert the Euler Angles into degrees
    EA = EA*(180/np.pi)

    return EA

# Set the Euler Angle to be converted
EA313 = np.array([9.0, 4.0, 0.0])

DCM = Euler313toDCM(EA313)
EP = DCMtoEP(DCM)
eulerAgain = DCMto313Euler(DCM)

print("This is the quaternion: ", EP)
print("These are the Euler angles: ", eulerAgain)

thing1 = [1, 1, 1]
thing2 = [2, 4, 6]

print("this is the result: ", np.multiply(thing1,thing2))

