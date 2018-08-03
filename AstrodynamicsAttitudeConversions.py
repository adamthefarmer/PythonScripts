

## Import libraries.
import math
import numpy as np
import matplotlib as mpl

# Sine and cosine shorthand functions 
def s(angle):
    result = np.sin(np.deg2rad(angle))
    return result

def c(angle):
    result = np.cos(np.deg2rad(angle))
    return result

## Conversion functions:

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

def EPtoDCM(EP):

    DCM = np.array([[EP[0]**2 + EP[1]**2 - EP[2]**2 - EP[3]**2, 2*(EP[1]*EP[2] + EP[0]*EP[3]), 2*(EP[1]*EP[3] - EP[0]*EP[2])],
                           [2*(EP[1]*EP[2] - EP[0]*EP[3]), EP[0]**2 - EP[1]**2 + EP[2]**2 - EP[3]**2, 2*(EP[2]*EP[3] + EP[0]*EP[1])],
                           [2*(EP[1]*EP[3] + EP[0]*EP[2]), 2*(EP[2]*EP[3] - EP[0]*EP[1]), EP[0]**2 - EP[1]**2 - EP[2]**2 + EP[3]**2]])


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

    #Returns the Euler Angles in radians
    EA[0] = np.arctan(DCM[2,0]/(-DCM[2,1]))
    EA[1] = np.arccos(DCM[2,2])
    EA[2] = np.arctan(DCM[0,2]/DCM[1,2])

    # Convert the Euler Angles into degrees
    EA = EA*(180/np.pi)

    return EA
