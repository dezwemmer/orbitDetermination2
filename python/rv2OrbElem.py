### RV2OrbElem method
# Description:
# This function is calculate the classical orbital elements from
# a position and velocity vector in the ECI frame (IJK)
#
# Inputs: SatState object with r and v defined
# Outputs: This will set the object's a,e,i,lan,ap

import constants as c
from math import degrees,acos
import numpy as np

def rv2OrbElem(obj):
    # obj.r = np.array([9031.5,-5316.9,-1647.2]) #test
    # obj.v = np.array([-2.8640,5.1112,-5.0805]) #test
    r = np.linalg.norm(obj.r)
    v = np.linalg.norm(obj.v)

    # Calculate total energy, E
    E = (v**2)/2 - c.muE/r

    # Calculate semimajor axis (a) and eccentricity (e)
    obj.a = -c.muE/(2*E)
    eVec = (1/c.muE)*((v**2 - c.muE/r)*obj.r - np.dot(obj.r,obj.v)*obj.v)
    obj.e = np.linalg.norm(eVec)

    # To calculate inclination, we need angular momentum vector
    h = np.cross(obj.r,obj.v)

    # ECI unit vectors
    I = np.array([1,0,0])
    J = np.array([0,1,0])
    K = np.array([0,0,1])
    
    # Calculate inclination (i)
    cosi = np.dot(K,h)/np.linalg.norm(h)
    obj.i = degrees(acos(cosi))
    print ("i: ",degrees(acos(cosi)))

    # Calculate Longitude of Ascending Node (lan)
    # calculate ascending node vector (n)
    nVec = np.cross(K,h)
    
    coslan = np.dot(I,nVec)/np.linalg.norm(nVec)
    if nVec[1] < 0: # quadrant check based on node vector
        obj.lan = degrees(2*np.pi - acos(coslan))
    else:
        obj.lan = degrees(acos(coslan))
    print("LAN: ", obj.lan)

    # Calculate Argument of Periapsis (ap)
    # Arg of Periapsis is the angle between node and eccentricity vectors (eVec)
    cosap = np.dot(nVec,eVec)/(np.linalg.norm(nVec)*obj.e)
    if eVec[2] < 0: # quadrant check based on eccentricity vector
        obj.ap = degrees(2*np.pi - acos(cosap))
    else:
        obj.ap = degrees(acos(cosap))
    print("Arg of Periapsis: ",obj.ap)

    # Calculate True Anomaly (ta)
    costa = np.dot(eVec,obj.r)/(obj.e*r)
    if np.dot(obj.r,obj.v):
        obj.ta = degrees(2*np.pi - acos(costa))
    else:
        obj.ta = degrees(acos(costa))
    print("True Anomaly: ",obj.ta)
