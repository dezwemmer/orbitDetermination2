### RV2OrbElem method
# Description:
# This function is calculate the classical orbital elements from
# a position and velocity vector in the ECI frame (IJK)
#
# Inputs: position and velocity vector 

import constants as c
from math import degrees,acos
import numpy as np

def rv2OrbElem(obj):
    obj.r = np.array([9031.5,-5316.9,-1647.2]) #test
    obj.v = np.array([-2.8640,5.1112,-5.0805]) #test
    r = np.linalg.norm(obj.r)
    v = np.linalg.norm(obj.v)

    # Calculate total energy, E
    E = (v**2)/2 - c.muE/r

    obj.a = -c.muE/(2*E)
    eVec = (1/c.muE)*((v**2 - c.muE/r)*obj.r - np.dot(obj.r,obj.v)*obj.v)
    obj.e = np.linalg.norm(eVec)

    # To calculate inclination, we need angular momentum vector
    h = np.cross(obj.r,obj.v)

    K = np.array([0,0,1])
    cosi = np.dot(K,h)/np.linalg.norm(h)
    print ("i",degrees(acos(cosi)))