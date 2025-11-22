### Gibbs method
# Description:
# This function is used to determine velocity of middle measurement
# when angle between measurements is large (5-20 degrees)
#
# Inputs: 3 range vectors in ECI frame

import constants as c
from math import asin,sqrt
import numpy as np

def gibbsAlgo(rVec1, rVec2, rVec3):
    z12 = np.cross(rVec1,rVec2)
    z23 = np.cross(rVec2,rVec3)
    z31 = np.cross(rVec3,rVec1)

    alphaCop = asin(np.dot(z23,rVec1)/(np.linalg.norm(z23)*np.linalg.norm(rVec1)))

    r1 = np.linalg.norm(rVec1)
    r2 = np.linalg.norm(rVec2)
    r3 = np.linalg.norm(rVec3)

    N = r1*z23 + r2*z31 + r3*z12
    D = z12 + z23 + z31
    S = (r2-r3)*rVec1 + (r3-r1)*rVec2 + (r1-r2)*rVec3

    B = np.cross(D,rVec2)

    Lg = sqrt(c.muE/(np.linalg.norm(N)*np.linalg.norm(D)))

    return (Lg/r2)*B + Lg*S