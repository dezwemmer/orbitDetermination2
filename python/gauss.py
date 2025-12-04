### Gauss Angles-Only Algorithm
# Description:
# This function is designed to take in 3 LOS unit vectors, 
# 3 Julian dates, and 3 site vectors and use Gauss-Angles only
# as defined by Vallado to calculate a state estimate (r,v) of
# the second measurement.

import constants as c
from gibbs import gibbsAlgo
from math import acos, degrees
import numpy as np
from rv2OrbElem import rv2OrbElem

def gaussAlgo(obj,L,jd,rSite_Eci):
    # calculate taus (change in times between observations)
    tau1 = (jd[0] - jd[1])*24*60*60
    tau3 = (jd[2] - jd[1])*24*60*60

    # Calculate parameters for estimating the middle range magnitude
    a1 = tau3 / (tau3 - tau1)
    a1u = tau3*((tau3 - tau1)**2 - tau3**2) / (6*(tau3 - tau1))
    a3 = - tau1 / (tau3 - tau1)
    a3u = - tau1*((tau3 - tau1)**2 - tau1**2) / (6*(tau3 - tau1))

    # Determinant of L matrix
    Ldet = np.linalg.det(L)

    # Calculate the Inverse of the L array (LOS vectors) and transpose for multiplication
    Linv = np.transpose(np.array([[  L[1][1]*L[2][2]-L[1][2]*L[2][1], -L[1][0]*L[2][2]+L[1][2]*L[2][0],  L[1][0]*L[2][1]-L[1][1]*L[2][0] ],
                        [ -L[0][1]*L[2][2]+L[0][2]*L[2][1],  L[0][0]*L[2][2]-L[0][2]*L[2][0], -L[0][0]*L[2][1]+L[0][1]*L[2][0] ],
                        [  L[0][1]*L[1][2]-L[0][2]*L[1][1], -L[0][0]*L[1][2]+L[0][2]*L[1][0],  L[0][0]*L[1][1]-L[0][1]*L[1][0] ]]) / Ldet)

    # Intermediate relationship and parameters for eighth-degree equation
    M = np.matmul(Linv,np.transpose(rSite_Eci))
    d1 = M[1][0]*a1 - M[1][1] + M[1][2]*a3
    d2 = M[1][0]*a1u + M[1][2]*a3u
    C = np.dot(L[:,1],rSite_Eci[1])
    coeffs8 = [1,0,(-(d1**2 + 2*C*d1 + np.linalg.norm(rSite_Eci[1])**2)),0,0,(-2*c.muE*(C*d2 + d1*d2)),0,0,(-(c.muE**2)*(d2**2))]
    roots8 = np.roots(coeffs8)
    realIdx = np.isreal(roots8)
    realRoots = roots8[realIdx].real

    # Use first root find an initial estimate of slant-range values
    u = c.muE / realRoots[0]**3
    cCoeffs = np.array([a1 + a1u * u, -1, a3 + a3u * u])
    #TODO: make slant range finder a function
    slantsInitial = np.divide(np.matmul(M,-cCoeffs),cCoeffs)

    print("Slant Ranges Initial Guess: :", slantsInitial)

    ## End of Gauss

    
    # Note:
    #   Herrick-Gibbs is used when measurements are close together (<=5 degrees)
    #   Gibbs is used when measurements are far apart (5-20 degrees)
    #   Larger than 20 degrees is too large of measurement for reliable accuracy
    
    # Use slant ranges to get position vectors
    r1 = slantsInitial[0]*L[:,0]+rSite_Eci[0]
    r2 = slantsInitial[1]*L[:,1]+rSite_Eci[1]
    r3 = slantsInitial[2]*L[:,2]+rSite_Eci[2]
    obj.r = r2

    # Determine angles between vectors. The size of angles determines the method
    # used for finding the velocity of the middle measurement.
    angle12 = degrees(acos(np.dot(r1,r2)/(np.linalg.norm(r1)*np.linalg.norm(r2))))
    angle23 = degrees(acos(np.dot(r2,r3)/(np.linalg.norm(r2)*np.linalg.norm(r3))))
    
    
    # Find middle velocity (v2) using either Gibbs of Herrick-Gibbs
    if (5 < angle12 <=20) or (5 < angle23 <=20):
        print("Measurement separation large...Using GIBBS method to solve for middle velocity.")
        obj.v = gibbsAlgo(r1,r2,r3)
    else:
        print("Measurement separation small...Using HERRICK-GIBBS method to solve for middle velocity.")
    
    rv2OrbElem(obj)
    # Find semiparameter (p) using RV2COE
    # TODO: make this a separate function
    
