### Gauss Angles-Only Algorithm
# Description:
# This function is designed to take in 3 LOS unit vectors, 
# 3 Julian dates, and 3 site vectors and use Gauss-Angles only
# as defined by Vallado to calculate a state estimate (r,v) of
# the second measurement.

import constants as c
import numpy as np

def gaussAlgo(L,jd,rSite_Eci):
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

    # Not entirely sure why this needs to be transposed
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

    print("Debug",rSite_Eci[1])

    # Use first root find an initial estimate of slant-range values
    u = c.muE / realRoots[0]**3
    cCoeffs = np.array([a1 + a1u * u, -1, a3 + a3u * u])
    slantsInitial = np.divide(np.matmul(M,-cCoeffs),cCoeffs)

    print("Slant Range Initial Guess: :", slantsInitial)