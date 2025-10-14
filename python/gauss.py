### Implementation of Angles-only Guass method using 3 JD/RA/Dec measurements

from math import cos,sin,radians, degrees, sqrt
import numpy as np

### Functions
def calcLosUnitVector(ra,dec):
    r = radians(ra)
    d = radians(dec)
    L = [cos(d) * cos(r), \
         cos(d) * sin(r), \
         sin(d)]
    return L

def gmst(jd):
    # Vallado page 188 (ex 3.5)
    Tuti = (jd - 2451545.0) / 36525.0
    gmst_sec = 67310.54841 + (876600*3600 + 8640184.812866)*Tuti + 0.093104*(Tuti**2) - (6.2e-6)*(Tuti**3)

    gmst_sec = gmst_sec % 86400.0
    gmst_deg = gmst_sec *(360.0 / 86400.0)
    return radians(gmst_deg)

def ecef2eci(ecef, gmstrad):
    # Passive rotation (rotates the coordinate system)
    rot = np.array([[cos(gmstrad), sin(gmstrad), 0],
                   [-sin(gmstrad), cos(gmstrad), 0],
                   [0, 0, 1]])
    eci_coords = np.dot(rot,ecef)
    return eci_coords   

###

### Constants   
rE = 6378.1366          # Semi-major axis of Earth [km]
eE = 0.0818191908426215
omegaE = 72.92115E-06
muE = 3.986e5           # gravitational constant [km^3/s^2]

### Inputs (replace this with a JSON or similiar input)
# Three different observations with a JD, Right Ascension (deg), and Declination (deg)
latSite = 40    # (deg)
lonSite = -110  # (deg) west is negative
altSite = 2000  # (m)
m1 = [2456159.986435, 0.939913, 18.667717]
m2 = [2456159.991991, 45.025748, 35.664741]
m3 = [2456159.994769, 67.886655, 36.996583]
m = [m1,m2,m3]
jd1 = 2456159.986435
jd2 = 2456159.991991
jd3 = 2456159.994769

### Calculate Site ECEF coordinates
# Calculate auxiliary terms
auxC = rE / sqrt(1 - (eE**2)*sin(radians(latSite))**2)
auxS = auxC * (1 - eE**2)

# Calculate vertical and horizontal components
rDelta = (auxC + altSite/1000) * cos(radians(latSite))
rKappa = (auxS + altSite/1000) * sin(radians(latSite))

# Calculate site position in ECEF
rSite_Ecef = [rDelta * cos(radians(lonSite)), rDelta * sin(radians(lonSite)), rKappa]

# To convert site coordinates from ECEF > ECI, we need the GMST at each JD
gmst1 = gmst(jd1)
gmst2 = gmst(jd2)
gmst3 = gmst(jd3)

# testgmst = gmst(2448855.009722)
# print(degrees(testgmst))

# Calculate site position in ECI at each epoch
# Note: We apply a negative angle here because we are doing a passive rotation, whereby
# we need to "undo" the rotation of the Earth.
rSite1_Eci = ecef2eci(np.array(rSite_Ecef),-gmst1)
rSite2_Eci = ecef2eci(np.array(rSite_Ecef),-gmst2)
rSite3_Eci = ecef2eci(np.array(rSite_Ecef),-gmst3)
# print("site eci1:",rSite1_Eci)
# print("site eci2:",rSite2_Eci)
# print("site eci3:",rSite3_Eci)


### Angles Only Gauss
# calculate taus (change in times between observations)
tau1 = (jd1 - jd2)*24*60*60
tau3 = (jd3 - jd2)*24*60*60

# Calculate parameters for estimating the middle range magnitude
a1 = tau3 / (tau3 - tau1)
a1u = tau3*((tau3 - tau1)**2 - tau3**2) / (6*(tau3 - tau1))
a3 = - tau1 / (tau3 - tau1)
a3u = - tau1*((tau3 - tau1)**2 - tau1**2) / (6*(tau3 - tau1))

# LOS vectors for 3 observation measurements
L = []
for i in range(3):
    los = calcLosUnitVector(m[i][1], m[i][2])
    # print(los)
    L.append(los)
L = np.transpose(L)

# Determinant of L matrix
Ldet = np.linalg.det(L)

# Not entirely sure why this needs to be transposed
Linv = np.transpose(np.array([[  L[1][1]*L[2][2]-L[1][2]*L[2][1], -L[1][0]*L[2][2]+L[1][2]*L[2][0],  L[1][0]*L[2][1]-L[1][1]*L[2][0] ],
                    [ -L[0][1]*L[2][2]+L[0][2]*L[2][1],  L[0][0]*L[2][2]-L[0][2]*L[2][0], -L[0][0]*L[2][1]+L[0][1]*L[2][0] ],
                    [  L[0][1]*L[1][2]-L[0][2]*L[1][1], -L[0][0]*L[1][2]+L[0][2]*L[1][0],  L[0][0]*L[1][1]-L[0][1]*L[1][0] ]]) / Ldet)

# Compile the 3 ECI site vectors into a single matrix
rSite_Eci = np.transpose(np.array([rSite1_Eci, rSite2_Eci, rSite3_Eci]))

# Intermediate relationship and parameters for eighth-degree equation
M = np.matmul(Linv,rSite_Eci)
d1 = M[1][0]*a1 - M[1][1] + M[1][2]*a3
d2 = M[1][0]*a1u + M[1][2]*a3u
C = np.dot(L[:,1],rSite2_Eci)
coeffs8 = [1,0,(-(d1**2 + 2*C*d1 + np.linalg.norm(rSite2_Eci)**2)),0,0,(-2*muE*(C*d2 + d1*d2)),0,0,(-(muE**2)*(d2**2))]
roots8 = np.roots(coeffs8)
realIdx = np.isreal(roots8)
realRoots = roots8[realIdx].real

# Use first root find an initial estimate of slant-range values
u = muE / realRoots[0]**3
cCoeffs = np.array([a1 + a1u * u, -1, a3 + a3u * u])
slantsInitial = np.divide(np.matmul(M,-cCoeffs),cCoeffs)

### [End of Gauss Angles Only]
print(L)
print(rSite1_Eci)
print(slantsInitial)

L = np.transpose(L)
print('r1 = ',slantsInitial[0]*L[0] + rSite1_Eci)

# TODO: fix L and other array transposes
# TODO: add iteration on slant range