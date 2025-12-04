# Orbit Determination Function
# Description: Calls in the input JSON file with observation data
#              and performs common calculations necessary for 
#              orbit determination methods 

import constants as c
from gauss import gaussAlgo
import json
from math import cos,sin,radians, degrees, sqrt
import numpy as np


### Supporting Functions
class SatState:
    def __init__(self,r=0.0,v=0.0,a=0.0,e=0.0,i=0.0,lan=0.0,ap=0.0):
        self.r = r      # km
        self.v = v      # km/s
        self.a = a      # km
        self.e = e
        self.i = i      # degrees
        self.lan = lan  # degrees
        self.ap = ap    # degrees


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

### OD Main Function
def orbitDetermination(inputFileName):
    print("       Orbit Determination Script")

    ### Read in Input JSON file
    try:
        with open("inputs/" + inputFileName, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        # If file doesn't exist, end program gracefully
        print("Error: Input file '",inputFileName, "' not found.  Quitting.")
        quit()
    except json.JSONDecodeError:
        print("Error: Could not decode from specified JSON.  Quitting.")
        quit()

    ### Define terms from JSON
    latSite = data["siteLat"]
    lonSite = data["siteLong"]
    altSite = data["siteAlt"]

    ### Calculate site ECEF coordinates
    # Calculate auxiliary terms
    auxC = c.rE / sqrt(1 - (c.eE**2)*sin(radians(latSite))**2)
    auxS = auxC * (1 - c.eE**2)

    # Calculate vertical and horizontal components
    rDelta = (auxC + altSite/1000) * cos(radians(latSite))
    rKappa = (auxS + altSite/1000) * sin(radians(latSite))

    # Calculate site position in ECEF
    rSite_Ecef = [rDelta * cos(radians(lonSite)), rDelta * sin(radians(lonSite)), rKappa]

    # To convert site coordinates from ECEF to ECI, we need the GMST at each of the 3 JDs
    gmst1 = gmst(data["jd"][0])
    gmst2 = gmst(data["jd"][1])
    gmst3 = gmst(data["jd"][2])

    ### Convert site position into ECI at each epoch (JD)
    # Note: We apply a negative angle here because we are doing a passive rotation, whereby
    # we need to "undo" the rotation of the Earth.
    rSite1_Eci = ecef2eci(np.array(rSite_Ecef),-gmst1)
    rSite2_Eci = ecef2eci(np.array(rSite_Ecef),-gmst2)
    rSite3_Eci = ecef2eci(np.array(rSite_Ecef),-gmst3)
    
    # Compile the 3 ECI site vectors into a single matrix
    rSite_Eci = np.array([rSite1_Eci, rSite2_Eci, rSite3_Eci])

    ### Calculate LOS vectors for 3 observations
    L = []
    for i in range(3):
        los = calcLosUnitVector(data["ra"][i], data["dec"][i])
        L.append(los)
    L = np.transpose(L)

    # Create a satellite state object
    svState = SatState(0.0,0.0)

    ### Call appropriate orbit determination algorithm
    if data["odtype"] == "gauss":
        print("Running Gauss Algorithm")
        gaussAlgo(svState,L,data["jd"],rSite_Eci)
        print("r2: ",svState.r)
        print("v2: ",svState.v)
    elif data["odtype"] == "double-r":
        print("Running Double-R Iteration Algorithm")
    else:
        print("Error: Unspecified or incorrect OD Type.")

if __name__ == "__od__":
    orbitDetermination()