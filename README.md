# Orbit Determination 2
Remake of a previous OD project.

Orbit propagation & conjunction under disturbances of earth oblateness and atmospheric drag perturbations.
Convert initial OE at epoch t0 into initial Cartesian state (r0, v0) in ECI.
Propagate state ahead to final epoch: JD 2454873.205555555
Orbit propagation is according to numerical integation of Newtonian EoM governing perturbed Keplerian motion.


Process:
1 Import TLE data from 2 satellites.
    a. Julian Date, RA, Declination, Local Sidereal Time (LST)
    b. 3 measurements each. 
2 Import site information (observation station was at U of Az ground station)
