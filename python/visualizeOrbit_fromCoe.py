import numpy as np
import matplotlib.pyplot as plt

def coe_to_r(a, e, i, raan, argp, nu):
    """
    Convert classical orbital elements to position vector in ECI frame
    Angles in degrees
    """
    # Convert degrees to radians
    i = np.radians(i)
    raan = np.radians(raan)
    argp = np.radians(argp)
    nu = np.radians(nu)

    # Distance from central body
    r_mag = a * (1 - e**2) / (1 + e * np.cos(nu))

    # Position in perifocal (PQW) frame
    r_pqw = np.array([
        r_mag * np.cos(nu),
        r_mag * np.sin(nu),
        0
    ])

    # Rotation matrices
    R3_W = np.array([
        [ np.cos(raan), -np.sin(raan), 0],
        [ np.sin(raan),  np.cos(raan), 0],
        [ 0,             0,            1]
    ])

    R1_i = np.array([
        [1, 0,           0],
        [0, np.cos(i), -np.sin(i)],
        [0, np.sin(i),  np.cos(i)]
    ])

    R3_w = np.array([
        [ np.cos(argp), -np.sin(argp), 0],
        [ np.sin(argp),  np.cos(argp), 0],
        [ 0,             0,            1]
    ])

    # Full rotation matrix: PQW → ECI
    Q = R3_W @ R1_i @ R3_w

    return Q @ r_pqw


def plot_orbit(a, e, i, raan, argp, num_points=500):
    nus = np.linspace(0, 360, num_points)
    orbit = np.array([coe_to_r(a, e, i, raan, argp, nu) for nu in nus])

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    ax.plot(orbit[:, 0], orbit[:, 1], orbit[:, 2], label="Orbit")
    ax.scatter(0, 0, 0, color="orange", s=100, label="Central Body")

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Orbital Visualization")
    ax.legend()
    ax.set_box_aspect([1, 1, 1])

    plt.show()


# -------------------------------
# Example usage (Earth-like orbit)
# -------------------------------
a = 7000        # km
e = 0.2
i = 45         # degrees
raan = 30      # degrees
argp = 40      # degrees

plot_orbit(a, e, i, raan, argp)
