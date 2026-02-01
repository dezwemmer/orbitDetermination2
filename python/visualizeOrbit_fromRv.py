import numpy as np
import matplotlib.pyplot as plt

def plot_eci_positions(r_eci):
    """
    r_eci: Nx3 array of ECI position vectors [km]
    """
    r_eci = np.array(r_eci)

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    ax.plot(r_eci[:, 0], r_eci[:, 1], r_eci[:, 2], label="Orbit")
    ax.scatter(0, 0, 0, color="orange", s=100, label="Central Body")

    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    ax.set_zlabel("Z (km)")
    ax.set_title("ECI Orbit Visualization")
    ax.legend()
    ax.set_box_aspect([1, 1, 1])

    plt.show()


# Example usage
r_eci = [
    [7000, 0, 0],
    [6900, 1000, 200],
    [6500, 2000, 600],
    [5800, 3000, 1200],
]

plot_eci_positions(r_eci)
