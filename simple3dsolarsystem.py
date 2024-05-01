import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# Constants
G = 6.67430e-11  # gravitational constant
dt = 3600 * 24    # time step (1 day)
num_steps = 365   # number of steps

# Define initial conditions for celestial bodies (mass, position, velocity)
bodies = {
    "Sun": {"mass": 1.989e30, "pos": np.array([0, 0, 0]), "vel": np.array([0, 0, 0])},
    "Mercury": {"mass": 3.301e23, "pos": np.array([0, 0.39e12, 0]), "vel": np.array([47000, 0, 0])},
    "Venus": {"mass": 4.867e24, "pos": np.array([0.723e12, 0, 0]), "vel": np.array([0, 35000, 0])},
    "Earth": {"mass": 5.972e24, "pos": np.array([1.0e12, 0, 0]), "vel": np.array([0, 29783, 0])},
    "Mars": {"mass": 6.39e23, "pos": np.array([1.524e12, 0, 0]), "vel": np.array([0, 24100, 0])},
}

# Simulation function
def simulate_step():
    for name1, body1 in bodies.items():
        force = np.array([0, 0, 0])
        for name2, body2 in bodies.items():
            if name1 != name2:
                r = body2["pos"] - body1["pos"]
                force += G * body1["mass"] * body2["mass"] * r / np.linalg.norm(r)**3
        body1["vel"] += force / body1["mass"] * dt
        body1["pos"] += body1["vel"] * dt

# Animation function
def update(frame):
    simulate_step()
    ax.clear()
    ax.set_title(f"Day {frame}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    for name, body in bodies.items():
        ax.scatter(body["pos"][0], body["pos"][1], body["pos"][2], label=name)
    ax.legend()

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

# Animate the simulation
ani = FuncAnimation(fig, update, frames=num_steps, interval=50)

plt.show()
