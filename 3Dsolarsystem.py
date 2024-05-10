import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define constants
G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
masses = {
    'Sun': 1.989e30,       # Mass of the Sun in kg
    'Mercury': 3.285e23,   # Mass of Mercury in kg
    'Venus': 4.867e24,     # Mass of Venus in kg
    'Earth': 5.972e24,     # Mass of Earth in kg
    'Mars': 6.39e23,       # Mass of Mars in kg
    'Jupiter': 1.898e27,   # Mass of Jupiter in kg
    'Saturn': 5.683e26,    # Mass of Saturn in kg
    'Uranus': 8.681e25,    # Mass of Uranus in kg
    'Neptune': 1.024e26    # Mass of Neptune in kg
}

# Initial positions and velocities (m and m/s)
positions = {
    'Sun': np.array([0.0, 0.0, 0.0]),
    'Mercury': np.array([0.39e12, 0.0, 0.0]),
    'Venus': np.array([0.72e12, 0.0, 0.0]),
    'Earth': np.array([1.0e12, 0.0, 0.0]),
    'Mars': np.array([1.52e12, 0.0, 0.0]),
    'Jupiter': np.array([5.2e12, 0.0, 0.0]),
    'Saturn': np.array([9.58e12, 0.0, 0.0]),
    'Uranus': np.array([19.18e12, 0.0, 0.0]),
    'Neptune': np.array([30.07e12, 0.0, 0.0])
}

velocities = {
    'Sun': np.array([0.0, 0.0, 0.0]),
    'Mercury': np.array([0.0, 4.79e4, 0.0]),
    'Venus': np.array([0.0, 3.5e4, 0.0]),
    'Earth': np.array([0.0, 2.98e4, 0.0]),
    'Mars': np.array([0.0, 2.41e4, 0.0]),
    'Jupiter': np.array([0.0, 1.3e4, 0.0]),
    'Saturn': np.array([0.0, 9.69e3, 0.0]),
    'Uranus': np.array([0.0, 6.8e3, 0.0]),
    'Neptune': np.array([0.0, 5.43e3, 0.0])
}

# Function to calculate orbital periods
def orbital_period(semimajor_axis, mass):
    return 2 * np.pi * np.sqrt(semimajor_axis**3 / (G * mass))

# Calculate orbital periods of planets
orbital_periods = {planet: orbital_period(np.linalg.norm(positions[planet]), masses[planet]) for planet in masses.keys()}

# Set simulation time to the multiple of the longest orbital period
t_max = max(orbital_periods.values()) * 2  # Set to 2 times the longest orbital period
dt = 3600 * 24 * 30 * 5  # Increased timestep to 5 months

# Initialize arrays to store positions
num_steps = int(t_max / dt)
positions_history = {planet: np.zeros((num_steps, 3)) for planet in masses}

# Simulation loop
for step in range(num_steps):
    for planet1 in masses:
        acceleration = np.array([0.0, 0.0, 0.0])
        for planet2 in masses:
            if planet1 != planet2:
                r = positions[planet2] - positions[planet1]
                r_mag = np.linalg.norm(r)
                acceleration += G * masses[planet2] / r_mag**3 * r
        velocities[planet1] += acceleration * dt
        positions[planet1] += velocities[planet1] * dt
        positions_history[planet1][step] = positions[planet1].copy()

# 3D plot showing the orbits
fig = plt.figure(figsize=(10, 8))
ax3d = fig.add_subplot(111, projection='3d')

for planet, pos_history in positions_history.items():
    ax3d.plot(pos_history[:, 0], pos_history[:, 1], pos_history[:, 2], label=planet)

ax3d.set_xlabel('X')
ax3d.set_ylabel('Y')
ax3d.set_zlabel('Z')
ax3d.set_title('Solar System Orbits (3D)')
ax3d.legend()

# 2D plot showing the orbits
fig2d, ax2d = plt.subplots(figsize=(8, 8))
for planet, pos_history in positions_history.items():
    ax2d.plot(pos_history[:, 0], pos_history[:, 1], label=planet)

ax2d.set_xlabel('X')
ax2d.set_ylabel('Y')
ax2d.set_title('Solar System Orbits (2D)')
ax2d.legend()

plt.show()
