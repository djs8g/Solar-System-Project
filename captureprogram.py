import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)

# Function to calculate gravitational acceleration
def gravitational_acceleration(r, M):
    # r: distance vector from the object to the larger body
    # M: mass of the larger body
    return -G * M / np.linalg.norm(r)**3 * r

# Function to simulate the trajectory of the object
def simulate_trajectory(r0, v0, M, dt, t_max):
    # r0: initial position vector of the object
    # v0: initial velocity vector of the object
    # M: mass of the larger body
    # dt: time step for numerical integration
    # t_max: maximum simulation time
    
    # Initialize arrays to store position and velocity
    num_steps = int(t_max / dt)
    t_values = np.linspace(0, t_max, num_steps)
    r_values = np.zeros((num_steps, 2))
    v_values = np.zeros((num_steps, 2))
    
    # Set initial conditions
    r_values[0] = r0
    v_values[0] = v0
    
    # Perform numerical integration using Euler's method
    for i in range(1, num_steps):
        # Calculate gravitational acceleration
        a = gravitational_acceleration(r_values[i-1], M)
        
        # Update velocity and position using Euler's method
        v_values[i] = v_values[i-1] + a * dt
        r_values[i] = r_values[i-1] + v_values[i] * dt
        
    return t_values, r_values, v_values

# Function to check for capture into a stable orbit
def check_for_capture(r_values, v_values):
    # Check if the object's trajectory crosses a threshold distance from the larger body
    # and if its velocity is within a certain range for capture
    
    # Define threshold distance and capture velocity range (for demonstration purposes)
    threshold_distance = 1.5e7  # Threshold distance from the larger body (m)
    capture_velocity_min = 0.9  # Minimum capture velocity (relative to circular orbit)
    capture_velocity_max = 1.1  # Maximum capture velocity (relative to circular orbit)
    
    # Calculate distance from the larger body at each time step
    distances = np.linalg.norm(r_values, axis=1)
    
    # Find the index where the object crosses the threshold distance
    capture_index = np.argmax(distances > threshold_distance)
    
    if capture_index == 0:
        return False  # Object did not cross threshold distance
    
    # Check if the velocity at capture is within the specified range
    v_capture = np.linalg.norm(v_values[capture_index])
    v_circular_orbit = np.sqrt(G * M / threshold_distance)  # Velocity for circular orbit
    if capture_velocity_min * v_circular_orbit <= v_capture <= capture_velocity_max * v_circular_orbit:
        return True  # Object captured into stable orbit
    else:
        return False  # Object did not meet capture criteria

# Main function to run the simulation
def main():
    # Initial conditions
    r0 = np.array([1.0e7, 0])  # Initial position vector of the object (m)
    v0 = np.array([0, 1.0e3])  # Initial velocity vector of the object (m/s)
    M = 5.972e24  # Mass of the larger body (kg)
    dt = 1000  # Time step for numerical integration (s)
    t_max = 100000  # Maximum simulation time (s)
    
    # Simulate trajectory
    t_values, r_values, v_values = simulate_trajectory(r0, v0, M, dt, t_max)
    
    # Check for capture
    captured = check_for_capture(r_values, v_values)
    
    # Plot trajectory
    plt.plot(r_values[:, 0], r_values[:, 1])
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Object Trajectory')
    plt.grid(True)
    plt.axis('equal')
    plt.show()
    
    if captured:
        print("Object captured into stable orbit")
    else:
        print("Object not captured into stable orbit")

# Run the simulation
main()
