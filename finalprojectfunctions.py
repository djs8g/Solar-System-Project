import numpy as np
import math
import random
import finalprojectconstants as fpconst
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython import display 

def Nbody_derivatives(pos, vel):
# given N bodies, implement the force equations as stated above.
    dpdt = vel
    dvdt = np.zeros(vel.shape)
    for i in range(N_bodies) :
        for j in range(N_bodies) :
            if i == j : 
                continue
            r = np.linalg.norm(pos[j]-pos[i])
            # 𝐯𝑖 = 𝑑𝐱𝑖/𝑑t
            # 𝐅𝑖 = 𝑚𝑖 * 𝑑𝐯𝑖/𝑑𝑡
            r_hat = (pos[j] - pos[i]) / r
            # 𝑟̂ = (𝐫𝑖 − 𝐫𝑗) / |𝐫𝑖 − 𝐫𝑗|
            r_hat=(pos[j] - pos[i])/r
            dvdt[i] += ( M[j] / r**2 ) * r_hat  
    return dpdt, dvdt

def Nbody_derivatives2(pos, vel):
    dpdt = vel
    dvdt = np.zeros(vel.shape)
    rvec = pos[np.newaxis,:,:] - pos[:,np.newaxis,:] 
    r = np.maximum(np.linalg.norm(rvec,axis=-1),1e-30)
    dvdt = -((M[np.newaxis,:]/(r*r*r))[:,:,np.newaxis]*rvec).sum(axis=1)
    return dpdt, dvdt

def Nbody_derivatives3(pos, vel):
    dpdt = vel
    dvdt = np.zeros_like(vel)
    for i in range(N_bodies) :
        for j in range(N_bodies):
            if i == j: 
                continue
            r = np.linalg.norm(pos[j]-pos[i])
            mass = M[j]
            rhat = (pos[j] - pos[i])/r
            dvdt[i] += -mass/(r*r)*rhat
        
    return dpdt, dvdt

def run_Nbody_rk2(tend, tframe, dt):
    y = initial_conditions()
    p = y[:,0:3]
    v = y[:,3:6]
    t = 0
    tnext = tframe
    positions = []
    while t<tend :
        while t < tnext :
            # compute using rk2
            delta_t = min(tnext-t,dt)
            dpdt, dvdt = Nbody_derivatives(p, v) 
            phalf, vhalf = p+dpdt*0.5*delta_t, v+dvdt*0.5*delta_t
            dpdt, dvdt = Nbody_derivatives(phalf, vhalf)
            p, v = p + dpdt*delta_t, v + dvdt*delta_t
            t += delta_t
        positions.append(p.copy())
        tnext += tframe
    return positions

def rk2animate(i, rk2positions):
    axrk2.clear()
    # Get the point from the points list at index i
    pos = rk2positions[i]
    axrk2.scatter(pos[0,0], pos[0,1], color='green', marker='o')
    axrk2.scatter(pos[1,0], pos[1,1], color='red', marker='o')
    # Set the x and y axis to display a fixed range
    axrk2.set_xlim([-5, 5])
    axrk2.set_ylim([-5, 5])

def run_Nbody_leapfrog(tend, tframe, dt):
    y = initial_conditions()
    p = y[:,0:3]
    v = y[:,3:6]
    t = 0
    tnext = tframe
    positions = []
    while t < tend:
        while t < tnext:
            # Leapfrog integration
            phalf = p + 0.5 * dt * v
            _, dvdt = Nbody_derivatives(phalf, v)
            v = v + dt * dvdt
            p = p + dt * v
            t += dt
        positions.append(p.copy())
        tnext += tframe
    return positions

# Animation function to use with leapfrog integrator
def lfanimate(i, lfpositions):
    axlf.clear()
    # Get the point from the points list at index i
    pos = lfpositions[i]
    axlf.scatter(pos[0,0], pos[0,1], color='green', marker='o')
    axlf.scatter(pos[1,0], pos[1,1], color='red', marker='o')
    # Set the x and y axis to display a fixed range
    axlf.set_xlim([-5, 5])
    axlf.set_ylim([-5, 5])

# Function to calculate gravitational force
def gravitational_force(mass1, mass2, distance):
    return fpconst.G * mass1 * mass2 / distance**2

# Calculate the gravitational force exerted on object 1 by object 2
def gravforce(obj1, obj2):
    # Calculate distance vector between object 1 and object 2
    r_vec = obj1.pos_vector - obj2.pos_vector
    # Calculate magnitude of distance vector
    r_mag = mag(r_vec)
    # Calcualte unit vector of distance vector
    r_hat = r_vec / r_mag
    # Calculate force magnitude
    force_mag = fpconst.G * obj1.mass * obj2.mass / r_mag**2
    # Calculate force vector
    force_vec = -force_mag * r_hat

    return force_vec

# Calculate the momenta influence of object 2 on object 1
def calculate_momenta(obj1, obj2, dt):
    obj1.force += gravforce(obj1, obj2)
    momentum_update += obj1.force * dt
    return momentum_update

# Update the position vector due to the influence of object 2 on object 1
def update_pos(obj1, obj2, dt):
    obj1.momentum += calculate_momenta(obj1, obj2, dt)
    pos_vector_update += obj1.momentum / obj1.mass * dt
    return pos_vector_update

# Function to update the positions and velocities of the planetary bodies
def update_positions_and_velocities(bodies, sun, time_step):
    num_bodies = len(bodies)
    num_dimensions = fpconst.dimensions  # Assuming 2D motion
    accelerations = np.zeros((num_bodies, num_dimensions))

    # Calculate accelerations due to interactions between planets
    for i in range(num_bodies):
        # Calculate interaction with the Sun
        bodies[i].pos_vector += update_pos(bodies[i], sun, time_step)
#        delta_x_pos = abs(bodies[i].getXPos() - sun.getXPos())
#        delta_y_pos = abs(bodies[i].getYPos() - sun.getYPos())
#        distance = math.sqrt(delta_x_pos**2 + delta_y_pos**2)
#        force_magnitude = gravitational_force(bodies[i].getMass(), sun.getMass(), distance)
#        force_direction_x, force_direction_y = calculate_force_direction(delta_x_pos, delta_y_pos, distance)
#        acceleration_x = force_direction_x * force_magnitude / bodies[i].getMass()
#        acceleration_y = force_direction_y * force_magnitude / bodies[i].getMass()
#        bodies[i].setAccelerationX(acceleration_x)
#        bodies[i].setAccelerationY(acceleration_y)
#        accelerations[i, 0] += acceleration_x
#        accelerations[i, 1] += acceleration_y
        
        # Calculate interactions with other planets
        for j in range(num_bodies):
            if i != j:
                bodies[i].pos_vector += update_pos(bodies[i], bodies[j], time_step)
#                delta_x_pos = abs(bodies[i].getXPos() - bodies[j].getXPos())
#                delta_y_pos = abs(bodies[i].getYPos() - bodies[j].getYPos())
#                distance = math.sqrt(delta_x_pos**2 + delta_y_pos**2)
#                force_magnitude = gravitational_force(bodies[i].getMass(), bodies[j].getMass(), distance)
#                force_direction_x, force_direction_y = calculate_force_direction(delta_x_pos, delta_y_pos, distance)
#                acceleration_x = force_direction_x * force_magnitude / bodies[i].getMass()
#                acceleration_y = force_direction_y * force_magnitude / bodies[i].getMass()
#                bodies[i].setAccelerationX(acceleration_x)
#                bodies[i].setAccelerationY(acceleration_y)
#                accelerations[i, 0] += acceleration_x
#                accelerations[i, 1] += acceleration_y

    # Update velocities and positions
    for i in range(num_bodies):
#        bodies[i].setXVel(bodies[i].getXVel() + accelerations[i, 0] * time_step)
#        bodies[i].setYVel(bodies[i].getYVel() + accelerations[i, 1] * time_step)
        new_x_pos = bodies[i].getXPos() + bodies[i].pos_vector[0]
        new_y_pos = bodies[i].getYPos() + bodies[i].pos_vector[1]
        new_z_pos = bodies[i].getZPos() + bodies[i].pos_vector[2]
        bodies[i].moveTo(new_x_pos, new_y_pos, new_z_pos)

    return bodies

# Function to calculate the direction of force between two bodies
def calculate_force_direction(delta_x, delta_y, distance):    
    # Handle the case when bodies are at the same position
    if distance == 0:
        return None

    # Normalize the vector to get direction
    direction_x = delta_x / distance
    direction_y = delta_y / distance

    # Return the direction vector as a tuple
    return direction_x, direction_y

# Function to calculate the direction of force between two bodies
def calculate_force_magnitude(delta_x, delta_y, distance, force_magnitude):
    # Calculate the direction vector
    direction_x, direction_y = calculate_force_direction(delta_x, delta_y, distance)

    # If bodies are at the same position, force magnitude is zero
    if distance == 0:
        return None

    # Calculate the force vector
    force_x = direction_x * force_magnitude
    force_y = direction_y * force_magnitude

    return force_x, force_y

# Randomize starting positions (along circular orbits)
def randomize_starting_positions(body):
    theta = random.uniform(0, 2 * math.pi)
    x_pos = body.orbital_radius * math.cos(theta)
    y_pos = body.orbital_radius * math.sin(theta)
#    positions = [x_pos, y_pos]
    vel = 2 * math.pi * body.orbital_radius / body.orbital_period
    x_vel = -vel * math.sin(theta)
    y_vel = vel * math.cos(theta)
#    velocities = [-vel * math.sin(theta), vel * math.cos(theta)]
    return x_pos, y_pos, x_vel, y_vel


# Define the differential equations for the motion of the planets
def planetary_motion(t, y):
    dydt = np.zeros_like(y)
    for i, (planet, _) in enumerate(y):
        dydt[i, :] = y[i, 2:]
        for j, (other_planet, _) in enumerate(y):
            if i != j:
                r = np.linalg.norm(y[i, :2] - y[j, :2])
                dydt[i, 2:] += -fpconst.G * fpconst.masses[other_planet] * (y[i, :2] - y[j, :2]) / r**3
    return dydt.flatten()

import itertools
import math
import matplotlib.pyplot as plt

from solarsystemvectorclass import Vector

class SolarSystem:
    def __init__(self, size, projection_2d=False):
        self.size = size
        self.projection_2d = projection_2d
        self.bodies = []

        self.fig, self.ax = plt.subplots(
            1,
            1,
            subplot_kw={"projection": "3d"},
            figsize=(self.size / 50, self.size / 50),
        )
        self.fig.tight_layout()
        if self.projection_2d:
            self.ax.view_init(10, 0)
        else:
            self.ax.view_init(0, 0)

    def add_body(self, body):
        self.bodies.append(body)

    def update_all(self):
        self.bodies.sort(key=lambda item: item.position[0])
        for body in self.bodies:
            body.move()
            body.draw()

    def draw_all(self):
        self.ax.set_xlim((-self.size / 2, self.size / 2))
        self.ax.set_ylim((-self.size / 2, self.size / 2))
        self.ax.set_zlim((-self.size / 2, self.size / 2))
        if self.projection_2d:
            self.ax.xaxis.set_ticklabels([])
            self.ax.yaxis.set_ticklabels([])
            self.ax.zaxis.set_ticklabels([])
        else:
            self.ax.axis(False)
        plt.pause(0.001)
        self.ax.clear()

    def calculate_all_body_interactions(self):
        bodies_copy = self.bodies.copy()
        for idx, first in enumerate(bodies_copy):
            for second in bodies_copy[idx + 1:]:
                first.accelerate_due_to_gravity(second)

class SolarSystemBody:
    min_display_size = 10
    display_log_base = 1.3

    def __init__(
        self,
        solar_system,
        mass,
        position=(0, 0, 0),
        velocity=(0, 0, 0),
    ):
        self.solar_system = solar_system
        self.mass = mass
        self.position = position
        self.velocity = Vector(*velocity)
        self.display_size = max(
            math.log(self.mass, self.display_log_base),
            self.min_display_size,
        )
        self.colour = "black"

        self.solar_system.add_body(self)

    def move(self):
        self.position = (
            self.position[0] + self.velocity[0],
            self.position[1] + self.velocity[1],
            self.position[2] + self.velocity[2],
        )

    def draw(self):
        self.solar_system.ax.plot(
            *self.position,
            marker="o",
            markersize=self.display_size + self.position[0] / 30,
            color=self.colour
        )
        if self.solar_system.projection_2d:
            self.solar_system.ax.plot(
                self.position[0],
                self.position[1],
                -self.solar_system.size / 2,
                marker="o",
                markersize=self.display_size / 2,
                color=(.5, .5, .5),
            )

    def accelerate_due_to_gravity(self, other):
        distance = Vector(*other.position) - Vector(*self.position)
        distance_mag = distance.get_magnitude()

        force_mag = self.mass * other.mass / (distance_mag ** 2)
        force = distance.normalize() * force_mag

        reverse = 1
        for body in self, other:
            acceleration = force / body.mass
            body.velocity += acceleration * reverse
            reverse = -1

class Sun(SolarSystemBody):
    def __init__(
        self,
        solar_system,
        mass=10_000,
        position=(0, 0, 0),
        velocity=(0, 0, 0),
    ):
        super(Sun, self).__init__(solar_system, mass, position, velocity)
        self.colour = "yellow"

class Planet(SolarSystemBody):
    colours = itertools.cycle([(1, 0, 0), (0, 1, 0), (0, 0, 1)])

    def __init__(
        self,
        solar_system,
        mass=10,
        position=(0, 0, 0),
        velocity=(0, 0, 0),
    ):
        super(Planet, self).__init__(solar_system, mass, position, velocity)
        self.colour = next(Planet.colours)


import math

class Vector:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"Vector({self.x}, {self.y}, {self.z})"

    def __str__(self):
        return f"{self.x}i + {self.y}j + {self.z}k"

    def __getitem__(self, item):
        if item == 0:
            return self.x
        elif item == 1:
            return self.y
        elif item == 2:
            return self.z
        else:
            raise IndexError("There are only three elements in the vector")

    def __add__(self, other):
        return Vector(
            self.x + other.x,
            self.y + other.y,
            self.z + other.z,
        )

    def __sub__(self, other):
        return Vector(
            self.x - other.x,
            self.y - other.y,
            self.z - other.z,
        )

    def __mul__(self, other):
        if isinstance(other, Vector):  # Vector dot product
            return (
                self.x * other.x
                + self.y * other.y
                + self.z * other.z
            )
        elif isinstance(other, (int, float)):  # Scalar multiplication
            return Vector(
                self.x * other,
                self.y * other,
                self.z * other,
            )
        else:
            raise TypeError("operand must be Vector, int, or float")

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return Vector(
                self.x / other,
                self.y / other,
                self.z / other,
            )
        else:
            raise TypeError("operand must be int or float")

    def get_magnitude(self):
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def normalize(self):
        magnitude = self.get_magnitude()
        return Vector(
            self.x / magnitude,
            self.y / magnitude,
            self.z / magnitude,
        )