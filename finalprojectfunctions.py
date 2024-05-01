import numpy as np
import finalprojectconstants as fpconst
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

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


'''
import pygame
import math
pygame.init()

WIDTH, HEIGHT = 800, 800
WIN = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Planet Simulator")

WHITE = (255, 255, 255)
YELLOW = (255, 255, 0)
BLUE = (100, 149, 237)
RED = (188, 39, 50)
DARK_GREY = (80, 78, 81)

FONT = pygame.font.SysFont("comicsans", 16)

class Planet:
    AU = 149.6e6 * 1000 # in meters
    G = 6.67428e-11
    SCALE = 250 / AU # 1 AU = 100 pixels
    TIMESTEP = 3600*24 # number of seconds in a day

    def __init__(self, x, y, radius, color, mass):
        self.x = x
        self.y = y 
        self.radius = radius
        self.color = color
        self.mass = mass
        self.orbit = []
        self.sun = False
        self.distance_to_sun = 0
        self.x_vel = 0
        self.y_vel = 0

    def draw(self, win):
        x = self.x * self.SCALE + WIDTH / 2
        y = self.y * self.SCALE + HEIGHT / 2
        
        if len(self.orbit) > 2:
            updated_points = []
            for point in self.orbit:
                x, y = point
                x = x * self.SCALE + WIDTH / 2
                y = y * self.SCALE + HEIGHT / 2
                updated_points.append(point)
            
            pygame.draw.lines(win, self.color, False, updated_points, 2)

        pygame.draw.circle(win, self.color, (x, y), self.radius)
        if not self.sun:
            distance_text = FONT.render(f"{round(self.distance_to_sun/1000, 1)}km", 1, WHITE)
            WIN.blit(distance_text, (x - distance_text.get_width()/2, y - distance_text.get_height()/2))
        
    def attraction(self, other):
        other_x, other_y = other.x, other.y
        distance_x = other_x - self.x
        distance_y = other_y - self.y
        distance = math.sqrt(distance_x ** 2 + distance_y ** 2)

        if other.sun:
            self.distance_to_sun = distance

        force = self.G * self.mass * other.mass / distance**2
        theta = math.atan2(distance_y, distance_x)
        force_x = math.cos(theta) * force
        force_y = math.sin(theta) * force
        return force_x, force_y

    def update_position(self, planets):
        total_fx = total_fy = 0
        for planet in planets:
            if self == planet:
                continue

            fx, fy = self.attraction(planet)
            total_fx += fx
            total_fx += fy

        self.x_vel += total_fx / self.mass * self.TIMESTEP
        self.y_vel += total_fy / self.mass * self.TIMESTEP

        self.x += self.x_vel * self.TIMESTEP
        self.y += self.y_vel * self.TIMESTEP
        self.orbit.append((self.x, self.y))

def main():
    run = True
    clock = pygame.time.Clock()

    sun = Planet(0, 0 , 30, YELLOW, 1.98892 * 10**30)
    sun.sun = True
    mercury = Planet(0.387 * Planet.AU, 0, 8, DARK_GREY, 3.30 * 10**23)
    mercury.y_vel = -47.4 * 1000
    venus = Planet(0.723 * Planet.AU, 0, 14, WHITE, 4.8685 * 10**24)
    venus.y_vel = -35.02 * 1000
    earth = Planet(-1 * Planet.AU, 0, 16, BLUE, 5.9742 * 10**24)
    earth.y_vel = 29.783 * 1000
    mars = Planet(-1.524 * Planet.AU, 0, 12, WHITE, 6.39 * 10**23)
    mars.y_vel = 24.077 * 1000

    planets = [sun, mercury, venus, earth, mars]
    
    while run:
        clock.tick(60)
        WIN.fill((0, 0, 0))
        pygame.display.update()
        
        for  event in pygame.event.get():
            if event.type == pygame.QUIT:
                run = False

        for planet in planets:
            planet.update_position(planets)
            planet.draw(WIN)
            
    pygame.quit()

main()