import turtle
import pygame as pg
from OpenGL.GL import *
import math
import random
import finalprojectconstants as const
import finalprojectfunctions as func

conversion_ratio = const.conversion_factor

class SolarSystem:
    def __init__(self, width, height):
        self.thesun = None
        self.planets = []
        self.ssturtle = turtle.Turtle()
        self.ssturtle.hideturtle()
        self.ssscreen = turtle.Screen()
        self.size_ratio = self.getMaxOrbitalRadius() / const.conversion_factor

        # Calculate screen width and height
#        screen_width = width / (const.screen_ratio * self.getMaxOrbitalRadius())
#        screen_height = height / (const.screen_ratio * self.getMaxOrbitalRadius())

#        self.ssscreen.setworldcoordinates(-width, -height, width, height)
        self.ssscreen.setup(width, height)
        self.ssscreen.tracer(const.tracer_value)

        # Add event listener for zooming
        self.ssscreen.onkey(self.zoom_in, 'Up')
        self.ssscreen.onkey(self.zoom_out, 'Down')
        self.ssscreen.listen()
    
    # Calculate maximum orbital radius
    def getMaxOrbitalRadius(self):
        max_orbital_radius = max([const.orbital_radii[name] for name in const.planets])
        print("Max Orbital Radius=", max_orbital_radius)
        return max_orbital_radius

    def addPlanet(self, aplanet):
        self.planets.append(aplanet)

    def addSun(self, asun):
        self.thesun = asun

    def showPlanets(self):
        for aplanet in self.planets:
            print(aplanet)

    def freeze(self):
        self.ssscreen.exitonclick()

    def movePlanets(self):
        G = const.G
        dt = const.time_step

        # Update positions and velocities of all planets, considering interactions
        self.planets = func.update_positions_and_velocities(self.planets, self.thesun, dt)

        self.ssscreen.update()  # Update the screen after moving the planets

    def zoom_in(self):
        self.ssscreen.setworldcoordinates(self.ssscreen.minx() * 0.9, self.ssscreen.miny() * 0.9,
                                           self.ssscreen.maxx() * 0.9, self.ssscreen.maxy() * 0.9)

    def zoom_out(self):
        self.ssscreen.setworldcoordinates(self.ssscreen.minx() * 1.1, self.ssscreen.miny() * 1.1,
                                           self.ssscreen.maxx() * 1.1, self.ssscreen.maxy() * 1.1)

class Sun:
    def __init__(self, name, diameter, mass, temp, size_ratio):
        self.name = name
        self.diameter = diameter
        self.mass = mass
        self.temp = temp
        self.size = self.diameter / size_ratio  # Ratio to adjust the size of the Sun
        self.x_pos = 0
        self.y_pos = 0
        self.z_pos = 0
        self.pos_vector = (self.x_pos, self.y_pos, self.z_pos)

        self.sturtle = turtle.Turtle()
        self.sturtle.shape("circle")
        self.sturtle.shapesize(self.size)  # Adjust size based on diameter and ratio
        self.sturtle.color("yellow")

    def getName(self):
        return self.name

    def getRadius(self):
        return self.diameter / 2.0

    def getMass(self):
        return self.mass

    def getTemperature(self):
        return self.temp

    def getVolume(self):
        volume = 4.0 / 3 * math.pi * getRadius(self)**3
        return volume

    def getSurfaceArea(self):
        surface_area = 4.0 * math.pi * getRadius(self)**2
        return surface_area

    def getDensity(self):
        density = self.mass / self.getVolume()
        return density

    def setName(self, newname):
        self.name = newname

    def __str__(self):
        return self.name

    def getXPos(self):
        return self.x_pos

    def getYPos(self):
        return self.y_pos

    def getZPos(self):
        return self.z_pos

class Planet:
    def __init__(self, name, mass, diameter, density, gravity, escape_velocity, rotational_period,
                 length_of_day, orbital_radius, perihelion, aphelion, orbital_period, orbital_velocity,
                 orbital_inclination, orbital_eccentricity, obliquity, mean_temp, surface_pressure,
                 num_moons, ring_system, global_magnetic_field, color, size_ratio):
        self.name = name
        self.mass = mass
        self.diameter = diameter
        self.density = density
        self.gravity = gravity
        self.escape_velocity = escape_velocity
        self.rotational_period = rotational_period
        self.length_of_day = length_of_day
        self.orbital_radius = orbital_radius
        self.perihelion = perihelion
        self.aphelion = aphelion
        self.orbital_period = orbital_period
        self.orbital_velocity = orbital_velocity
        self.orbital_inclination = orbital_inclination
        self.orbital_eccentricity = orbital_eccentricity
        self.obliquity = obliquity
        self.mean_temp = mean_temp
        self.surface_pressure = surface_pressure
        self.num_moons = num_moons
        self.ring_system = ring_system
        self.global_magnetic_field = global_magnetic_field
        self.color = color
        self.size_ratio = size_ratio
        self.size = diameter / size_ratio  # Adjusted based upon size of the planet relative to the size of the Sun and the largest orbit

        # Initialize the poisitions
        self.init_x_pos, self.init_y_pos, self.init_xvel, self.init_yvel = func.randomize_starting_positions(self)
#        self.init_x_pos = self.getInitXPosition()
        self.x_pos = self.init_x_pos
#        self.init_y_pos = self.getInitYPosition()
        self.y_pos = self.init_y_pos
        self.init_z_pos = 0
        self.z_pos = self.init_z_pos
        self.pos_vector = (self.x_pos, self.y_pos, self.z_pos)

        # Initialize the velocities
#        self.init_xvel = self.getInitXVelocity()
        self.x_vel = self.init_xvel
#        self.init_yvel = self.getInitYVelocity()
        self.y_vel = self.init_yvel
        self.init_zvel = 0
        self.z_vel = self.init_zvel

        # Initialize the accelerations
        self.acceleration_x = 0
        self.acceleration_y = 0
        self.acceleration_z = 0

        # Initialize momenta
        self.momentum = 0

        # Initialize force
        self.force = 0

        # Create a Turtle object for the planet
        self.pturtle = turtle.Turtle()
        self.pturtle.up()
        self.pturtle.color(self.color)
        self.pturtle.shape("circle")
        self.pturtle.shapesize(self.size)  # Adjust size based on diameter and ratio
#        self.pturtle.goto(self.x_pos, self.y_pos, self.z_pos)
        self.pturtle.goto(self.x_pos, self.y_pos)
        self.pturtle.down()

    def getName(self):
        return self.name

    def getRadius(self):
        return self.radius

    def getMass(self):
        return self.mass

#    def getInitXPosition(self):
        # Randomly generate x position within the range of orbital radius
#        return random.uniform(-self.orbital_radius, self.orbital_radius)

#    def getInitYPosition(self):
        # Randomly generate y position within the range of orbital radius
#        return random.uniform(-self.orbital_radius, self.orbital_radius)

    def getInitXVelocity(self):
        # Calculate the initial x velocity based on orbital velocity
        orbital_speed = self.orbital_velocity
        orbital_radius = math.sqrt(self.init_x_pos ** 2 + self.init_y_pos ** 2)
        return -orbital_speed * (self.init_x_pos / orbital_radius)

    def getInitYVelocity(self):
        # Calculate the initial y velocity based on orbital velocity
        orbital_speed = self.orbital_velocity
        orbital_radius = math.sqrt(self.init_x_pos ** 2 + self.init_y_pos ** 2)
        return -orbital_speed * (self.init_y_pos / orbital_radius)
    
    def getOrbitalRadius(self):
        return self.orbital_radius

    def getVolume(self):
        volume = 4.0 / 3 * math.pi * self.radius**3
        return volume

    def getSurfaceArea(self):
        surface_area = 4.0 * math.pi * self.radius**2
        return surface_area

    def getDensity(self):
        density = self.mass / self.getVolume()
        return density

    def setName(self, newname):
        self.name = newname

    def show(self):
        print(self.name)   

    def __str__(self):
        return self.name

    def moveTo(self, newx, newy):
        self.x = newx
        self.y = newy
        self.z = newz
        self.pturtle.goto(newx, newy)

    def getXPos(self):
        return self.x_pos

    def getYPos(self):
        return self.y_pos

    def getZPos(self):
        return self.z_pos

    def getXVel(self):
        return self.x_vel

    def getYVel(self):
        return self.y_vel

    def getZVel(self):
        return self.z_vel

    def setXVel(self, new_xvel):
        self.x_vel = new_xvel

    def setYVel(self, new_yvel):
        self.y_vel = new_yvel

    def setZVel(self, new_zvel):
        self.z_vel = new_zvel

    def getAccelerationX(self):
        return self.acceleration_x

    def getAccelerationY(self):
        return self.acceleration_y

    def getAccelerationZ(self):
        return self.acceleration_z

    def setAccelerationX(self, new_acceleration_x):
        self.acceleration_x = new_acceleration_x

    def setAccelerationY(self, new_acceleration_y):
        self.acceleration_y = new_acceleration_y

    def setAccelerationY(self, new_acceleration_z):
        self.acceleration_z = new_acceleration_z

def createSolarSystem():
    solarsystem = SolarSystem(const.num_pixels, const.num_pixels)

    sun = Sun(const.nameofStar, const.diameterofSun, const.massofSun, const.tempofSun, solarsystem.size_ratio)
    solarsystem.addSun(sun)

    # Iterate over the planet names and assign attributes
    for name in const.planets:
        planet = Planet(
            name,
            const.masses[name],
            const.diameter[name],
            const.density[name],
            const.gravity[name],
            const.escape_velocity[name],
            const.rotational_period[name],
            const.length_of_day[name],
            const.orbital_radii[name],
            const.perihelion[name],
            const.aphelion[name],
            const.orbital_period[name],
            const.orbital_velocity[name],
            const.orbital_inclination[name],
            const.orbital_eccentricity[name],
            const.obliquity[name],
            const.mean_temp[name],
            const.surface_pressure[name],
            const.num_moons[name],
            const.ring_system[name],
            const.global_magnetic_field[name],
            const.color[name],
            size_ratio = solarsystem.size_ratio
        )
        solarsystem.addPlanet(planet)

    return solarsystem

def AnimateSolarSystem():
    SolarSystem = createSolarSystem()
    SolarSystem.ssscreen.bgcolor("black")  # Set background color to black

    numTimePeriods = 2000
    for amove in range(numTimePeriods):
        SolarSystem.movePlanets()
        print("Loop Check!")
        print("amove=",amove)
        time.sleep(0.1)  # Adjust the delay time as needed

    SolarSystem.freeze()

AnimateSolarSystem()

turtle.done()