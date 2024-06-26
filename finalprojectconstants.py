# Define gravitational constant (m^3 kg^-1 s^-2)
G = 6.67430e-11

# Define other constants
minute = 60 # in seconds
hour = 60 * minute # in seconds
day = 24 * hour # in seconds
year = day * 365 # in seconds
dimensions = 2
time_step = hour
num_pixels = 1200
tracer_value = 1

# Define the constants for the Sun
massofSun = 1.989e30 # kg
diameterofSun = 1400000 # km
tempofSun = 5600.0 # degrees C at the surface
nameofStar = "Sun"

# Define screen ratio for animation for km to pixels to larger than the maximum orbit
screen_ratio = 1.2 

# Calculate conversion factor for km to pixels
conversion_factor = num_pixels / screen_ratio

planets = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]

# Colors for planets
color = {
    "Mercury": 'grey',
    "Venus": 'olive',
    "Earth": 'blue',
    "Mars": 'red',
    "Jupiter": 'goldenrod',
    "Saturn": 'burlywood',
    "Uranus": 'light blue',
    "Neptune": 'navy'
}

# Mass (kg) - This is the mass of the planet in kilograms.
masses = {
    "Mercury": 3.285e23,
    "Venus": 4.867e24,
    "Earth": 5.972e24,
    "Mars": 6.39e23,
    "Jupiter": 1.898e27,
    "Saturn": 5.683e26,
    "Uranus": 8.681e25,
    "Neptune": 1.024e26
}

# Diameter (km) - The diameter of the planet at the equator, the distance through the center of the planet from one point on the equator to the opposite side, in kilometers.
diameter = {
    "Mercury": 4879.0,
    "Venus": 12104.0,
    "Earth": 12756.0,
    "Mars": 6792.0,
    "Jupiter": 142984.0,
    "Saturn": 120536.0,
    "Uranus": 51118.0,
    "Neptune": 49528.0
}

# Density (kg/m3) - The average density (mass divided by volume) of the whole planet (not including the atmosphere for the terrestrial planets) in kilograms per cubic meter.
density = {
    "Mercury": 5429,
    "Venus": 5243,
    "Earth": 5514,
    "Mars": 3934,
    "Jupiter": 1326,
    "Saturn": 687,
    "Uranus": 1270,
    "Neptune": 1638
}

# Gravity (m/s2 or ft/s2) - The gravitational acceleration on the surface at the equator in meters per second squared, including the effects of rotation. For the gas giant planets the gravity is given at the 1 bar pressure level in the atmosphere.
gravity = {
    "Mercury": 3.7,
    "Venus": 8.9,
    "Earth": 9.8,
    "Mars": 3.7,
    "Jupiter": 23.1,
    "Saturn": 9.0,
    "Uranus": 8.7,
    "Neptune": 11.0
}

# Escape Velocity (km/s) - Initial velocity, in kilometers per second, needed at the surface (at the 1 bar pressure level for the gas giants) to escape the body's gravitational pull, ignoring atmospheric drag.
escape_velocity = {
    "Mercury": 4.3,
    "Venus": 10.4,
    "Earth": 11.2,
    "Mars": 5.0,
    "Jupiter": 59.5,
    "Saturn": 35.5,
    "Uranus": 21.3,
    "Neptune": 23.5
}

# Rotational Period (hours) - This is the time it takes for the planet to complete one rotation relative to the fixed background stars (not relative to the Sun) in hours. Negative numbers indicate retrograde (backwards relative to the Earth) rotation.
rotational_period = {
    "Mercury": 1407.6,
    "Venus": -5832.5,
    "Earth": 23.9,
    "Mars": 24.6,
    "Jupiter": 9.9,
    "Saturn": 10.7,
    "Uranus": -17.2,
    "Neptune": 16.1
}

# Length of Day (hours) - The average time in hours for the Sun to move from the noon position in the sky at a point on the equator back to the same position.
length_of_day = {
    "Mercury": 4222.6,
    "Venus": 2802.0,
    "Earth": 24.0,
    "Mars": 24.7,
    "Jupiter": 9.9,
    "Saturn": 10.7,
    "Uranus": 17.2,
    "Neptune": 16.1
}

# Distance from Sun (km) - This is the average distance from the planet to the Sun in kilometers, also known as the semi-major axis. All planets have orbits which are elliptical, not perfectly circular, so there is a point in the orbit at which the planet is closest to the Sun, the perihelion, and a point furthest from the Sun, the aphelion. The average distance from the Sun is midway between these two values. The average distance from the Earth to the Sun is defined as 1 Astronomical Unit (AU).
orbital_radii = {
    "Mercury": 57.9e6,
    "Venus": 108.2e6,
    "Earth": 149.6e6,
    "Mars": 228e6,
    "Jupiter": 778e6,
    "Saturn": 1432e6,
    "Uranus": 2867e6,
    "Neptune": 4515e6
}

# Perihelion (km) - The closest point in a planet's orbit about the Sun, see "Distance from Sun" above.
perihelion = {
    "Mercury": 46.0e6,
    "Venus": 107.5e6,
    "Earth": 147.1e6,
    "Mars": 206.7e6,
    "Jupiter": 740.6e8,
    "Saturn": 1357.6e6,
    "Uranus": 2732.7e6,
    "Neptune": 4471.1e6
}

# Aphelion (km) - The furthest point in a planet's orbit about the Sun, see "Distance from Sun" above.
aphelion = {
    "Mercury": 69.8e6,
    "Venus": 108.9e6,
    "Earth": 152.1e6,
    "Mars": 249.3e6,
    "Jupiter": 816.4e6,
    "Saturn": 1506.5e6,
    "Uranus": 3001.4e6,
    "Neptune": 4558.9e6
}

# Orbital Period (days) - This is the time in Earth days for a planet to orbit the Sun from one vernal equinox to the next. Also known as the tropical orbit period, this is equal to a year on Earth.
orbital_period = {
    "Mercury": 88.0,
    "Venus": 224.7,
    "Earth": 365.2,
    "Mars": 687.0,
    "Jupiter": 4331.0,
    "Saturn": 10747.0,
    "Uranus": 30589.0,
    "Neptune": 59800.0
}

# Orbital Velocity (km/s) - The average velocity or speed of the planet as it orbits the Sun, in kilometers per second.
orbital_velocity = {
    "Mercury": 47.4,
    "Venus": 35.0,
    "Earth": 29.8,
    "Mars": 24.1,
    "Jupiter": 13.1,
    "Saturn": 9.7,
    "Uranus": 6.8,
    "Neptune": 5.4
}

# Orbital Inclination (degrees) - The angle in degrees at which a planets orbit around the Sun is tilted relative to the ecliptic plane. The ecliptic plane is defined as the plane containing the Earth's orbit, so the Earth's inclination is 0.
orbital_inclination = {
    "Mercury": 7.0,
    "Venus": 3.4,
    "Earth": 5.1,
    "Mars": 1.8,
    "Jupiter": 1.3,
    "Saturn": 2.5,
    "Uranus": 0.8,
    "Neptune": 1.8
}

# Orbital Eccentricity - This is a measure of how far a planet's orbit about the Sun (or the Moon's orbit about the Earth) is from being circular. The larger the eccentricity, the more elongated is the orbit, an eccentricity of 0 means the orbit is a perfect circle. There are no units for eccentricity.
orbital_eccentricity = {
    "Mercury": 0.206,
    "Venus": 0.007,
    "Earth": 0.017,
    "Mars": 0.094,
    "Jupiter": 0.049,
    "Saturn": 0.052,
    "Uranus": 0.047,
    "Neptune": 0.010
}

# Obliquity to Orbit (degrees) - The angle in degrees the axis of a planet (the imaginary line running through the center of the planet from the north to south poles) is tilted relative to a line perpendicular to the planet's orbit around the Sun, north pole defined by right hand rule. Note: Venus rotates in a retrograde direction, opposite the other planets, so the tilt is almost 180 degrees, it is considered to be spinning with its "top", or north pole pointing "downward" (southward). Uranus rotates almost on its side relative to the orbit, Pluto is pointing slightly "down". The ratios with Earth refer to the axis without reference to north or south.
obliquity = {
    "Mercury": 0.034,
    "Venus": 177.4,
    "Earth": 23.4,
    "Mars": 25.2,
    "Jupiter": 3.1,
    "Saturn": 26.7,
    "Uranus": 97.8,
    "Neptune": 28.3
}

# Mean Temperature (C) - This is the average temperature over the whole planet's surface (or for the gas giants at the one bar level) in degrees C (Celsius or Centigrade). For Mercury and the Moon, for example, this is an average over the sunlit (very hot) and dark (very cold) hemispheres and so is not representative of any given region on the planet, and most of the surface is quite different from this average value. As with the Earth, there will tend to be variations in temperature from the equator to the poles, from the day to night sides, and seasonal changes on most of the planets.
mean_temp = {
    "Mercury": 167.0,
    "Venus": 464.0,
    "Earth": 15.0,
    "Mars": -65.0,
    "Jupiter": -110.0,
    "Saturn": -140.0,
    "Uranus": -195.0,
    "Neptune": -200.0
}

# Surface Pressure (bars) - This is the atmospheric pressure (the weight of the atmosphere per unit area) at the surface of the planet in bars. Note: The surfaces of Jupiter, Saturn, Uranus, and Neptune are deep in the atmosphere and the location and pressures are not known.
surface_pressure = {
    "Mercury": 0,
    "Venus": 92,
    "Earth": 1,
    "Mars": 0.01,
    "Jupiter": float("NaN"),
    "Saturn": float("NaN"),
    "Uranus": float("NaN"),
    "Neptune": float("NaN")
}

# Number of Moons - This gives the number of IAU officially confirmed moons orbiting the planet. New moons are still being discovered. 
num_moons = {
    "Mercury": 0,
    "Venus": 0,
    "Earth": 1,
    "Mars": 2,
    "Jupiter": 95,
    "Saturn": 146,
    "Uranus": 28,
    "Neptune": 16
}

# Ring System? - This tells whether a planet has a set of rings around it, Saturn being the most obvious example.
ring_system = {
    "Mercury": "No",
    "Venus": "No",
    "Earth": "No",
    "Mars": "No",
    "Jupiter": "Yes",
    "Saturn": "Yes",
    "Uranus": "Yes",
    "Neptune": "Yes"
}

# Global Magnetic Field? - This tells whether the planet has a measurable large-scale magnetic field. Mars and the Moon have localized regional magnetic fields but no global field.
global_magnetic_field = {
    "Mercury": "Yes",
    "Venus": "No",
    "Earth": "Yes",
    "Mars": "No",
    "Jupiter": "Yes",
    "Saturn": "Yes",
    "Uranus": "Yes",
    "Neptune": "Yes"
}
