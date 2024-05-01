import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

# Constants
G = 6.67430e-11  # gravitational constant
dt = 3600 * 24    # time step (1 day)
num_steps = 365   # number of steps

# Define initial conditions for celestial bodies (mass, position, velocity)
bodies = {
    "Sun": {"mass": 1.989e30, "pos": [0, 0, 0], "vel": [0, 0, 0], "color": (1, 1, 0)},  # yellow
    "Mercury": {"mass": 3.301e23, "pos": [0, 0.39e12, 0], "vel": [47000, 0, 0], "color": (0.7, 0.7, 0.7)},  # gray
    "Venus": {"mass": 4.867e24, "pos": [0.723e12, 0, 0], "vel": [0, 35000, 0], "color": (0.8, 0.4, 0)},  # orange
    "Earth": {"mass": 5.972e24, "pos": [1.0e12, 0, 0], "vel": [0, 29783, 0], "color": (0, 0.5, 1)},  # blue
    "Mars": {"mass": 6.39e23, "pos": [1.524e12, 0, 0], "vel": [0, 24100, 0], "color": (1, 0, 0)}  # red
}

# Function to simulate one step
def simulate_step():
    for name1, body1 in bodies.items():
        force = [0, 0, 0]
        for name2, body2 in bodies.items():
            if name1 != name2:
                r = [body2["pos"][i] - body1["pos"][i] for i in range(3)]
                dist = sum(x**2 for x in r) ** 0.5
                force = [force[i] + G * body1["mass"] * body2["mass"] * r[i] / dist**3 for i in range(3)]
        body1["vel"] = [body1["vel"][i] + force[i] / body1["mass"] * dt for i in range(3)]
        body1["pos"] = [body1["pos"][i] + body1["vel"][i] * dt for i in range(3)]

# Pygame initialization
pygame.init()
display = (800, 600)
pygame.display.set_mode(display, DOUBLEBUF | OPENGL)

# OpenGL initialization
gluPerspective(45, (display[0] / display[1]), 0.1, 10000.0)
glTranslatef(0.0, 0.0, -5000)

# Main loop
for _ in range(num_steps):
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            quit()

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    # Draw celestial bodies
    for name, body in bodies.items():
        glPushMatrix()
        glColor3fv(body["color"])
        glTranslatef(*body["pos"])
        glutSolidSphere(1000000, 20, 20)  # adjust the size of the spheres according to the scale
        glPopMatrix()

    pygame.display.flip()
    pygame.time.wait(10)  # adjust the speed of the simulation

    # Simulate one step
    simulate_step()
