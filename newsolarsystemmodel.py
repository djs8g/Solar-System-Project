import turtle, math, random

# Configure the variables
sizePlanet = [4, 0.4, 0.9, 1.0, 0.6, 2.5, 2.0, 1.6, 1.5]
sizeOrbit = [2, 60, 80, 100, 120, 170, 230, 300, 500]
orbitSpeed = 5000
planetColours = [['yellow', 'orange'],['grey', 'white'],['orange', 'yellow'],['blue', 'green'],['red', 'orange'],['orange', 'brown'],['yellow', 'white'],['turquoise', 'blue'],['blue', 'turquoise']]

# Configure the background
screen = turtle.Screen()
screen.bgcolor('black')

Sun = turtle.Turtle()
Mercury = turtle.Turtle()
Venus = turtle.Turtle()
Earth = turtle.Turtle()
Mars = turtle.Turtle()
Jupiter = turtle.Turtle()
Saturn = turtle.Turtle()
Uranus = turtle.Turtle()
Neptune = turtle.Turtle()

# Configure the planet parameters
Bodies = [Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune]

for i in range(len(Bodies)):
    Bodies[i].speed(orbitSpeed/sizeOrbit[i])
    Bodies[i].shape('circle')
    Bodies[i].shapesize(stretch_wid=sizePlanet[i])
    Bodies[i].fillcolor(planetColours[i][0])
    Bodies[i].pencolor(planetColours[i][1])

    # Configure the planet positions
    Bodies[i].left(random.randint(0,360))
    Bodies[i].up()
    Bodies[i].forward(sizeOrbit[i])
    Bodies[i].down()
    Bodies[i].left(90)

# Orbit the planets
while True:
    for i in range(len(Bodies)):
        if i != 0:
            Bodies[i].forward((sizeOrbit[i]*math.pi/0.4)/sizeOrbit[i])
            Bodies[i].left(450/sizeOrbit[i])


turtle.done()