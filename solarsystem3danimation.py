from vpython import*
scene = canvas()

G = 6.67e-11
AU = 1.496e+11

sun = sphere(pos=vector(0,0,0), radius=AU/20, color=color.yellow, m=1988500e24)
mercury = sphere(pos=vector(0.387*AU,0,0), radius=AU/40, make_trail=True, color=color.brown)
mercury.v = vector(0,47.36e3,0)
venus = sphere(pos=vector(0.72*AU,0,0), radius=AU/40, make_trail=True, color=color.orange)
venus.v = vector(0,35.02e3,0)
earth = sphere(pos=vector(AU,0,0), radius=AU/2, make_trail=True, texture=textures.earth)
earth.v = vector(0,29.78e3,0)

planets=[mercury,venus,earth]
t = 0
dt = 3600

while t < 3600*20*100:
    rate(1000)
    for p in planets:
        p.a = -G*sun.m*norm(p.pos)/mag(p.pos)**2
        p.v = p.v + p.a*dt
        p.pos = p.pos + p.v*dt
    t = t + dt
    