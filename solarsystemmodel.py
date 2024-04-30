G = 6.67e-11
AU = 1.496e+11

sun = sphere(pos=vector(0,0,0), radius=AU/20, color=color.yellow, m=1988500e24)
mercury = sphere(pos=vector(57.909e9,0,0), radius=AU/40, make_trail=True)
mercury.v = vector(0,47.36e3,0)

t = 0
dt = 3600

while t < 3600*20*100:
    rate(1000)
    mercury.a = -G*sun.m*norm(mercury.pos)/mag(mercury.pos)**2
    mercury.v = mercury.v + mercury.a*dt
    mnercury.pos = mercury.pos + mercury.v*dt
    t = t + dt
