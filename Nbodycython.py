import numpy as np
import array
import math 
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython import display
from nbody_cython import cython_Nbody_derivatives

def cython_derivs(pos,vel) :
    return cython_Nbody_derivatives(pos,vel,M)

def run_Nbody_rk2(tend, tframe, dt, p0, v0, derivatives):
    p, v = p0.copy(), v0.copy()
    t = 0
    tnext = tframe
    positions = [p.copy()]
    while t<tend :
        while t < tnext :
            delta_t = min(tnext-t,dt)
            dpdt, dvdt = derivatives(p,v) 
            phalf, vhalf = p+dpdt*0.5*delta_t, v+dvdt*0.5*delta_t
            dpdt, dvdt = derivatives(phalf, vhalf)
            p, v = p + dpdt*delta_t, v + dvdt*delta_t
            t += delta_t

        positions.append(p.copy())
        tnext += tframe
    return positions
    
def animate(i, positions):
    ax.clear()
    # Get the point from the points list at index i
    pos = positions[i]
    ax.scatter(pos[:,0], pos[:,1], color='green', marker='o')
    # Set the x and y axis to display a fixed range
    ax.set_xlim([-10,10])
    ax.set_ylim([-10,10])

positions = run_Nbody_rk2(frames*tframe, tframe, dt, p0, v0, cython_derivs)

# confirm it matches earlier results
fig, ax = plt.subplots()
ani = FuncAnimation(fig, lambda i : animate(i, positions), frames=frames, interval=50, repeat=False)
video = ani.to_html5_video()
ani.save('NbodyRK2Cython.gif')
html = display.HTML(video)
display.display(html)
plt.close()
