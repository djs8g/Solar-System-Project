import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython import display 

def Nbody_derivatives(pos, vel):
    dpdt = vel
    dvdt = np.zeros(vel.shape)
    for i in range(N_bodies):
        for j in range(N_bodies):
            if i == j: 
                continue
            r = np.linalg.norm(pos[j]-pos[i])
            mass = M[j]
            rhat = (pos[j] - pos[i])/r
            dvdt[i] -= mass/(r*r)*rhat
        
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

def initial_conditions() : 
    pos = np.random.random([N_bodies,3])
    vel = np.random.random([N_bodies,3])

    return pos, vel
    
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

N_massive = 20
N_bodies = N_massive
M = np.ones(N_bodies)

tframe = 0.01
dt = 0.001
frames = 100

# use the same initial conditions for all future runs
p0, v0 = initial_conditions()

positions = run_Nbody_rk2(frames*tframe, tframe, dt, p0, v0, Nbody_derivatives)
positions2 = run_Nbody_rk2(frames*tframe, tframe, dt, p0, v0, Nbody_derivatives2)
positions3 = run_Nbody_rk2(frames*tframe, tframe, dt, p0, v0, Nbody_derivatives3)

fig, ax = plt.subplots()
ani = FuncAnimation(fig, lambda i : animate(i, positions), frames=len(positions), interval=50, repeat=False)
video = ani.to_html5_video()
ani.save('NbodyRK2-0.gif')
html = display.HTML(video)
display.display(html)
#plt.close()

fig, ax = plt.subplots()
ani = FuncAnimation(fig, lambda i : animate(i, positions2), frames=len(positions2), interval=50, repeat=False)
video = ani.to_html5_video()
ani.save('NbodyRK2-1.gif')
html = display.HTML(video)
display.display(html)
#plt.close()

fig, ax = plt.subplots()
ani = FuncAnimation(fig, lambda i : animate(i, positions3), frames=frames, interval=50, repeat=False)
video = ani.to_html5_video()
ani.save('NbodyRK2-2.gif')
html = display.HTML(video)
display.display(html)
plt.close()
