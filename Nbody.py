import numpy as np
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

# some simple initial conditions
def initial_conditions(): 
    pos_and_vel = np.zeros([N_bodies,6])

    pos_and_vel[0,0] = -1
    pos_and_vel[0,4] = -0.25
    pos_and_vel[1,0] = 1
    pos_and_vel[1,4] = 0.25
    return pos_and_vel

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

N_bodies = 2
M = np.array([1.0, 0.25])

frames = 100
tframe = 0.25
dt = 0.025

# Call the RK2 integrator
rk2positions = run_Nbody_rk2(frames*tframe, tframe, dt)

# Call the leapfrog integrator
lfpositions = run_Nbody_leapfrog(frames * tframe, tframe, dt)

# Use the RK2 function in FuncAnimation
figrk2, axrk2 = plt.subplots()
anirk2 = FuncAnimation(figrk2, lambda i : rk2animate(i, rk2positions), frames=len(rk2positions), interval=50, repeat=False)
video = anirk2.to_html5_video()
anirk2.save('NbodyRK2.gif')
html = display.HTML(video)
display.display(html)
#plt.close()

# Use the RK2 function in FuncAnimation
for pos in rk2positions:
    axrk2.scatter(pos[0,0], pos[0,1], color='green', marker='o')
    axrk2.scatter(pos[1,0], pos[1,1], color='red', marker='o')

plt.xlabel('x')
plt.ylabel('y')
plt.show()

# Use the lfanimate function in FuncAnimation
figlf, axlf = plt.subplots()
anilf = FuncAnimation(figlf, lambda i : lfanimate(i, lfpositions), frames=len(lfpositions), interval=50, repeat=False)
video = anilf.to_html5_video()
anilf.save('NbodyLF.gif')
html = display.HTML(video)
display.display(html)
#plt.close()

# Plot the positions for leapfrog integration
for pos in lfpositions:
    axlf.scatter(pos[0,0], pos[0,1], color='green', marker='o')
    axlf.scatter(pos[1,0], pos[1,1], color='red', marker='o')

plt.xlabel('x')
plt.ylabel('y')
plt.show()

plt.close()
