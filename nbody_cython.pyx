import numpy as np
import array
import math 
from cpython cimport array

def cython_Nbody_derivatives(pos, vel, M):
    cdef int N_bodies = M.size
    cdef int i = 0, j = 0 
    cdef double [:,:] p = pos
    cdef double [:] mass = M
    rhat = np.zeros(3)
    cdef double [:] rh = rhat
    cdef double r = 0., r2 = 0.
    dpdt = vel
    dvdt_arr = np.zeros(vel.shape)
    cdef double [:,:] dvdt = dvdt_arr
    for i in range(N_bodies):
        for j in range(N_bodies):
            if i == j: 
                continue
            r2 = 0.
            for k in range(3) : 
                r2 += p[j,k]*p[j,k] + p[i,k]*p[i,k]
            r = math.sqrt(r2)
            for k in range(3) : 
                rh[k] = (p[j,k] - p[i,k])/r
                dvdt[i,k] += -mass[j]/(r*r)*rh[k]
        
    return dpdt, dvdt_arr
