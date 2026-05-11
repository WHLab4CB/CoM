#!/usr/bin/env python3


import numpy as np
from committor_numpy import NumpyCommittor

muller_params = {
    'A': [-200, -100, -170, 15],
    'a': [-1, -1, -6.5, 0.7],
    'b': [0, 0, 11, 0.6],
    'c': [-10, -10, -6.5, 0.7],
    'x0': [1, 0, -0.5, -1],
    'y0': [0, 0.5, 1.5, 1]
}

def muller_potential(position):
    x,y = position
    V = 0.0
    for i in range(4):
        dx = x - muller_params['x0'][i]
        dy = y - muller_params['y0'][i]
        V += muller_params['A'][i]*np.exp(muller_params['a'][i]*dx**2+ muller_params['b'][i]*dx*dy + muller_params['c'][i]*dy**2)

    return V

def muller_gradient(position):
    x,y = position
    grad = np.zeros(2)
    for i in range(4):
        dx = x - muller_params['x0'][i]
        dy = y - muller_params['y0'][i]
        exp_term = np.exp(muller_params['a'][i]*dx**2 + muller_params['b'][i]*dx*dy + muller_params['c'][i]*dy**2)

        dVdx = muller_params['A'][i]*exp_term * (2*muller_params['a'][i]*dx + muller_params['b'][i]*dy)

        dVdy = muller_params['A'][i]*exp_term * (muller_params['b'][i]*dx + 2*muller_params['c'][i]*dy)

        grad += np.array([dVdx, dVdy])
    return grad

def exit_cond(index,x):
    hit = 0
    q_curr = committor.forward(x)
    if index==0:
        if q_curr>q_series[index+1]:
            hit = 1
        elif q_curr<epsilon:
            hit = 2
    elif index==q_len-1:
        if q_curr>1-epsilon:
            hit = 1
        elif q_curr<q_series[index-1]:
            hit = 2
    else:
        if q_curr>q_series[index+1]:
            hit = 1
        elif q_curr<q_series[index-1]:
            hit = 2

    return hit

q_series = [5e-4,1e-3,1e-2,3e-2,0.1,0.5,0.9,0.95,0.98,0.995]
q_len = len(q_series)
temp = 10.0
mass = 1.0
gamma = 10.0
dt = 1e-4
epsilon = 1e-12
committor = NumpyCommittor("../data/network_params.npz")

Kmat = np.zeros((q_len+2, q_len+2))
Kmat[0,1] = 1.
Kmat[q_len+1,0] = 1.
tave = []

fname = f"../data/reactant_samples.dat"
coords = np.loadtxt(fname,usecols=(0,1))
#vels = np.loadtxt(fname,usecols=(2,3))
num_coords = len(coords)
tlist = []
for j in range(num_coords):
    x = coords[j]
    v = np.sqrt(temp/mass)*np.random.randn(2)
    #v = vels[j]
    nstep = 0
    q_curr = committor.forward(x)
    if q_curr > q_series[0]:
        continue
    while True:
        force = -muller_gradient(x)
        x += v*dt
        v = v - gamma*dt/mass*v + force/mass*dt + np.sqrt(2*temp*gamma*dt)/mass*np.random.randn(2)
        nstep += 1
        q_curr = committor.forward(x)
        if q_curr > q_series[0]:
            tlist.append(nstep*dt)
            break
tave.append(np.average(tlist))


for i in range(q_len):
    fname = f"../data/committor_{q_series[i]:.2e}_samples.dat"
    coords = np.loadtxt(fname,usecols=(0,1))
    #vels = np.loadtxt(fname,usecols=(2,3))
    num_coords = len(coords)
    tlist = []
    for j in range(num_coords):
        x = coords[j]
        v = np.sqrt(temp/mass)*np.random.randn(2)
        #v= vels[j]
        nstep = 0
        hit = exit_cond(i,x)
        if hit>0:
            continue
        while True:
            force = -muller_gradient(x)
            x += v*dt
            v = v - gamma*dt/mass*v + force/mass*dt + np.sqrt(2*temp*gamma*dt)/mass*np.random.randn(2)
            nstep += 1
            hit = exit_cond(i,x)
            if hit > 0:
                if hit==1:
                    Kmat[i+1,i+2] += 1
                else:
                    Kmat[i+1,i] += 1
                tlist.append(nstep*dt)
                break
    tave.append(np.average(tlist))

fname = f"../data/product_samples.dat"
coords = np.loadtxt(fname,usecols=(0,1))
#vels = np.loadtxt(fname,usecols=(2,3))
num_coords = len(coords)
tlist = []
for j in range(num_coords):
    x = coords[j]
    v = np.sqrt(temp/mass)*np.random.randn(2)
    #v = vels[j] 
    nstep = 0
    q_curr = committor.forward(x)
    if q_curr < q_series[-1]:
        continue
    while True:
        force = -muller_gradient(x)
        x += v*dt
        v = v - gamma*dt/mass*v + force/mass*dt + np.sqrt(2*temp*gamma*dt)/mass*np.random.randn(2)
        nstep += 1
        q_curr = committor.forward(x)
        if q_curr < q_series[-1]:
            tlist.append(nstep*dt)
            break
tave.append(np.average(tlist))

for i in range(q_len+2):
    if i==0 or i==q_len+1:
        continue
    rowsum = sum(Kmat[i,:])
    for j in range(q_len+2):
        Kmat[i,j] = Kmat[i,j]/rowsum


fout1 = f"../data/lifetime.dat"
np.savetxt(fout1, tave, fmt='%.2f')
fout2 = f"../data/Kmat.dat.simulate"
np.savetxt(fout2, Kmat, fmt='%.2f')
