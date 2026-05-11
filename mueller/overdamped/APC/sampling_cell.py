#!/usr/bin/env python3
### w
### 2025.10.29

import numpy as np
import math
import random
import struct
from typing import List, Tuple
import sys

### output
f_sampling = open("../data/samples_cell.txt","w")
f_xt_cell = open("../data/xt_cell.txt","w")

### params
ank_list = np.loadtxt("./ank.dat")
xt_list = np.loadtxt("../data/xt_ini.txt",usecols=(0,1,2))


n_ank = len(ank_list)
ank_id = ank_list[:,0]
x_ank = ank_list[:,1:] 

kBT_sam = 20.0
kBT = 10.0
gamma = 10
dt_sampling = 5e-4
k_con = 8000.0
dt = 1e-3
save_interval = 200
target_samples = 100
n_steps = 10
n_trajectories = 5

reactant_center = np.array([-0.27, 1.73])
product_center = np.array([0.84, 0.00])
radius = 0.1

sampling_x_range = (-1.7, 1.2)
sampling_y_range = (-0.4, 2.1)


terms = [
    {'A': -200.0, 'a': -1.0, 'b': 0.0, 'c': -10.0, 'x0': 1.0, 'y0': 0.0},
    {'A': -100.0, 'a': -1.0, 'b': 0.0, 'c': -10.0, 'x0': 0.0, 'y0': 0.5},
    {'A': -170.0, 'a': -6.5, 'b': 11.0, 'c': -6.5, 'x0': -0.5, 'y0': 1.5},
    {'A':   15.0, 'a':  0.7, 'b': 0.6, 'c':  0.7, 'x0': -1.0, 'y0': 1.0}
]

def compute_potential(x, y, terms):
    V = 0.0
    for term in terms:
        dx = x - term['x0']
        dy = y - term['y0']
        exponent = term['a'] * dx**2 + term['b'] * dx * dy + term['c'] * dy**2
        V += term['A'] * np.exp(exponent)
    return V

def compute_gradient(x, y, terms):
    grad_x, grad_y = 0.0, 0.0
    for term in terms:
        dx = x - term['x0']
        dy = y - term['y0']
        exponent = term['a'] * dx**2 + term['b'] * dx * dy + term['c'] * dy**2
        prefactor = term['A'] * np.exp(exponent)
        grad_x += prefactor * (2 * term['a'] * dx + term['b'] * dy)
        grad_y += prefactor * (term['b'] * dx + 2 * term['c'] * dy)
    return grad_x, grad_y

def cell_bias_force(x, x_ank, ank1,k_con):
    ### Calculate bias force for sampling in cells
    n_ank = len(x_ank)
    # Calculate distance between this point and initial anchor
    x_diff0 = x[0] - x_ank[ank1, 0]
    x_diff1 = x[1] - x_ank[ank1, 1]
    dx1 = math.sqrt(x_diff0**2 + x_diff1**2)
    
    # Calculate bias force and harmonic restraints
    u2x_bias = [0.0, 0.0]
    
    # Harmonic wall restraints
    for i in range(n_ank):
        if i == ank1:
            continue
            
        x_diff2 = x[0] - x_ank[i, 0]
        x_diff3 = x[1] - x_ank[i, 1]
        dx3 = math.sqrt(x_diff2**2 + x_diff3**2)
        dx2 = dx1 - dx3

        if dx2 > 0.0:
            u2x_bias[0] += 2.0 * k_con * dx2 * (x_diff0 / dx1 - x_diff2 / dx3)
            u2x_bias[1] += 2.0 * k_con * dx2 * (x_diff1 / dx1 - x_diff3 / dx3)
    return u2x_bias

def in_reactant(x, y):
    dx = x - reactant_center[0]
    dy = y - reactant_center[1]
    return dx**2 + dy**2 <= radius**2

def in_product(x, y):
    dx = x - product_center[0]
    dy = y - product_center[1]
    return dx**2 + dy**2 <= radius**2

def run_sampling_cell(x_ank, initial_pos, this_ank, kBT_sam, dt, gamma, save_interval, target_samples, x_boundary, y_boundary):
    current_pos = np.array(initial_pos, dtype=np.float64)
    noise_std = np.sqrt(2 * kBT_sam * dt / gamma)
    x_sampling = []
    for step in range(200000):
        grad_x, grad_y = compute_gradient(*current_pos, terms)
        u2x_bias = cell_bias_force(current_pos, x_ank, this_ank,k_con)
        xt_old = current_pos.copy()

        current_pos[0] += (-(grad_x + u2x_bias[0]) / gamma) * dt + np.random.normal(0, noise_std)
        current_pos[1] += (-(grad_y + u2x_bias[1]) / gamma) * dt + np.random.normal(0, noise_std)
        
        # Check boundaries
        if (current_pos[0] <= x_boundary[0] or current_pos[0] >= x_boundary[1] or 
            current_pos[1] <= y_boundary[0] or current_pos[1] >= y_boundary[1]):
            current_pos = xt_old.copy()
            
        if (step+1) % save_interval == 0:
            #current_ank = which_near(x_ank, current_pos)
            #if current_ank == this_ank:
            #    x_sampling.append(current_pos.copy())
            x_sampling.append(current_pos.copy())
        if len(x_sampling) >= target_samples:
            break
    return x_sampling

def run_trajectory(initial_pos, kBT, gamma, dt, n_steps, terms):
    current_pos = np.array(initial_pos, dtype=np.float64)
    if in_reactant(*current_pos) or in_product(*current_pos):
        return current_pos.tolist()
    
    noise_std = np.sqrt(2 * kBT * dt / gamma)
    for _ in range(n_steps):
        grad_x, grad_y = compute_gradient(*current_pos, terms)
        current_pos[0] += (-grad_x / gamma) * dt + np.random.normal(0, noise_std)
        current_pos[1] += (-grad_y / gamma) * dt + np.random.normal(0, noise_std)
        if in_reactant(*current_pos) or in_product(*current_pos):
            break
    return current_pos.tolist()

def which_near(x_ank, xt):
    ### Determine which anchor is nearest
    n_ank = len(x_ank)
    near_ank = -9999
    dr_min = 0.0
    
    for i in range(n_ank):
        dx2 = math.sqrt((xt[0] - x_ank[i, 0])**2 + (xt[1] - x_ank[i, 1])**2)
        if i == 0:
            dr_min = dx2
            near_ank = i
        else:
            if dx2 < dr_min:
                dr_min = dx2
                near_ank = i              
    return int(near_ank)

### run
initial_configs = []
for i in range(n_ank):
    ini_ank1 = int(xt_list[i][0])
    xt = xt_list[i][1:]
    this_ank = which_near(x_ank, xt)
    if this_ank != ini_ank1:
        print("Warning: initial position is", this_ank, ", not in target anchor",ini_ank1,"!")
    
    sampling_list = run_sampling_cell(x_ank, xt, ini_ank1, kBT_sam, dt_sampling, gamma, save_interval, target_samples, sampling_x_range, sampling_y_range)
    this_ank_1 = which_near(x_ank, sampling_list[-1])
    #print(f'{this_ank_1:.0f} {sampling_list[-1][0]:.6f} {sampling_list[-1][1]:.6f}',file=f_xt_cell)
    print(f'{ini_ank1:.0f} {sampling_list[-1][0]:.6f} {sampling_list[-1][1]:.6f}',file=f_xt_cell)
   

    ### without E_cut 
    for j in range(len(sampling_list)):
        #if not in_reactant(sampling_list[j][0], sampling_list[j][1]) and not in_product(sampling_list[j][0], sampling_list[j][1]):
        initial_configs.append((sampling_list[j][0], sampling_list[j][1]))
        this_ank_2 = which_near(x_ank, sampling_list[j])
        print(f'{this_ank:.0f} {sampling_list[j][0]:.6f} {sampling_list[j][1]:.6f} {this_ank_2:.0f} ',file=f_sampling)

with open('../data/trajectories.txt', 'w') as f:
    for init in initial_configs:
        endpoints = []
        for _ in range(n_trajectories):
            endpoint = run_trajectory(init, kBT, gamma, dt, n_steps, terms)
            endpoints.extend(endpoint)
        line = [f"{coord:.6f}" for coord in list(init) + endpoints]
        f.write(' '.join(line) + '\n')


