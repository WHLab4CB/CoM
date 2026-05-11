#!/usr/bin/env python3


import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import lil_matrix, diags
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys


def cal_near_ank(colvar,ank_id_list,colvar_ank_list):
    min_c = 999999
    min_v = -1

    n_ank = len(colvar_ank_list)
    for i in range(n_ank):
        dist_i = np.sqrt((colvar[0]-colvar_ank_list[i][0])**2+(colvar[1]-colvar_ank_list[i][1])**2)
        if i == 0:
            min_c = dist_i
            min_v = i
        else:
            if (dist_i < min_c):
                min_c = dist_i
                min_v = int(ank_id_list[i])
    return min_v

def load_data(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            parts = list(map(float, line.strip().split()))
            pairs = np.array([(parts[i],parts[i+1]) for i in range(0,len(parts),2)])
            data.append(pairs)
    return data

### input
data = load_data('../data/trajectories_all.txt')
n = len(data)
k = len(data[0])

first_pairs = np.array([line[0] for line in data])
states = np.concatenate([line[1:] for line in data], axis=0)
n_states = n*(k-1)

### calculate which near ank ------------------------------------
ank_list_0 = np.loadtxt("./ank.dat")
ank_id_list = ank_list_0[:,0]
colvar_ank_list = ank_list_0[:,1:]
n_ank = len(colvar_ank_list)

states_cell_list = []
for i in range(n_states):
    this_cell = cal_near_ank(states[i],ank_id_list,colvar_ank_list)
    states_cell_list.append(this_cell)
states_cell_list = np.array(states_cell_list)
### -------------------------------------------------------------

### calculate C_values ----------------------------------------------
kd_tree = cKDTree(first_pairs)
n_neighbors = 10
sigma = 0.1
#prob = 1/(n_neighbors*(k-1))

P = lil_matrix((n_states,n_states))

for s_idx in range(n_states):
    point = states[s_idx]
    distances, indices = kd_tree.query(point, k=n_neighbors)
    weight = []
    for j in range(n_neighbors):
        weight.append(np.exp(-(distances[j]/sigma)**2))
    tot_weight = sum(weight)
    counter = 0
    for j in indices:
        start = j*(k-1)
        end = start + (k-1)
        cols = np.arange(start,end)
        for col in cols:
            P[s_idx,col] += weight[counter]/tot_weight/(k-1)
            #P[s_idx,col] += prob
        counter += 1

P = P.tocsr()

reactant_center = np.array([-0.27,1.73])
product_center = np.array([0.84,0.00])
radius = 0.1

reactant_dist = np.linalg.norm(states - reactant_center, axis=1)
product_dist = np.linalg.norm(states - product_center, axis=1)

in_reactant = reactant_dist <= radius
in_product = product_dist <= radius
mask = np.logical_or(in_reactant, in_product)

diag_data = np.where(mask, 0.0, 1.0)
diag_mat = diags(diag_data)
P = diag_mat @ P

diag_data = np.where(in_product, 1.0, 0.0)
diag_mat = diags(diag_data)
P = diag_mat + P

ep = np.zeros(n_states)
ep[in_product] = 1.0

tolerance = 1e-8

C = P @ ep
while True:
    Cold = C.copy()
    C = P @ C
    #Cnorm = np.linalg.norm(C-Cold)
    #if Cnorm < tolerance:
    #    break

    max_diff = np.max(np.abs(C - Cold))
    if max_diff < tolerance:
        break

output_data = np.column_stack((states, C, states_cell_list))
np.savetxt('../data/C_values.dat', output_data, fmt='%.6f %.6f %.6e %.0f')

