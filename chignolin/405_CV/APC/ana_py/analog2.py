#!/usr/bin/env python3


import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import lil_matrix, diags
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import copy
import sys


def cal_near_ank(colvar,ank_id_list,colvar_ank_list):
    min_c = 999999
    min_v = -1

    n_ank = len(colvar_ank_list)
    for i in range(n_ank):
        dist_i = np.abs(colvar-colvar_ank_list[i])
        if i == 0:
            min_c = dist_i
            min_v = int(ank_id_list[i])
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
            #pairs = np.array([(parts[i],parts[i+1]) for i in range(0,len(parts),2)])
            #pairs = np.array([(parts[i],parts[i+1],parts[i+2],parts[i+3],parts[i+4]) for i in range(0,len(parts),5)])
            pairs = np.array([(*parts[i:i+406],) for i in range(0, len(parts), 406)])
            data.append(pairs)
    return data


### input
f_traj_name = sys.argv[1]
data = load_data(f_traj_name)
n = len(data)
k = len(data[0])

first_pairs = np.array([line[0] for line in data])
states = np.concatenate([line[1:] for line in data], axis=0)

states_coords = states[:,:405] 
mask_values = states[:,405] 
n_states = n*(k-1)
hbond1 = states[:,74]
hbond2 = states[:,78]
hbond3 = states[:,191]
hbond4 = states[:,195]

### calculate which near ank ------------------------------------
ank_list_0 = np.loadtxt("./ank.txt")
ank_id_list = ank_list_0[:,0]
colvar_ank_list = ank_list_0[:,1]
n_ank = len(colvar_ank_list)

states_cell_list = []
for i in range(n_states):
    this_cell = cal_near_ank(mask_values[i],ank_id_list,colvar_ank_list)
    states_cell_list.append(this_cell)
states_cell_list = np.array(states_cell_list)
### -------------------------------------------------------------

### calculate C_values ----------------------------------------------
kd_tree = cKDTree(first_pairs[:,:405])
n_neighbors = 50
sigma = 4
reactant_center = 15.0
product_center = 40.0

P = lil_matrix((n_states,n_states))

for s_idx in range(n_states):
    point = states_coords[s_idx]
    distances, indices = kd_tree.query(point, k=n_neighbors)
    weight = []
    for j in range(n_neighbors):
        weight.append(np.exp(-(distances[j]/sigma)**2))
    tot_weight = sum(weight)
    if tot_weight == 0:
        print("Warning:",s_idx,distances,weight)
    counter = 0
    for j in indices:
        start = j*(k-1)
        end = start + (k-1)
        cols = np.arange(start,end)
        for col in cols:
            P[s_idx,col] += weight[counter]/tot_weight/(k-1)
        counter += 1

P = P.tocsr()

in_reactant = mask_values <= reactant_center
in_product =  mask_values >= product_center
mask = np.logical_or(in_reactant, in_product)

diag_data = np.where(mask, 0.0, 1.0)
diag_mat = diags(diag_data)
P = diag_mat @ P

diag_data = np.where(in_product, 1.0, 0.0)
diag_mat = diags(diag_data)
P = diag_mat + P

ep = np.zeros(n_states)
ep[in_product] = 1.0

tolerance = 1e-6

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

output_data = np.column_stack((hbond1,hbond2,hbond3,hbond4,mask_values,C,states_cell_list))
np.savetxt('./C_h4.dat', output_data, fmt='%.6f %.6f %.6f %.6f %.6f %.6e %.0f')


output_data = np.column_stack((states,C,states_cell_list))
states_fmt = ' '.join(['%.6f'] * 406)   
c_fmt = '%.6e'                          
cell_fmt = '%.0f'                      

fmt_str = f"{states_fmt} {c_fmt} {cell_fmt}"


np.savetxt('./C_values.dat', output_data, fmt=fmt_str)



