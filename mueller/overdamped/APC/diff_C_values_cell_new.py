#!/usr/bin/env python3
### read points_cell from C_values

import numpy as np
import sys

### input
f_name_1 = sys.argv[1]
f_name_2 = sys.argv[2]

### output 
f_out = open("../data/aver_diff_c_value.txt","w")

###calculate which near ank
ank_list = np.loadtxt("./ank.dat")
n_ank = len(ank_list)
epsilon = 1e-15


### data pre-processing
cell_id_list = np.loadtxt(f_name_1,usecols=3)
c_value_1 = np.loadtxt(f_name_1,usecols=2)
n_value = len(c_value_1)
c_value_2 = np.loadtxt(f_name_2,usecols=2)
### check data
if len(c_value_2) < n_value:
    print(f"Warning: the number of {f_name_2}'s data is ({len(c_value_2)}), less than {f_name_1} ({n_value})")
    sys.exit()
c_value_2 = c_value_2[:n_value]

### calculate the percentage change in the C_value within each cell
c_value_percentage = [[] for _ in range(n_ank)]
for i in range(n_value):
    this_ank = int(cell_id_list[i])
    term_1 = np.log10(c_value_1[i]+epsilon) + np.log10(1.0-c_value_1[i]+epsilon)
    term_2 = np.log10(c_value_2[i]+epsilon) + np.log10(1.0-c_value_2[i]+epsilon)
    term = np.abs((term_2-term_1)/term_1)
    c_value_percentage[this_ank].append(term)
    

aver_c_value_percentage = np.zeros(n_ank)
for i in range(n_ank):
    print(i,len(c_value_percentage[i]))
    if len(c_value_percentage[i])==0:
        aver_c_value_percentage[i] = 0.0
    else:
        aver_c_value_percentage[i] = np.average(c_value_percentage[i])
    print(f"{i:.0f} {aver_c_value_percentage[i]:.6f}",file=f_out)


