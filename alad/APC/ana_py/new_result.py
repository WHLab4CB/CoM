#!/usr/bin/env python3

### This script id used for collect results from run_traj

import numpy as np
import math
import copy
import os
import subprocess
import sys

### ----- setup parameters ------------
i_cycle = int(sys.argv[1])
n_cell = 64
n_point = 100
n_traj = 5
c_diff_cut = 0.1

if i_cycle < 2:
    cell_list = [i for i in range(0, n_cell)]
    print("Updata cell points:")
    print(cell_list)

else:
    diff_c_list = np.loadtxt("../../run_traj_"+str(int(i_cycle-1))+"/data0/aver_diff_c_value.txt")
    cell_list = []
    for i in range(len(diff_c_list)):
        if diff_c_list[i][1] > c_diff_cut:
            cell_list.append(int(diff_c_list[i][0])) 
    print("Updata cell points:")
    print(cell_list)

### ----------------- read the output of trajectory ----------------------------------
f_traj_out = open("./run_traj.out","w")
#for k in range(n_cell):
for k in cell_list:
    for i in range(n_point):
        for j in range(n_traj):
            name_log = "../cell_"+str(k)+"/point_"+str(i)+"/traj_"+str(j)+"/traj.out"  ###path to your file
        #name_log = "./step_6.log"  
            if os.path.exists(name_log):
                with open(name_log, 'r', encoding='utf-8') as f_in:
                    lines = f_in.readlines()        
                if lines:
                    if j == 0:
                        first_line = lines[0].strip().split()
                        print(first_line[2],first_line[3],end=" ",file=f_traj_out)
                        last_line = lines[-1].strip().split()
                        print(last_line[2],last_line[3],end=" ",file=f_traj_out)
                    else:
                        if j == (n_traj-1):
                            last_line = lines[-1].strip().split()
                            print(last_line[2],last_line[3],file=f_traj_out)
                        else:
                            last_line = lines[-1].strip().split()
                            print(last_line[2],last_line[3],end=" ",file=f_traj_out)

                else:
                    print("lines not found!")
            else:
                print("Warning:",name_log,"is not found!")

