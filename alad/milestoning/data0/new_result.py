#!/usr/bin/env python3

### This script id used for collect results from cm

import numpy as np
import math
import copy
import os
import subprocess

### ----- setup parameters ------------
n_mile = 12
n_traj = 10


### ----------------- read the output of trajectory forward in time ----------------------------------
f_forward_out = open("./cm.out","w")
for i in range(n_mile):
    for j in range(n_traj):
        name_log = "../mile_"+str(i)+"/traj_"+str(j)+"/shooting.log"  ###path to your file
        if os.path.exists(name_log):
            proc_temp = subprocess.Popen('cat ' + name_log + ' | grep final_result',stdout=subprocess.PIPE,shell=True)
            line_temp = proc_temp.stdout.readlines()
            if len(line_temp) > 0:
                line_list = (line_temp[0].decode('utf_8')).split()
                #print("line_list:",line_list)
                print(line_list[2],line_list[3],line_list[4],line_list[5],j,file=f_forward_out)
            else:
                print("Warning: final_result not found in mile_",i,"/traj_",j,"/shooting.log!",sep="")
        else:
            print("Warning:",name_log,"is not found!")
f_forward_out.close()
