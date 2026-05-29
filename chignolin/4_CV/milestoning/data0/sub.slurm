#!/usr/bin/bash

#SBATCH -J shooting
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=16

n_traj=10

for ((i=0; i<$n_traj; i++))
do
    cd ./traj_$i
    charmrun namd2 +auto-provision +isomalloc_sync shooting.namd >& shooting.log 
    cd ../
done

