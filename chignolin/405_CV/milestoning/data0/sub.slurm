#!/usr/bin/bash

#SBATCH -J shooting
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=16


charmrun namd2 +auto-provision +isomalloc_sync shooting.namd >& shooting.log 

