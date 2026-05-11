#!/usr/bin/bash

#SBATCH -J APC
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=2

current_dir="CoM/alad/APC"
#echo "$current_dir"

n_ank=64
i=0
#j=0
n_points=100
n_traj=5


for ((j=0; j<$n_ank; j++))
do 
    mkdir -p $current_dir/sample_$i/cell_$j/
    cp $current_dir/sample_data/sample_cell.namd $current_dir/sample_$i/cell_$j/
    cp $current_dir/sample_data/Aladipep.pdb $current_dir/sample_$i/cell_$j/
    cp $current_dir/sample_data/Aladipep.psf $current_dir/sample_$i/cell_$j/
    cp $current_dir/sample_data/par_all22_prot_nocmap.prm $current_dir/sample_$i/cell_$j/
    cp $current_dir/sample_data/colvar_sample_cell.conf $current_dir/sample_$i/cell_$j/
    cp $current_dir/sample_data/ank.txt $current_dir/sample_$i/cell_$j/
    if [ $i -eq 0 ]; then
        cp $current_dir/sample_ini/cell_$j/ini.coor $current_dir/sample_$i/cell_$j/
        cp $current_dir/sample_ini/cell_$j/sample_cell.tcl $current_dir/sample_$i/cell_$j/
    else
        let "i_1=i-1"
        cp $current_dir/sample_${i_1}/cell_$j/sample_cell.coor $current_dir/sample_$i/cell_$j/ini.coor
        cp $current_dir/sample_ini/cell_$j/sample_cell.tcl $current_dir/sample_$i/cell_$j/
    fi

    cd $current_dir/sample_$i/cell_$j/
    charmrun namd2 +auto-provision +isomalloc_sync sample_cell.namd >& sample_cell.log
    cat sample_cell.log | grep final_result > sample_cell.out

    ### run_traj
    for ((k=0; k<$n_points; k++))
    do
        let "sample_value_int=k*400+400"

        mkdir -p $current_dir/run_traj_$i/cell_$j/point_$k/

        cp "$current_dir/sample_$i/cell_$j/sample_cell.${sample_value_int}.coor" "$current_dir/run_traj_$i/cell_$j/point_$k/sample_cell.coor" 
        cp $current_dir/run_traj_data/Aladipep.pdb $current_dir/run_traj_$i/cell_$j/point_$k
        cp $current_dir/run_traj_data/Aladipep.psf $current_dir/run_traj_$i/cell_$j/point_$k
        cp $current_dir/run_traj_data/par_all22_prot_nocmap.prm $current_dir/run_traj_$i/cell_$j/point_$k
        cp $current_dir/run_traj_data/colvar_run.conf $current_dir/run_traj_$i/cell_$j/point_$k
        cp $current_dir/run_traj_data/run_traj.tcl $current_dir/run_traj_$i/cell_$j/point_$k
    
        for ((l=0; l<$n_traj; l++))
        do
            mkdir -p $current_dir/run_traj_$i/cell_$j/point_$k/traj_$l
            cp $current_dir/run_traj_data/run_traj.namd $current_dir/run_traj_$i/cell_$j/point_$k/traj_$l
        done
    done

    for ((l=0; l<$n_traj; l++))
    do
        for ((k=0; k<$n_points; k++))
        do
            cd $current_dir/run_traj_$i/cell_$j/point_$k/traj_$l/
            charmrun namd2 +auto-provision +isomalloc_sync run_traj.namd >& run_traj.log
            cat run_traj.log | grep final_result > traj.out
        done
    done
done

###analysis
mkdir -p $current_dir/run_traj_$i/ana_data
cp $current_dir/ana_py/ank.txt $current_dir/run_traj_$i/ana_data
cd $current_dir/run_traj_$i/ana_data
python3 ../../ana_py/new_result.py $i
if [ $i -eq 0 ]; then
    cp ./run_traj.out ./run_traj_all.out
else
    cat ../../run_traj_$i_1/ana_data/run_traj_all.out ./run_traj.out > ./run_traj_all.out
fi
python3 ../../ana_py/analog2.py ./run_traj_all.out
if [ $i -ne 0 ]; then
    python3 ../../ana_py/diff_C_values_cell_new.py ../../run_traj_$i_1/ana_data/C_values.dat ./C_values.dat
fi


