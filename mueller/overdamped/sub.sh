#!/usr/bin/bash

#SBATCH -J APC
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=4



mkdir -p data
### test_0
echo " test_0"
mkdir -p test_0
cp ./APC/ank.dat ./data/xt_ini.txt
cd APC
python3 ./sampling_cell.py 
cp ../data/trajectories.txt ../data/trajectories_all.txt
python3 ./analog2.py
cd ..

cd plot
python3 ./plot_c_mueller.py
cd ..

cp ./data/C_values.dat ./test_0
cp ./data/xt_ini.txt ./test_0
cp ./data/xt_cell.txt ./test_0
cp ./data/trajectories.txt ./test_0
cp ./plot/c_mueller.png ./test_0


### test_1 
i=1
echo " test_$i"
let "i_1=$i-1"
mkdir -p ./test_$i
cp ./test_$i_1/xt_cell.txt ./data/xt_ini.txt

cd APC
python3 ./sampling_cell.py 
cat ../test_$i_1/trajectories.txt ../data/trajectories.txt > ../data/trajectories_all.txt
python3 ./analog2.py
python3 ./diff_C_values_cell_new.py ../test_$i_1/C_values.dat ../data/C_values.dat
cp ../data/aver_diff_c_value.txt ../test_$i/
cd ..

cd plot
python3 ./plot_c_mueller.py
cd ..

cp ./data/C_values.dat ./test_$i
cp ./data/xt_ini.txt ./test_$i
cp ./data/xt_cell.txt ./test_$i
cp ./data/trajectories_all.txt ./test_$i/trajectories.txt
cp ./plot/c_mueller.png ./test_$i

for ((i=2; i<10; i++))
do
    echo " test_$i"
    let "i_1=$i-1"
    cp ./test_$i_1/xt_cell.txt ./data/xt_ini.txt

    cd APC
    python3 ./sampling_cell_diff.py 
    python_output=$(python ./sampling_cell_diff.py)
    exit_code=$?
    if [ $exit_code -eq 1 ]; then
        cd ..
        break
    else
        cat ../test_$i_1/trajectories.txt ../data/trajectories.txt > ../data/trajectories_all.txt
        python3 ./analog2.py
        python3 ./diff_C_values_cell_new.py ../test_$i_1/C_values.dat ../data/C_values.dat
        cd ..

        cd plot
        python3 ./plot_c_mueller.py
        cd ..
        
        mkdir -p ./test_$i
        
        cp ./data/aver_diff_c_value.txt ./test_$i/
        cp ./data/C_values.dat ./test_$i
        cp ./data/xt_ini.txt ./test_$i
        cp ./data/xt_cell.txt ./test_$i
        cp ./data/trajectories_all.txt ./test_$i/trajectories.txt
        cp ./plot/c_mueller.png ./test_$i
    fi
done

cd APC
    python3 ./nn.py
cd ..

cd plot
    python3 ./plot_iso_c.py
cd ..



