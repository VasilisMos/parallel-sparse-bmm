#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --partition=batch
#SBATCH --time=00:01:00

OUTPUT=../logs/console_output.txt
num_iters=4

cd ..
# # Test Multithreaded
module purge
module load gcc/8.2.0

# make sequential
# make multithreaded
make clear_times

# iters=( $(seq 1 ${num_iters} ) )
iters=( 2 4 5 8 10)

./seq.out


for j in "${iters[@]}"
do  
    ./multithreaded.out $j
    echo $j
done