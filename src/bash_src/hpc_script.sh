#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=batch
#SBATCH --time=00:01:00

OUTPUT=../logs/console_output.txt
num_iters=4

# Test Multithreaded
module purge
module load gcc/8.2.0

make sequential >> OUTPUT
make multithreaded >> OUTPUT
make clear_times >> OUTPUT

for j in "${iters[@]}"
do  
    ./seq.out >> OUTPUT
    ./multithreaded.out 8  >> $OUTPUT
done