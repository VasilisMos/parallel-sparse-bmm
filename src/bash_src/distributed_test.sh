#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --partition=batch
#SBATCH --time=00:02:00

OUTPUT=./logs/console_output.txt
num_iters=4
dec_exponent=$((10 ** 5))
MPI_PROCS=8

cd ..
module purge
module load gcc/9.2.0 openmpi/3.1.5  # version capable of being used alongside CUDA
module load octave

dims=( 2 3 4 5 6 8 )
iters=( $(seq 1 ${num_iters} ) )

make clear_times >> $OUTPUT
make distributed >> $OUTPUT
 
for n in "${dims[@]}"
do
    dim_real=$(($n * $dec_exponent))
    make data N=$dim_real >> $OUTPUT
    for j in "${iters[@]}"
    do  
        fast_test_distributed MPI_PROCS=$MPI_PROCS  >> $OUTPUT
    done
done

make clean >> $OUTPUT



