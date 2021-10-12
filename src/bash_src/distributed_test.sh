#!/bin/bash
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --partition=batch
#SBATCH --time=00:05:00

OUTPUT=../logs/console_output.txt
num_iters=4
dec_exponent=$((12 * 10 ** 4))
MPI_PROCS=5

cd ..
module purge
module load gcc/8.2.0 openmpi  # version capable of being used alongside CUDA
module load octave

make distributed >> $OUTPUT
make clear_times >> $OUTPUT

dims=( 2 3 4 5 6 8 )
iters=( $(seq 1 ${num_iters} ) )

for n in "${dims[@]}"
do
    dim_real=$(($n * $dec_exponent))
    echo $dim_real
    make data N=$dim_real >> $OUTPUT
    
    for j in "${iters[@]}"
    do  
        make ultra_fast_test_distributed MPI_PROCS=$MPI_PROCS  >> $OUTPUT
    done
done

make clean >> $OUTPUT



