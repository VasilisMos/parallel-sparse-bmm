#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=batch
#SBATCH --time=00:02:00

OUTPUT=../logs/console_output.txt
num_iters=4
dec_exponent=$((10 ** 5))
THREAD_NUM=8

# cd ..
module purge
module load gcc/8.2.0
module load octave

dims=( 2 3 4 5 6 8 )
iters=( $(seq 1 ${num_iters} ) )

make clear_times >> $OUTPUT
make multithreaded >> $OUTPUT
 
for n in "${dims[@]}"
do
    dim_real=$(($n * $dec_exponent))
    make data N=$dim_real >> $OUTPUT
    for j in "${iters[@]}"
    do  
        ./multithreaded.out ${THREAD_NUM} >> $OUTPUT
    done
done

make clean >> $OUTPUT