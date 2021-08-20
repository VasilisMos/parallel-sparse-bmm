#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=batch
#SBATCH --time=00:02:00

OUTPUT=./logs/console_output.txt
num_iters=4
dec_exponent=$((10 ** 4))

cd ..
module purge
module load gcc/7.3.0
module load octave

dims=( 2 3 4 5 6 8 )
iters=( $(seq 1 ${num_iters} ) )

make clear_times >> $OUTPUT
make sequential >> $OUTPUT
 
for n in "${dims[@]}"
do
    dim_real=$(($n * $dec_exponent))
    make data N=$dim_real >> $OUTPUT
    for j in "${iters[@]}"
    do  
        ./seq.out >> $OUTPUT
    done
done

make clean >> $OUTPUT