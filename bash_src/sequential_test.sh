#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=batch
#SBATCH --time=00:02:00

OUTPUT=./logs/console_output.txt
num_iters=4
dec_exponent=10000

cd ..
module purge
module load gcc/7.3.0
# module load matlab/R2021a
module load octave

dims=( 2 3 4 5 6 8 )
iters=( $(seq 1 ${num_iters} ) )

make clear_times >> $OUTPUT
make sequential >> $OUTPUT
make data
 
for n in "${dims[@]}"
do
    dim_real=$(($n * $dec_exponent))
    for j in "${iters[@]}"
    do  
        # matlab -nodisplay -r "cd('./matlab'); generate_datasets(${dim_real}); exit" >> $OUTPUT
        ./seq.out >> $OUTPUT
    done
done

make clean >> $OUTPUT