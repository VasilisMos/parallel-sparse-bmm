#!/bin/bash

OUTPUT=../logs/console_output.txt
dec_exponent=$((10 ** 6))

dims=( 2 3 4 5 6 )
threads=( 2 3 4 5 6 8 )
iters=( $(seq 1 ${num_iters} ) )

make clear_times >> $OUTPUT
make multithreaded >> $OUTPUT
 
for n in "${dims[@]}"
do
    dim_real=$(($n * $dec_exponent))
    make data N=$dim_real >> $OUTPUT
    ./seq.out
    make ultra_fast_test_distributed2 MPI_PROCS=4
    
    for j in "${threads[@]}"
    do  
        ./multithreaded.out ${j}
        make ultra_fast_test_distributed MPI_PROCS=${j}
    done
done