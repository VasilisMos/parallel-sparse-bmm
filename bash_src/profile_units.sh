#!/bin/bash

num_iters=4
num_procs=( 2 3 4 5 6 8 )
iters=( $(seq 1 ${num_iters} ) )

make distributed >> temp.txt
make data >> temp.txt
make clear_times >> temp.txt

for i in "${num_procs[@]}"
do
    for j in "${iters[@]}"
    do
        make ultra_fast_test_distributed MPI_PROCS=$i >> temp.txt
        # printf "${j}\n"
    done
done

make performance 

make clean && make clear_times >> temp.txt

printf "\n" && rm temp.txt