#!/bin/bash

num_iters=4
num_procs=( 2 3 4 5 6 8 )
iters=( $(seq 1 ${num_iters} ) )
OUTPUT=../logs/console_output.txt

rm $OUTPUT

make clear_times >> $OUTPUT
make distributed >> $OUTPUT
make data >> $OUTPUT

for i in "${num_procs[@]}"
do
    for j in "${iters[@]}"
    do
        make ultra_fast_test_distributed MPI_PROCS=$i >> $OUTPUT
        # printf "${j}\n"
    done
done

make performance 

make clean && make clear_times >> $OUTPUT

printf "\n"