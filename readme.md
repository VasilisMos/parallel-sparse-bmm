# Boolean Matrix Multiplication using Parallel Programming in C++

This project implements a parallelization of Boolean Matrix Multiplication **A*B=C** in `C++`, where **A,B,C** are Sparse Matrices in CSC (*Compressed sparse column*) format. The parallelization is achieved with in two stages (Multithreading:OpenMP, Distributed Programming:MPI)

Software prerequisites to run this Project:

1. `Octave`
2. `g++` compiler
3. `OpenMPI` Library

## Generate Matrices A,B,C

To generate Sparse Matrices **A,B,C**, simply open an `Octave` command prompt and move to the ./matlab directory of this Project and type:

    >> generate_datasets;

Alternatively, simply execute from project ./src directory:

    $make data
    $make data N=5e6

You can also change the dimensions and the sparsity level  of produced matrices at the top of file ./src/bash_src/set_parameters.sh.

## Execute Boolean Matrix Multiplication

On a linux terminal (again on the project ./src folder) run (for testing and purge, respectively):

    make all
    make filter N=5e6 D=2
    
    ./seq.out
    ./multithreaded.out 4
    mpirun -np X ./distributed.out %or make all && make ultra_fast_test_distributed MPI_PROCS=X
    mpirun -np 4 ./distributed2.out %or make all && make ultra_fast_test_distributed2 MPI_PROCS=4
    mpirun -np 4 ./hybrid.out X
    mpirun -np X ./hybrid2.out Y

    make clean
    make clear_times

## Results Validation

The validator at ./src/matlab/validate_bmm_results.m reads **A,B,C_result** from their respective files and checks whether **A*B = C_result** ('*' sign for boolean matrix multiplication operator).

On `Octave`, move to project ./src directory and execute on command line:

    >> validate_bmm_results_octave;

## References

[1] MATLAB ".mtx" files I/O: <https://math.nist.gov/MatrixMarket/mmio/matlab/mmiomatlab.html>

[2] Louridas - Four Russians Algorithm https://louridas.github.io/rwa/assignments/four-russians/

[3] PetterS - SuiteSparse Matrix Package https://github.com/PetterS/SuiteSparse
