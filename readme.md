# Boolean Matrix Multiplication using Parallel Programming in C++

This project implements a parallelization of Boolean Matrix Multiplication **A*B=C** in `C++`, where **A,B,C** are Sparse Matrices in CSC (*Compressed sparse column*) format. The parallelization is achieved with in two stages. The first, second,

Software prerequisites to run this Project:

1. `Octave` (also works with `MATLAB`)
2. `g++` compiler
3. `Python3`

## How Computation is done

TODO.

## Generate Matrices A,B,C

To generate Sparse Matrices **A,B,C**, simply open a `MATLAB\Octave` command prompt and move to the ./matlab directory of this Project and type:

    >> generate_datasets;

Alternatively, if `Octave` is installed, simply execute from project directory:

    $make data

You can change the dimensions and the sparsity level  of produced matrices at the top of file ./matlab/generateDatasets.m.

## Execute Boolean Matrix Multiplication

To run the project on Windows with `gcc` run the following command on a cmd at the project folder:

    .\make.ps1 test

Clean project executables, results with:

    .\make.ps1 clean

Similarly, on a linux terminal (again on the project main folder) run (for testing and purge, respectively):

    make test
    make test_omp  
    make test_omp NUM_THREADS=x #to specify number x of spawned threads
    make test_distributed  
    make test_distributed MPI_PROCS=x #to specify proccess number equal to x

    make clean

## Perfomance - Speedup

Here boolean matrix multiplication execution times are presented for each version (MATLAB implementation is the baseline). For the measurements, random generated sparse matrices where used with sparsity 1%. 

| NxN | 1000x1000 | 5000x5000 | 25000x25000 |
| --- | ----------- | ------------- | ------------- |
| Sequential - Double Precision (MATLAB) |   |  |  |
| Sequential - Boolean (C++) |  |  |   |
| OpenMP (8 cores) |   |  |  |
| MPI (8 nodes) |  |  |  |
| MPI (8 nodes) + OpenMP (4 cores) per node |  |  |  |

## Results Validation

To validate the produced results from the BMM type on a `MATLAB` command prompt at the location ./matlab:

    >> valiadate_bmm_results;

The validator at ./matlab/validate_bmm_results.m reads **A,B,C_result** from their respective files and checks whether **A*B = C_result** ('*' sign for boolean matrix multiplication operator).

Similarly, on `Octave`, move to project directory and execute on command line:

    >> validate_bmm_results_octave;

## References

[1] MATLAB ".mtx" files I/O: <https://math.nist.gov/MatrixMarket/mmio/matlab/mmiomatlab.html>

[2] Louridas - Four Russians Algorithm https://louridas.github.io/rwa/assignments/four-russians/

[3] PetterS - SuiteSparse Matrix Package https://github.com/PetterS/SuiteSparse