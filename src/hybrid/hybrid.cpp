#include <mpi.h>
#include "bmm_hybrid.hpp"

void *run_bmm(int proc_num, int total_procs, int threads){
    if(proc_num) run_bmm_slave(proc_num, total_procs, threads);
    else run_bmm_master(total_procs, threads);
}

int main(int argc, char* argv[]){

    if(!(argc==2)) { cout << "Not correct args" << endl; exit(1); }

    int threads = atoi(argv[1]);

    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    run_bmm(world_rank,world_size, threads);

    // Finalize the MPI environment.
    // printf("P%d Exiting Gracefully\n",world_rank);
    MPI_Finalize();
}