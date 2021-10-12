#include <mpi.h>
#include "bmm_hybrid.hpp"

void *run_bmm(int proc_num, int total_procs){
    if(proc_num) run_bmm_slave(proc_num, total_procs);
    else run_bmm_master(total_procs);
}

int main(int argc, char* argv[]){

    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    run_bmm(world_rank,world_size);

    // Finalize the MPI environment.
    // printf("P%d Exiting Gracefully\n",world_rank);
    MPI_Finalize();
}