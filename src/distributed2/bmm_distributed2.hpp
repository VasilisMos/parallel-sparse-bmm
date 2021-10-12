#ifndef BMM_DISTRIBUTED_HPP
#define BMM_DISTRIBUTED_HPP

#include "../sparse/bmm_blocking.hpp"
#include "../headers/my_mpi_comms.hpp"
#include <mpi.h>
#include <math.h>

#define MASTER 0
#define MPI_DEST p

#define TAG_PARAM 0
#define TAG_COL 1
#define TAG_ROW 2

#define TIC tic();
#define TOC toc();

#define LEFT ((proc_num>0) ? proc_num-1 : total_procs-1)
#define RIGHT ((proc_num != total_procs-1) ? proc_num+1 : 0)

csc *block_bmm(csc **A,csc **B, csc *C,csc ** temps, int nb,int pid);

void run_bmm_master(int total_procs);
void run_bmm_slave(int proc_num, int total_procs);
void ring(csc **Abl, csc**Bbl,csc **Cbl,csc **temps, int proc_num,int total_procs);

inline void neighbors_check(int proc_num, int left, int right){
    printf("P%d, My left is %d, My right is %d \n", proc_num, left, right);
}

#endif
