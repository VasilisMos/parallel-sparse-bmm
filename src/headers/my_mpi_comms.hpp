#ifndef MY_MPI_COMMS_HPP
#define MY_MPI_COMMS_HPP

#include <mpi.h>
#include "../sparse/bmm_blocking.hpp"

void send_block(csc* block,int nb, int receiver);
void isend_block(csc* block, int nb, int receiver, MPI_Request *reqs);
csc **rcv_blocks(int total_procs, int sender);
inline void isend_v(int *v, int count, int dest, int tag, MPI_Request *req ){ MPI_Isend(v, count, MPI_INT, dest, tag, MPI_COMM_WORLD, req); }

inline void init_buffer(int *buf, csc *A){
    *(buf) = A->rowS;
    *(buf+1) = A->nnz;
}

#endif