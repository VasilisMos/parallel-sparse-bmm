#include "my_mpi_comms.hpp"

#define TAG_PARAM 0
#define TAG_COL 1
#define TAG_ROW 2

void send_block(csc* block,int nb, int receiver){

    int buf[2] = {block->rowS, block->nnz};
    MPI_Send(buf, 2, MPI_INT, receiver, TAG_PARAM, MPI_COMM_WORLD);
    MPI_Send(block->col_ptr, block->rowS + 1, MPI_INT, receiver, TAG_COL, MPI_COMM_WORLD);
    MPI_Send(block->row, block->nnz, MPI_INT, receiver, TAG_ROW, MPI_COMM_WORLD);
}

void isend_block(csc* block, int nb, int receiver, MPI_Request *reqs){
    int buf[2];

    init_buffer(buf, block);

    isend_v(buf, 2, receiver, TAG_PARAM, reqs);
    isend_v(block->col_ptr, buf[0] + 1, receiver, TAG_COL, reqs + 1);
    isend_v(block->row, buf[1], receiver, TAG_ROW, reqs + 2);
}

csc **rcv_blocks(int total_procs, int sender){

    csc ** M = (csc**)malloc( total_procs * sizeof(csc*) );
    int buf[2];

    for (int i = 0; i < total_procs; i++)
    {
        MPI_Recv(buf, 2, MPI_INT, sender, TAG_PARAM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        M[i] = initCsc(buf[0], buf[0], buf[1]);

        MPI_Recv(M[i]->col_ptr, M[i]->rowS + 1, MPI_INT, sender, TAG_COL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(M[i]->row, M[i]->nnz, MPI_INT, sender, TAG_ROW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return M;
}




// init_buffer(buf, Abl_send[j]);
// isend_v(buf, 2, RIGHT, TAG_PARAM, reqs + (3 * j));
// isend_v(Abl_send[j]->col_ptr, buf[0] + 1, RIGHT, TAG_COL, reqs + (3 * j + 1));
// isend_v(Abl_send[j]->row, buf[1], RIGHT, TAG_ROW, reqs + (3 * j + 2));

// void block_bmm(csc **A,csc **B, csc *C,csc ** temps, int nb,int pid){

//     for (int s = 0; s < nb; s++) // C(p,q) = A(p,:)*B(:,q)
//         bmm(A[s], B[s], temps[s]);

//     C = initCsc( B[0]->rowS, B[0]->colS,  A[0]->nnz + B[0]->nnz );
//     merge_blocks(C, temps, nb);
// }