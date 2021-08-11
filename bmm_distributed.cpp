#include "bmm_distributed.hpp"
#define RING_ENABLE 1
#define curr_pos ( ((proc_num-i) % total_procs) < 0 ? ((proc_num-i) % total_procs) + total_procs : ((proc_num-i) % total_procs) )
MPI_Request reqs[ 3 * 4 ];

csc *block_bmm(csc **A,csc **B,csc ** temps, int nb,int pid){

    for (int s = 0; s < nb; s++) // C(p,q) = A(p,:)*B(:,q)
        bmm(A[s], B[s], temps[s]);

    csc *C = initCsc( B[0]->rowS, B[0]->colS, 5 * (A[0]->nnz + B[0]->nnz) );
    merge_blocks(C, temps, nb);

    return C;
}

void run_bmm_master(int total_procs){

    csc **Abl, **Bbl, *C, **Cbl = (csc **)malloc(total_procs * sizeof(csc *));
    struct timespec initiation, finishing;
    struct timespec t0, t1;
    csc **temps = (csc **)malloc(total_procs * sizeof(csc *));

    /* Read datasets */

    csc **Abl_total, **Bbl_total, **Cbl_total;

    csc *A = (csc *)parse_data(fname1, CSC);
    csc *B = (csc *)parse_data(fname2, CSC);
    C = initCsc(A->rowS, B->colS, 3 * (A->nnz + B->nnz));  print_version(A, B, C); printf("Num of Processes:%d\n", total_procs); 
    initiation = TIC

    t0 = TIC
    Abl_total = create_blocks(A, total_procs); // Matrix A cut in blocks
    Bbl_total = create_blocks(B, total_procs); // Matrix B cut in blocks
    t1 = TOC 
    cout << "Block Creation:"; time_elapsed(t1,t0);
    Cbl_total = (csc **)malloc(total_procs * total_procs * sizeof(csc *));
    for (int s = 0; s < total_procs; s++)   temps[s] = initCsc(Abl_total[0]->rowS, Bbl_total[s]->colS, 3 * (Abl_total[s]->nnz + Bbl_total[s]->nnz)); 

    Abl = (csc **)malloc(total_procs * sizeof(csc *));
    Bbl = (csc **)malloc(total_procs * sizeof(csc *));

    // Get the 1st block row from A and 1st block column from B for proccess 0
    for (int s = 0; s < total_procs; s++)
    {
        Bbl[s] = Bbl_total[s * total_procs];
        Abl[s] = Abl_total[s];
    }

    //Distribute Datasets to processes
    for (int p = 1; p < total_procs; p++)
    {
        for (int row = 0; row < total_procs; row++)
            send_block(Bbl_total[row * total_procs + p], total_procs, MPI_DEST);

        for (int col = 0; col < total_procs; col++)
            send_block(Abl_total[p * total_procs + col], total_procs, MPI_DEST);
    }

    // Cbl[0] = block_bmm_unoptimized(Abl, Bbl, total_procs, MASTER);
    Cbl[0] = block_bmm(Abl, Bbl, temps, total_procs, MASTER);


    MPI_Barrier(MPI_COMM_WORLD); t0 = TIC 
    ring(Abl, Bbl,Cbl, temps, MASTER, total_procs); t1 = TOC

    cout << "Ring:"; time_elapsed(t1,t0);

    // Gather Results
    for (int j = 0; j < total_procs; j++)
        Cbl_total[j * total_procs] = Cbl[j];

    for (int p = 1; p < total_procs; p++)
    {
        csc **temp = rcv_blocks(total_procs, p);

        for (int j = 0; j < total_procs; j++){
            Cbl_total[j * total_procs + p] = temp[j]; 
        }          
    }

    //Merge results
    t0 = TIC
    unify_blocks(C, Cbl_total, total_procs); finishing = TOC 

    cout << "Unify Block: "; time_elapsed(finishing,t0);
    cout << "Total time: "; time_elapsed(finishing, initiation);

    //Finalize
    write_times(A->rowS,initiation,finishing,DISTRIBUTED,total_procs);
    write_mtx_csc(C, fname3); 
}

void run_bmm_slave(int proc_num, int total_procs){

    //Get Dataset from MASTER Proccess
    csc **Bbl = rcv_blocks(total_procs,MASTER);
    csc **Abl = rcv_blocks(total_procs,MASTER);
    csc **Cbl = (csc**)malloc( total_procs * sizeof(csc*) ), **temps = (csc **)malloc(total_procs * sizeof(csc *));

    for (int s = 0; s < total_procs; s++)   temps[s] = initCsc(Abl[s]->rowS, Bbl[s]->colS, 3 * (Abl[s]->nnz + Bbl[s]->nnz)); 

    // Cbl[proc_num] = block_bmm_unoptimized(Abl, Bbl, total_procs, MASTER);
    Cbl[proc_num] = block_bmm(Abl, Bbl, temps, total_procs, MASTER);
    MPI_Barrier(MPI_COMM_WORLD);

    ring(Abl, Bbl,Cbl, temps, proc_num, total_procs);

    // Send results to MASTER proccess
    for(int s=0;s<total_procs;s++)
        send_block(Cbl[s], total_procs, MASTER);
}

void ring(csc **Abl, csc**Bbl,csc **Cbl,csc **temps, int proc_num,int total_procs){

    // printf("P%d Starting ring\n",proc_num);

    int buf_rcv[2];
    csc *Abl_send[total_procs], *Abl_received[total_procs];
    csc *Abl_test[total_procs]; for(int i=0;i<total_procs;i++) Abl_test[i] = initCsc(Abl[i]->rowS,Abl[i]->rowS,Abl[i]->nnz);

    for (int j = 0; j < total_procs; j++)  Abl_send[j] = Abl[j];          

    for (int i = 1; i < total_procs; i++)
    {
        for (int j = 0; j < total_procs; j++)
        {
            isend_block(Abl_send[j], total_procs, RIGHT, reqs + (3*j) );

            MPI_Recv(buf_rcv, 2, MPI_INT, LEFT, TAG_PARAM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            Abl_received[j] = initCsc(buf_rcv[0], buf_rcv[0], buf_rcv[1]);
            MPI_Recv(Abl_received[j]->col_ptr, buf_rcv[0] + 1, MPI_INT, LEFT, TAG_COL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(Abl_received[j]->row, buf_rcv[1], MPI_INT, LEFT, TAG_ROW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Barrier(MPI_COMM_WORLD);
        }

        Cbl[curr_pos] = block_bmm(Abl_received, Bbl,temps, total_procs, proc_num);

        for (int j = 0; j < total_procs; j++)
        {
            destroyCsc(Abl_send[j]);
            Abl_send[j] = Abl_received[j];
            Abl_received[j] = NULL;
        }
    }

    // printf("P%d Ending ring\n",proc_num);
    MPI_Barrier(MPI_COMM_WORLD);
}

// End of File




