#include "bmm_distributed2.hpp"
#define RING_ENABLE 1
#define curr_pos ( ((proc_num-i) % total_procs) < 0 ? ((proc_num-i) % total_procs) + total_procs : ((proc_num-i) % total_procs) )
MPI_Request reqs[ 3 * 4 ];

csc *block_bmm(csc **A,csc **B,csc ** temps, int nb,int pid){

    csc *C = initCsc( B[0]->rowS, B[0]->colS, 5 * (A[0]->nnz + B[0]->nnz) );

    for (int s = 0; s < nb; s++) // C(p,q) = A(p,:)*B(:,q)
        bmm(A[s], B[s], temps[s]);

    merge_blocks(C, temps, nb);

    return C;
}

void run_bmm_master(int total_procs){

    int sqrt_total_procs = (int) sqrt(total_procs);
    int nb = sqrt_total_procs;

    csc **Abl, **Bbl, *C, **Cbl = (csc **)malloc( nb * sizeof(csc *) );
    struct timespec initiation, finishing, t0, t1;
    csc **temps = (csc **)malloc( nb * sizeof(csc *) );

    /* Read datasets */

    csc **Abl_total, **Bbl_total, **Cbl_total;

    csc *A = (csc *)parse_data(fname1, CSC);
    csc *B = (csc *)parse_data(fname2, CSC);
    C = initCsc(A->rowS, B->colS, 3 * (A->nnz + B->nnz));  print_version(A, B, C); printf("Num of Processes:%d\n", total_procs); 
    initiation = TIC

    t0 = TIC
    Abl_total = create_blocks(A, nb); // Matrix A cut in blocks
    Bbl_total = create_blocks(B, nb); // Matrix B cut in blocks
    t1 = TOC 
    cout << "Block Creation:"; time_elapsed(t1,t0);


    Cbl_total = (csc **)malloc( nb * nb * sizeof(csc *));
    for (int s = 0; s < nb; s++)   temps[s] = initCsc(Abl_total[0]->rowS, Bbl_total[s]->colS, 3 * (Abl_total[s]->nnz + Bbl_total[s]->nnz)); 

    Abl = (csc **)malloc( nb * sizeof(csc *));
    Bbl = (csc **)malloc( nb * sizeof(csc *));

    // Get the 1st block row from A and 1st block column from B for proccess 0
    for (int s = 0; s < nb; s++)
    {
        Bbl[s] = Bbl_total[s * nb];
        Abl[s] = Abl_total[s];
    }

    t0 = TIC

    //Distribute Datasets to processes
    for (int p = 1; p < total_procs; p++)
    {
        int i = p / nb;
        int j = p % nb;

        // printf("Proc:%d sqrt()=%d, corresponds to (i,j)=(%d,%d)\n",p,sqrt_total_procs,i,j);
        
        for (int row = 0; row < nb; row++)
            send_block(Bbl_total[row * nb + j], sqrt_total_procs, MPI_DEST);

        for (int col = 0; col < nb; col++)
            send_block(Abl_total[i * nb + col], sqrt_total_procs, MPI_DEST);
    }

    t1 = TOC

    cout << "Sending Operations..."; time_elapsed(t1,t0);

    // cout << "Blocks Sent from root" << endl;

    Cbl_total[0] = block_bmm(Abl, Bbl, temps, nb, MASTER);

    MPI_Barrier(MPI_COMM_WORLD);

    for (int p = 1; p < total_procs; p++)
    {
        int i = p / sqrt_total_procs;
        int j = p % sqrt_total_procs;

        csc *temp = rcv_block(p);

        Cbl_total[ i * sqrt_total_procs + j ] = temp;
    }

    //Merge results
    t0 = TIC
    unify_blocks(C, Cbl_total, nb); finishing = TOC 

    cout << "Unify Block: "; time_elapsed(finishing,t0);
    cout << "Total time: "; time_elapsed(finishing, initiation);

    //Finalize
    write_times(A->rowS,initiation,finishing,DISTRIBUTED,total_procs);
    write_mtx_csc(C, fname3); 
}

void run_bmm_slave(int proc_num, int total_procs){

    int sqrt_total_procs = (int) sqrt(total_procs);
    int nb = sqrt_total_procs;

    //Get Dataset from MASTER Proccess
    csc **Bbl = rcv_blocks2(nb,MASTER);
    csc **Abl = rcv_blocks2(nb,MASTER);
    csc * Cbl = (csc*)malloc(sizeof(csc));
    csc **temps = (csc **)malloc(nb * sizeof(csc *));

    //cout << proc_num << ",Received Abl, Bbl " << endl;

    for (int s = 0; s < nb; s++)   temps[s] = initCsc(Abl[s]->rowS, Bbl[s]->colS, 3 * (Abl[s]->nnz + Bbl[s]->nnz)); 

    Cbl = block_bmm(Abl, Bbl, temps, nb, MASTER);
    MPI_Barrier(MPI_COMM_WORLD);

    // Send results to MASTER proccess
    send_block(Cbl, nb, MASTER);
}

// End of File
