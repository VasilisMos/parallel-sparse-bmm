#include "masked_bmm_multithread.hpp"

#define fnameF "../datasets/test/filter.mtx"

csc *masked_block_bmm(csc **A,csc **B,csc ** temps, csc *F, int nb,int pid){

    for (int s = 0; s < nb; s++) // C(p,q) = A(p,:)*B(:,q)
        masked_bmm(A[s], B[s], temps[s], F);

    csc *C = initCsc( B[0]->rowS, B[0]->colS, 5 * (A[0]->nnz + B[0]->nnz) );
    merge_blocks(C, temps, nb);

    return C;
}

csc *masked_bmm_multithreads(csc **Abl_total, csc **Bbl_total, csc **F_total, int total_procs){

    int b = Abl_total[0]->rowS, nb = total_procs, n = nb* b;
    int C_nnz = 10 * nb * (nnz(Abl_total[0]) + nnz(Bbl_total[0]));

    csc *C = initCsc(n,n, C_nnz);
    csc **Cbl_total = (csc**)malloc( total_procs * total_procs * sizeof(csc*) );
    omp_set_num_threads(total_procs);

    struct timespec t1,t2;

    cout << "Parallel Region:"; t1 = tic();

#pragma omp parallel shared(Abl_total, Bbl_total,Cbl_total,F_total) num_threads(total_procs)
{
    int proc_num = omp_get_thread_num();

    csc **Abl = (csc **)malloc( total_procs * sizeof(csc*) );
    csc **Bbl = (csc **)malloc( total_procs * sizeof(csc*) );
    csc **Fbl = (csc **)malloc( total_procs * sizeof(csc*) );
    csc **temps = (csc **)malloc( total_procs * sizeof(csc*) );

    for(int i=0;i<total_procs;i++){
        Bbl[i] = Bbl_total[ i * total_procs + proc_num ];
        Fbl[i] = F_total[ i * total_procs + proc_num ];
        // temps[i] = initCsc(Bbl[i]->rowS, Bbl[i]->colS, 3 * (Bbl[i]->nnz + Bbl[i]->nnz) );
    }

    for (int i = 0; i < total_procs; i++){
        for (int j = 0; j < total_procs; j++)
            Abl[j] = Abl_total[ i * total_procs + j];

        temps[i] = initCsc(Bbl[i]->rowS, Bbl[i]->colS, 3 * (Bbl[i]->nnz + Bbl[i]->nnz) );
        Cbl_total[ i * total_procs + proc_num ] = masked_block_bmm(Abl, Bbl, temps, Fbl[i], total_procs, proc_num);
    }
}
    t2 = toc(); time_elapsed(t2,t1);
    
    cout << "Unify blocks:";  t1 = tic();
    unify_blocks(C, Cbl_total, total_procs);
    t2 = toc(); time_elapsed(t2,t1);

    return C;
}

void run_bmm_masked_mulithreaded(int total_procs){

    csc **Abl_total, **Bbl_total, **F_total, *C;
    struct timespec initiation, finishing,t0,t1;

    csc *A = (csc *)parse_data(fname1, CSC);
    csc *B = (csc *)parse_data(fname2, CSC);
    csc *F = (csc *)parse_data(fnameF, CSC);
    C = initCsc(A->rowS, B->colS, 6 * (A->nnz + B->nnz));  print_version(A, B, C); printf("Num of Processes:%d\n", total_procs); 
    initiation = TIC 

    t0 = TIC
    Abl_total = create_blocks(A, total_procs); // Matrix A cut in blocks
    Bbl_total = create_blocks(B, total_procs); // Matrix B cut in blocks
    F_total = create_blocks(F, total_procs);   // Filter M cut in blocks
    t1 = TOC 

    C = masked_bmm_multithreads(Abl_total, Bbl_total, F_total, total_procs);

    finishing = TOC 
    cout << "Block Creation:"; time_elapsed(t1,t0);
    cout << "Total time: "; time_elapsed(finishing, initiation);
    cout << "C nnz=" << C->nnz << endl;

    //Finalize
    write_times(A->rowS,initiation,finishing,MULTITHREADED,total_procs);
    write_mtx_csc(C, fname3); 

    destroyCsc(A); destroyCsc(B); destroyCsc(C);
}