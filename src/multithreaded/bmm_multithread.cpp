#include "bmm_multithread.hpp"

#define A(p,q) Abl[( (p) * (nb) + (q) )]
#define B(p,q) Bbl[( (p) * (nb) + (q) )]
#define C(p,q) Cbl[( (p) * (nb) + (q) )]
#define Aps A(p,s)
#define Bsq B(s,q)
#define Cpq C(p,q)

csc ** create_blocks_parallel(csc *A,int nb, int total_threads){

    int n = A->colS;
    int nnz = A->nnz;
    int b = n/nb;


    //printf("Block Creation, nb = %d, b = %d\n",nb,b);
    csc **Abl = (csc**)malloc( nb * nb * sizeof(csc*));


    //Initialize all matrix sub-blocks
#pragma omp parallel for shared(Abl,nb,b) num_threads(total_threads)
    for (int q = 0; q < nb; q++)
        for (int p = 0; p < nb; p++){
            A(p,q) = initCsc(b,b,3*nnz/nb);
            A(p,q)->nnz = 0;
        }

#pragma omp parallel for shared(Abl) num_threads(total_threads)
    for (int j=0; j<n;j++) // Scan A(:,j)
    {
        int pend = A->col_ptr[j+1];
        int q = j/b; int p;
        //printf("Column %d\n",j);
        for(int pA = A->col_ptr[j]; pA < pend; pA++ )
        {
            int k = A->row[pA];

            //printf("A(%d,%d)->col_ptr[%d]++\n",k/b,q,j - (q*b));
            A(k/b,q)->col_ptr[j - (q*b)+1]++;
            A(k/b,q)->row[(A(k/b,q)->nnz)++] = k - (k/b)*b;
        }
    }

#pragma omp parallel for shared(Abl) num_threads(total_threads)
    /*Normalize col_ptr values and parameters for each matrix subblock*/
    for (int p = 0; p < nb; p++)
    {
        for (int q = 0; q < nb; q++)
        {
            for (int j = 1; j < b; j++)
            {
                A(p,q)->col_ptr[j] += A(p,q)->col_ptr[j-1];
            }
            A(p,q)->col_ptr[b] = A(p,q)->nnz;
            A(p,q)->rowS = b;
            A(p,q)->colS = b;
        }
    }

    return Abl;
}

csc *block_bmm(csc **A,csc **B,csc ** temps, int nb,int pid){

    for (int s = 0; s < nb; s++) // C(p,q) = A(p,:)*B(:,q)
        bmm(A[s], B[s], temps[s]);

    csc *C = initCsc( B[0]->rowS, B[0]->colS, 5 * (A[0]->nnz + B[0]->nnz) );
    merge_blocks(C, temps, nb);

    return C;
}

void run_bmm_multithread(int total_procs){

    /* Read datasets */

    csc **Abl_total, **Bbl_total, *C;
    struct timespec initiation, finishing,t0,t1;

    csc *A = (csc *)parse_data(fname1, CSC);
    csc *B = (csc *)parse_data(fname2, CSC);
    C = initCsc(A->rowS, B->colS, 6 * (A->nnz + B->nnz));  print_version(A, B, C); printf("Num of Processes:%d\n", total_procs); 
    initiation = TIC

    t0 = TIC
    Abl_total = create_blocks_parallel(A, total_procs,total_procs); // Matrix A cut in blocks
    Bbl_total = create_blocks_parallel(B, total_procs,total_procs); // Matrix B cut in blocks
    t1 = TOC 

    C = bmm_multithreads(Abl_total, Bbl_total, total_procs);

    finishing = TOC 
    cout << "Block Creation:"; time_elapsed(t1,t0);
    cout << "Total time: "; time_elapsed(finishing, initiation);
    cout << "C nnz=" << C->nnz << endl;

    //Finalize
    write_times(A->rowS,initiation,finishing,MULTITHREADED,total_procs);
    // write_mtx_csc(C, fname3); 

    destroyCsc(A); destroyCsc(B); destroyCsc(C);
}

csc *bmm_multithreads(csc **Abl_total, csc **Bbl_total, int total_procs){

    int b = Abl_total[0]->rowS, nb = total_procs, n = nb* b;
    int C_nnz = 10 * nb * (nnz(Abl_total[0]) + nnz(Bbl_total[0]));

    csc *C = initCsc(n,n, C_nnz);
    csc **Cbl_total = (csc**)malloc( total_procs * total_procs * sizeof(csc*) );
    omp_set_num_threads(total_procs);

    struct timespec t1,t2;

    t1 = tic();
#pragma omp parallel shared(Abl_total, Bbl_total,Cbl_total) num_threads(total_procs)
// for(int proc_num=0;proc_num<total_procs;proc_num++)
{
    int proc_num = omp_get_thread_num();

    csc **Abl = (csc **)malloc( total_procs * sizeof(csc*) );
    csc **Bbl = (csc **)malloc( total_procs * sizeof(csc*) );
    csc **temps = (csc **)malloc( total_procs * sizeof(csc*) );

    for(int i=0;i<total_procs;i++){
        Bbl[i] = Bbl_total[ i * total_procs + proc_num ];
        temps[i] = initCsc(Bbl[i]->rowS, Bbl[i]->colS, 3 * (Bbl[i]->nnz + Bbl[i]->nnz) );
    }
       
    for (int i = 0; i < total_procs; i++){
        for (int j = 0; j < total_procs; j++)
            Abl[j] = Abl_total[ i * total_procs + j];

        Cbl_total[ i * total_procs + proc_num ] = block_bmm(Abl, Bbl, temps, total_procs, proc_num);
    }    
}   
    t2 = toc();
    cout << "Parallel Region:"; time_elapsed(t2,t1);
    
    t1 = tic();
    unify_blocks(C, Cbl_total, total_procs);
    t2 = toc();

    cout << "Unify blocks:"; time_elapsed(t2,t1);

    return C;
}
