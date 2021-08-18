#include "./headers/testing.hpp"
#include "bmm_multithread.hpp"
#include "./headers/parameters.hpp"
#include <string.h>
#include <omp.h>

#define GET_INPUT (argc==1 ? 4 : atoi(argv[1]))

#define A(p,q) Abl[( (p) * (nb) + (q) )]
#define B(p,q) Bbl[( (p) * (nb) + (q) )]
#define C(p,q) Cbl[( (p) * (nb) + (q) )]
#define Aps A(p,s)
#define Bsq B(s,q)
#define Cpq C(p,q)

int total_procs;

void bmm_blocking_mult(csc *A, csc *B, csc *C){

    /* Blocking Variables */
    int nb = BLOCKING_FACTOR; //Blocking Factor
    int n = A->rowS;
    int b = n/nb;

    csc ** Abl = create_blocks(A,nb); // Matrix A cut in blocks
    csc ** Bbl = create_blocks(B,nb); // Matrix B cut in blocks
    csc ** Cbl = create_blocks(C,nb); // Empty matrix to store C

    csc **temps = (csc**)malloc( nb * sizeof(csc*) );

    /* Classical bmm useful Variables */
    int nnz = 0;

// #pragma omp parallel for
    for (int q = 0; q < nb; q++)          // Calculate C(:,q)
    {

        for(int i=0;i<nb;i++) temps[i] = initCsc(b,b,A->nnz + B->nnz);

        for (int p = 0; p < nb; p++)      // Being at C(p,q) = 0;
        {
            for (int s = 0; s < nb; s++) // C(p,q) = A(p,:)*B(:,q)
            {
                temps[s] = initCsc(b,b, (A->nnz + B->nnz)/nb/nb*2);
                bmm(Aps, Bsq, temps[s]);

            }
            merge_blocks(Cpq,temps,nb);
        }
    }

    unify_blocks(C, Cbl, nb);
}

int main(int argc, char *argv[]){

    total_procs = GET_INPUT;

    bmm_wrapper(bmm_blocking, fname1, fname2, fname3, SEQUENTIAL, 1);
    run_bmm_multithread(total_procs);
    // bmm_wrapper(bmm_blocking_mult, fname1, fname2, fname3, MULTITHREADED, total_procs);

    printf("Main is exiting successfully\n");
}