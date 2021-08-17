#include "./headers/testing.hpp"
#include "bmm_multithread.hpp"
#include "./headers/parameters.hpp"
#include <string.h>
#include <omp.h>

#define TIC tic();
#define TOC toc();
#define GET_INPUT (argc==1 ? 4 : atoi(argv[1]))

void bmm_omp(csc *A, csc *B, csc *C){
    int Bncol = A->colS, Anrow = A->rowS;

    int pb = 0 ;
    int cnz = 0 ;
    int mark = 0 ;
    
    int pbend; int pcstart, pcmax, pa, paend, k, i;

    int *Ap = A->col_ptr;
    int *Bp = B->col_ptr;
    int *Cp = C->col_ptr;

    int *Ai = A->row;
    int *Bi = B->row;
    int *Ci = C->row;

    int *Flag = (int *)malloc( Anrow * sizeof(int) );

    for (int j = 0 ; j < Bncol ; j++)
    {
        /* Compute nnz (C (:,j)) */
        mark-- ;                        /* Flag [0..n-1] != mark is now true */ 
        pb = Bp [j] ;
        pbend = Bp [j+1] ;
        pcstart = cnz ;
        pcmax = cnz + Anrow ;
        Cp [j] = cnz ;
        /* cnz += nnz (C (:,j)), stopping early if nnz(C(:,j)) == Anrow */
        for ( ; pb < pbend && cnz < pcmax ; pb++)
        {
            k = Bi [pb] ;               /* nonzero entry B(k,j) */
            paend = Ap [k+1] ;
            for (pa = Ap [k] ; pa < paend ; pa++)
            {
                i = Ai [pa] ;           /* nonzero entry A(i,k) */
                if (Flag [i] != mark)
                {
                    /* C(i,j) is a new nonzero */
                    Flag [i] = mark ;   /* mark i as appearing in C(:,j) */
                    Ci[cnz] = i;
                    cnz++ ;
                }
            }
        }
    }

    Cp[Bncol] = cnz;
    C->nnz = cnz;

    for(int i=0;i<Bncol;i++){
        if(Cp[i+1]-Cp[i])
            mergesort<int>(Ci,Cp[i],Cp[i+1]-1);
    }

    free(Flag);
}

int main(int argc, char *argv[]){
    // run_bmm_multithread(GET_INPUT);

    bmm_wrapper(bmm_omp, fname1, fname2, fname3);

    printf("Main is exiting successfully\n");
}