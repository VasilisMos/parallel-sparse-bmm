#include "masked_bmm.hpp"


void masked_bmm(csc *A, csc *B, csc *C, csc *F){

    int Bncol = A->colS, Anrow = A->rowS;

    int pb = 0 ;
    int cnz = 0 ;
    int mark = 0 ;
    
    int pbend; int pcstart, pcmax, pa, paend, k, i;

    int *Ap = A->col_ptr;
    int *Bp = B->col_ptr;
    int *Cp = C->col_ptr;
    int *Fp = F->col_ptr;

    int *Ai = A->row;
    int *Bi = B->row;
    int *Ci = C->row;
    int *Fi = F->row;

    int *Flag = (int *)malloc( Anrow * sizeof(int) );
    int *mask_flag = (int *)malloc( Anrow * sizeof(int) );

    int pm_end;

    for (int j = 0 ; j < Bncol ; j++)
    {
        pm_end = Fp[j+1];
        for(int i = Fp[j] ; i < pm_end ; i++)
            mask_flag[ Fi[i] ] = j;

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
                if ((mask_flag[i] == j) && Flag [i] != mark)
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

    /* Sort the row indices for each column
     * (they are unsorted in general) 
     */
    for(int i=0;i<Bncol;i++){
        if(Cp[i+1]-Cp[i])
            mergesort<int>(Ci,Cp[i],Cp[i+1]-1);
    }

    free(Flag);
    free(mask_flag);
}