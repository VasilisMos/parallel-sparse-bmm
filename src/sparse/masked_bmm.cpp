#include "masked_bmm.hpp"

#define A(p,q) Abl[( (p) * (nb) + (q) )]
#define B(p,q) Bbl[( (p) * (nb) + (q) )]
#define C(p,q) Cbl[( (p) * (nb) + (q) )]
#define F(p,q) Fbl[( (p) * (nb) + (q) )]
#define Aps A(p,s)
#define Bsq B(s,q)
#define Cpq C(p,q)
#define Fpq F(p,q)

void masked_bmm_opt(csc *A, csc *B, csc *C, csc *F){

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
    // int *mask_flag = (int *)calloc( Anrow * sizeof(int), sizeof(int) );
    int *mask_flag = (int *)malloc( Anrow * sizeof(int) );

    int pm_end;

    for (int j = 0 ; j < Bncol ; j++)
    {
        pm_end = Fp[j+1];
        for(int i = Fp[j] ; i < pm_end ; i++)
            mask_flag[ Fi[i] ] = j+1;

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
                if ((mask_flag[i] == j+1) && Flag [i] != mark)
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
    // int *mask_flag = (int *)calloc( Anrow * sizeof(int), sizeof(int) );
    int *mask_flag = (int *)malloc( Anrow * sizeof(int) );

    int pm_end;

    for (int j = 0 ; j < Bncol ; j++)
    {
        pm_end = Fp[j+1];
        for(int i = Fp[j] ; i < pm_end ; i++)
            mask_flag[ Fi[i] ] = j+1;

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
                if ((mask_flag[i] == j+1) && Flag [i] != mark)
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

void masked_bmm_blocking(csc *A, csc *B, csc *C,csc *F){

    /* Blocking Variables */
    int nb = BLOCKING_FACTOR; //Blocking Factor
    int n = A->rowS;
    int b = n/nb;

    csc ** Abl = create_blocks(A,nb); // Matrix A cut in blocks
    csc ** Bbl = create_blocks(B,nb); // Matrix B cut in blocks
    csc ** Cbl = create_blocks(C,nb); // Empty matrix to store C
    csc ** Fbl = create_blocks(F,nb); // Matrix F cut in blocks

    csc **temps = (csc**)malloc( nb * sizeof(csc*) );

    // for(int s=0;s<nb;s++) temps[s] = initCsc(b,b,(A->nnz + B->nnz)/nb/nb*3);

    for (int q = 0; q < nb; q++)          // Calculate C(:,q)
    {
        for (int p = 0; p < nb; p++)      // Being at C(p,q) = 0;
        {
            for (int s = 0; s < nb; s++)  // C(p,q) = A(p,:)*B(:,q)
            {
                temps[s] = initCsc(b,b,(A->nnz + B->nnz)/nb/nb*3);
                masked_bmm(Aps, Bsq, temps[s], Fpq);
            }
            merge_blocks(Cpq,temps,nb);
        }
    }

    unify_blocks(C, Cbl, nb);
}