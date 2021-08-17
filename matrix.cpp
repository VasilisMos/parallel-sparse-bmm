#include "matrix.hpp"

#define A(x,y) A[ (x) * n2 + (y) ]
#define B(x,y) B[ (x) * n3 + (y) ]
#define C(x,y) C[ (x) * n3 + (y) ]

/*
 * Boolean Matrix Multiplication w/o blocking
 * Tested and Working on MATLAB Benchmark with
 * 'sequential1.cpp'
 */
void bmm(csc *A,csc *B, csc *C){
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

    /* Sort the row indices for each column
     * (they are unsorted in general) 
     */
    for(int i=0;i<Bncol;i++){
        if(Cp[i+1]-Cp[i])
            mergesort<int>(Ci,Cp[i],Cp[i+1]-1);
    }

    free(Flag);
}

void bmm(csr *A,csc *B,csr *C){
    int n1 = A->rowS, n2 = A->colS, n3 = B->colS;
    int nnz = 0, difA,difB;
    int *ptrA,*ptrB;
    int c1,c2;

    C->r_p[0] = 0;
    
    for(int i=0;i<n1;i++){
        difA = A->r_p[i+1]-A->r_p[i];
        C->r_p[i+1] = C->r_p[i];
        if(!difA) continue; 
        

        for(int j=0;j<n3;j++){
            
            difB = B->col_ptr[j+1]-B->col_ptr[j];     if(!difB) continue;

            ptrA = A->col+A->r_p[i];
            ptrB = B->row+B->col_ptr[j];

            for(c1 = 0,c2=0; c1<difA && c2<difB; ){
                if(ptrA[c1] == ptrB[c2]) {
                    C->r_p[i+1]++;
                    C->col[nnz++] = j;
                    break;
                }
                else if(ptrA[c1] > ptrB[c2]) c2++;
                else c1++;
            }
        }
    }


    C->nnz = nnz;
    //printf("C dims = (%d,%d), nnz = %d\n", C->rowS, C->colS,nnz);
}

/*void bmm_slow(csr *A,csc *B,csr *C){
    int n1 = A->rowS, n2 = A->colS, n3 = B->colS;
    int nnz = 0;

    C->r_p[0] = 0;
    
    for(int i=0;i<n1;i++){
        int difA = A->r_p[i+1]-A->r_p[i];
        if(!difA) continue;
        C->r_p[i+1] = C->r_p[i];

        for(int j=0;j<n3;j++){
            int difB = B->col_ptr[j+1]-B->col_ptr[j]; 
            if(!difB) continue;

            int bit = sparse_and(A->col+A->r_p[i], B->row+B->col_ptr[j], difA ,difB);
            if(bit){
                C->r_p[i+1]++;
                C->col[nnz++] = j;
            }
        }
    }
    C->nnz = nnz;
}*/

int sparse_and_naive(int *A, int *B, int n1, int n2){
    int flag = 0;

    for(int c1 = 0,c2=0; c1<n1 && c2<n2; ){
        if(A[c1] == B[c2]) flag=1;

        if(A[c1] > B[c2]) c2++;
        else c1++;
    }

    return flag;
}

int sparse_and(int *A, int *B, int n1, int n2){

    for(int c1 = 0,c2=0; c1<n1 && c2<n2; ){
        if(A[c1] == B[c2]) return 1;

        if(A[c1] > B[c2]) c2++;
        else c1++;
    }

    return 0;
}

void bmm_wrapper(void (*bmm_implementation)(csc*, csc*, csc*), char *f1, char *f2, char *fout ){
    csc *A = (csc*)parse_data(f1, CSC);
    csc *B = (csc*)parse_data(f2, CSC);
    csc *C = initCsc(A->rowS,B->colS, 2*(A->nnz + B->nnz));
    print_version(A,B,C);

    tic();
    bmm_implementation(A,B,C);
    toc();

    write_mtx_csc(C, fout);
    destroyCsc(A); destroyCsc(B); destroyCsc(C);
}

/*
 * Boolean matrix multiplication on dense matrices A,B that produces a dense
 * C matrix as a result. We assume that A is n1xn2 and B is n2xn3 (so C is n1xn3)
 */
void bmm(int *A, int *B, int *C, int n1, int n2, int n3){
    
    for(int i=0;i<n1;i++){
        for(int j=0;j<n3;j++){
            C(i,j) = 0;
            for(int k=0;k<n2;k++)
                C(i,j) = C(i,j) || (A(i,k) * B(k,j)); 

        }
    }
}

/*
 * Boolean matrix multiplication on dense matrices A,B that produces a dense
 * C matrix as a result. We assume that A is n1xn2 and B is n2xn3 (so C is n1xn3)
 * 
 * It is output sensitive, so if a c(i,j) becomes 1 at the middle of the computation,
 * the remaining v( A(i,:) ^ B(:,j) ) are not calculated, as they won't affect the result
 */
void bmm2(int *A, int *B, int *C, int n1, int n2, int n3){
    for(int i=0;i<n1;i++){
        for(int j=0;j<n3;j++){
            C(i,j) = 0;
            for(int k=0;k<n2;k++){
                C(i,j) = C(i,j) || (A(i,k) * B(k,j)); 
                if(C(i,j)) break;
            }
        }
    }
}


/*
 * Counts the nnz(A*B), where A,B sparse matrices
 * source: https://github.com/PetterS/SuiteSparse/blob/master/MATLAB_Tools/SSMULT/ssmult_saxpy.c
 */
int count_nnz(csc *A, csc *B, csc *C){

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

    //TODO
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
                    cnz++ ;
                }
            }
        }
    }

    Cp[Bncol] = cnz;

    free(Flag);
    return cnz;
}
