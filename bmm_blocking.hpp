#ifndef BMM_BLOCKING_HPP
#define BMM_BLOCKING_HPP

#include "matrix.hpp"
#include "./headers/parameters.hpp"

void bmm_blocking(csc *A, csc *B, csc *C,int PRINT=0);
csc ** create_blocks(csc *A,int nb);
void merge_blocks(csc *dest, csc **temps, int nb);
csc *unify_blocks(csc *C,csc **Cbl, int nb);

inline void print_block(csc *B,int p, int q){
    printf("Block (%d,%d), nnz = %d\n",p,q,B->nnz);
    printCsc(B);
}


#endif