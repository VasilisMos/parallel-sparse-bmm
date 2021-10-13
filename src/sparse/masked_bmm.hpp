#ifndef MASKED_BMM_HPP
#define MASKED_BMM_HPP

#include "matrix.hpp"
#include "../headers/parameters.hpp"
#include "bmm_blocking.hpp"



void masked_bmm(csc *A, csc *B, csc *C, csc *F);
void masked_bmm_blocking(csc *A, csc *B, csc *C,csc *F);

inline Filter *csc2Filt(csc *A){
    Filter *f = (Filter*)malloc(sizeof(Filter));

    f->col_ptr = A->col_ptr;
    f->row = A->row;
    f->val = (int*)malloc(A->nnz * sizeof(int));
    f->nnz = A->nnz;
    f->rowS = A->rowS;
    f->colS = A->colS;

    for(int i=0;i<f->nnz;i++)
        f->val[i] = 1;

    return f;
}

#endif