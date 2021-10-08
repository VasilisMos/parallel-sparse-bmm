#ifndef FOUR_RUSSIANS_UTILS_HPP
#define FOUR_RUSSIANS_UTILS_HPP

#include "../sparse/matrix.hpp"
#include "../headers/parameters.hpp"
#include <math.h>

inline void loop_logic(int *bp, int *j, int *k){
    if((*bp)==1){
        *bp=*j+1;
        (*k)++;
    }
    else
        (*bp)--;
}

int num(int *Ai, int n, int offset,int dim);

inline void zero_padding(csc *A, int n){
    int index_max = ((int)ceil(n/log2(n))) * ((int)floor(log2(n)));

    if(index_max>n) //Needs Padding 
        A->colS = A->colS + index_max-n;
}

inline void zero_padding(csr *A, int n){
    int index_max = ((int)ceil(n/log2(n))) * ((int)floor(log2(n)));

    if(index_max>n) //Needs Padding 
        A->rowS = A->rowS + index_max-n;
}





#endif