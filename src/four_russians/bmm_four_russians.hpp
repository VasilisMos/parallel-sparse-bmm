#ifndef BMM_FOUR_RUSSIANS_HPP
#define BMM_FOUR_RUSSIANS_HPP

#include "../sparse/matrix.hpp"
#include "../headers/parameters.hpp"
#include <math.h>


csc *bmm_four_russians(csc *A, csc *B, csc *C);

csc *get_vertical_chunk(csc *A, int start, int end);

int num(int *Ai, int n, int offset);
void or_rows(csr *A, int i1, csr *Rs, int j1, int i_targ, int offset);

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

csr *csc2csr(csc *A);


template <typename T>
inline T *my_malloc(int size){

  T* res = (T*)malloc(size * sizeof(T));
  
  if (res == NULL){
     printf("Error in malloc\n");
     exit(1);
  }
  
  return res;
}

#endif