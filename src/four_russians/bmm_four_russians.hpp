#ifndef BMM_FOUR_RUSSIANS_HPP
#define BMM_FOUR_RUSSIANS_HPP

#include "four_russians_utils.hpp"

csc *bmm_four_russians(csc *A, csc *B, csc *C);

csc *get_vertical_chunk(csc *A, int start);
csr *get_horizontal_chunk(csr *B, int start);

void or_rows(csr *A, int i1, csr *Rs, int j1, int i_targ, int offset);
void or_rows2(csr *A, int i1, csr *Rs, int j1, int i_targ, int offset);
csr *or_csr(csr *A, csr *B);
int vector_or(int *r1, int size1, int *r2, int size2, int *dest);

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