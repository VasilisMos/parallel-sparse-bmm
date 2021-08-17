#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <typeinfo>
#include "sparse.hpp"

#include "headers/sorting.hpp"
#include "headers/utils.hpp"
#include "headers/my_time.hpp"

#define max(x,y) ( ((x)>(y)) ? (x) : (y) )
#define min(x,y) ( ((x)<(y)) ? (x) : (y) )

using namespace std;

void bmm_wrapper(void (*bmm_implementation)(csc*, csc*, csc*), char *f1, char *f2, char *fout, int TYPE, int procs );

void bmm(csr *A,csc *B,csr *C);
void bmm(csc *A, csc*B, csc *C);

void bmm(int *A, int *B, int *C, int n1, int n2, int n3);
void bmm2(int *A, int *B, int *C, int n1, int n2, int n3);

int count_nnz(csc *A,csc *B, csc *C);

int sparse_and(int *A, int *B, int n1, int n2);
int sparse_and_naive(int *A, int *B, int n1, int n2);

inline int *rand_bool(int row_S,int col_S){
    int *A = (int *)malloc( row_S * col_S * sizeof(int) );

    srand(time(NULL));

    for(int i=0;i<row_S*col_S;i++)
        A[i] = rand()%2;

    return A;
}

template <typename T>
inline T* zeros(int row_S,int col_S){
    T *A = (T *)malloc( row_S * col_S * sizeof(T) );

    for(int i=0;i<row_S*col_S;i++)
        A[i] = (T) 0;

    return A;
}

template <typename T>
inline bool assertMatrices(T* A, T*B, int n1,int n2,bool print){
    bool flag = 1;

    for(int i=0;i<n1*n2;i++)
        if( A[i] != B[i] ) { flag = 0; break; }
    
    const char *msg = flag ? "Matrices are equal" : "Matrices are NOT equal";

    if(print) cout << msg << endl;

    return flag;
}

#endif










