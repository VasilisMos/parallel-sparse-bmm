#ifndef SPARSE_HPP
#define SPARSE_HPP

#include "../headers/mmio.hpp"
#include <stdlib.h>

#define COO 0
#define CSR 1
#define CSC 2

typedef struct{
    int * row_coo;
    int * col_coo;
    int nnz; //number of nonzero elements
    int rowS; //row_cardinality
    int colS; //column_cardinality
}
coo;

typedef struct {
    int *r_p;
    int *col;
    int nnz; //number of nonzero elements
    int rowS; //row_cardinality
    int colS; //column_cardinality
}
csr;

typedef struct {
    int *row;
    int *col_ptr;
    int nnz; //number of nonzero elements
    int rowS; //row_cardinality
    int colS; //column_cardinality
}
csc;


/*Manipulate Csr and Coo form matrices*/
coo *readCoo(char filename[]);
void *parse_data(char * fname, int type);
csr *Dense2Csr(int *A, int n/*n = nodes*/);

csr *coo2csr(coo *in);
csc *coo2csc(coo *in);

coo *csr2coo(csr *A_csr);
coo *csc2coo(csc *A_csc);
void coo2csc(
  int       * const row,       /*!< CSC row start indices */
  int       * const col,       /*!< CSC column indices */
  int const * const row_coo,   /*!< COO row indices */
  int const * const col_coo,   /*!< COO column indices */
  int const         nnz,       /*!< Number of nonzero elements */
  int const         n,         /*!< Number of rows/columns */
  int const         isOneBased /*!< Whether COO is 0- or 1-based */
);
void write_mtx_csr(csr *A_csr, char *path);
void write_mtx_csc(csc *A_csc, char *path);

void destroyCoo(coo *a);
void destroyCsr(csr *a);
void destroyCsc(csc *a);
coo *initCoo(int rowS,int colS, int nnz);
csr *initCsr(int rowS,int colS, int nnz);
csc *initCsc(int rowS,int colS, int nnz);

void printCoo(coo *a);
void printCsr(csr *a);
void printCsc(csc *a);

inline void print_version(csc *A, csc *B, csc *C){
    printf("\n\nBoolean Matrix Multiplication |M=%d|K=%d|N=%d|nnz(A)=%d|nnz(B)=%d|\n",A->rowS,A->colS,B->colS,A->nnz,B->nnz);
}

inline void csc_info(csc *A){
    printf("CSC Matrix dim = (%d,%d), nnz = %d\n",A->rowS, A->colS, A->nnz);
}

inline int nnz(csc* A) { return A->nnz; }
inline int nnz(csr* A) { return A->nnz; }

#endif