#ifndef TEST_HPP
#define TEST_HPP

#include "../matrix.hpp"
#include "my_time.hpp"
#include "parameters.hpp"

/*inline void test_bmm_sparse(){

    csr *A = (csr*)parse_data(fname1, CSR);
    csc *B = (csc*)parse_data(fname2, CSC);

    csr *C = initCsr(A->rowS, B->colS, 2*A->nnz + 2*B->nnz);

    std::cout << "Boolean Matrix Multiplication" << std::endl; tic();
    bmm(A,B,C); toc();
    write_mtx_csr(C, fname3);

    destroyCsr(A); destroyCsc(B); destroyCsr(C);
}*/

inline void test_bmm_sequential_dense(){
    int n = 1000;
    int n1 = n, n2 = n, n3 = n;

    int *A = rand_bool(n1,n2);
    int *B = rand_bool(n2,n3);
    int *C = zeros<int>(n1,n3);
    int *C2 = zeros<int>(n1,n3);

tic(); bmm(A,B,C,n1,n2,n3); toc();
tic(); bmm2(A,B,C2,n1,n2,n3); toc();

    assertMatrices<int>(C,C2,n1,n3,true);

    free(A); free(B); free(C); free(C2);
}

inline void test_COO_CSR_conversions(char *fn1){
    coo *A_coo = readCoo(fn1);
    csr *A_csr = coo2csr(A_coo);
    csc *A_csc = coo2csc(A_coo);

    cout << "Print COO:" << endl; printCoo(A_coo);
    cout << "Print CSR:" << endl; printCsr(A_csr);
    cout << "Print CSC:" << endl; printCsc(A_csc);

    destroyCoo(A_coo);  destroyCsr(A_csr); destroyCsc(A_csc);
}

inline void test_bmm_sparse_fast(){

    csc *A = (csc*)parse_data(fname1, CSC);
    csc *B = (csc*)parse_data(fname2, CSC);
    csc *C = initCsc(A->rowS,B->colS, 2*(A->nnz + B->nnz));
    
    print_version(A,B,C);
    tic(); bmm(A,B,C); toc(); std::cout << std::endl;

    write_mtx_csc(C, fname3);
    destroyCsc(A); destroyCsc(B); destroyCsc(C);
}

#endif