#include "../headers/testing.hpp"
#include "../sparse/masked_bmm.hpp"
#include <string.h>

#define fname_Filter "../datasets/test/filter.mtx"

void test_bmm_sparse_masked(){

    struct timespec t1,t2;

    csc *A = (csc*)parse_data(fname1, CSC);
    csc *B = (csc*)parse_data(fname2, CSC);
    csc *C = initCsc(A->rowS,B->colS, 2*(A->nnz + B->nnz));
    csc *F = (csc*)parse_data(fname_Filter, CSC);
    
    print_version(A,B,C);
    t1 = tic(); masked_bmm(A,B,C,F); t2 = toc(); time_elapsed(t2,t1);
    std::cout << std::endl;

    write_times(A->rowS, t1,t2, SEQUENTIAL,1);
    write_mtx_csc(C, fname3);
    destroyCsc(A); destroyCsc(B); destroyCsc(C);
}

int main(){
    test_bmm_sparse_masked();
}