#include "../headers/testing.hpp"
#include "../sparse/bmm_blocking.hpp"
#include <string.h>

/* Breaks the original A sparse matrix in square submatrices 
 * Assumption dim(A)/BLOCKING_FACTOR is integer
 * TEST PASSED
*/
void test_blocking_creation(int nb){
    cout << "---------------------------BLOCK CREATION-------------" << endl;
    csc *A = (csc*)parse_data(fname1, CSC); print_version(A,A,A);
    tic(); csc ** A_bl = create_blocks(A,nb);  toc();
    cout << "---------------------------BLOCK CREATION - DONE-------------" << endl;

}

void test_bmm_sparse_blocking(){

    struct timespec t1,t2;

    csc *A = (csc*)parse_data(fname1, CSC);
    csc *B = (csc*)parse_data(fname2, CSC);
    csc *C = initCsc(A->rowS,B->colS, 2*(A->nnz + B->nnz));
    
    print_version(A,B,C);
    t1 = tic(); bmm_blocking(A,B,C); t2 = toc(); time_elapsed(t2,t1);
    std::cout << std::endl;

    write_times(A->rowS, t1,t2, SEQUENTIAL,1);
    write_mtx_csc(C, fname3);
    destroyCsc(A); destroyCsc(B); destroyCsc(C);
}

int main(int argc, char * argv[]){
    
    test_bmm_sparse_fast();
    test_bmm_sparse_blocking();
    // test_blocking_creation(BLOCKING_FACTOR);

    cout << "Main is exiting successfully" << endl;  return 0;
}