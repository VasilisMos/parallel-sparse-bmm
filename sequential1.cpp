#include "./headers/testing.hpp"
#include "bmm_blocking.hpp"
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

void profile_blocking_creation(){

    int nb[11] = {1,2,5,10,20,25};
    int n = 6;
    csc *A = (csc*)parse_data(fname1, CSC); print_version(A,A,A);

    struct timespec start,end;

    for(int i=0;i<6;i++){
        printf("iter %d/%d, nb = %d\n",i+1,n,nb[i]);
        
        start = tic();
        csc ** A_bl = create_blocks(A,nb[i]);
        end = toc(); cout << endl;

        // time_elapsed(end, start);

        for(int j=0;j<(nb[i]-1)*(nb[i]-1);j++){
            destroyCsc(A_bl[i]);
        }

        free(A_bl);
    }

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

    test_bmm_sparse_blocking();
    
    // test_bmm_sparse_fast();
    // test_blocking_creation(BLOCKING_FACTOR);
    // profile_blocking_creation();

    cout << "Main is exiting successfully" << endl;  return 0;
}