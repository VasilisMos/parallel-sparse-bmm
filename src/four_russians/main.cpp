#include "bmm_four_russians.hpp"

void test_bmm_four_russians(){

    struct timespec t1,t2;

    csc *A = (csc*)parse_data(fname1, CSC);
    csc *B = (csc*)parse_data(fname2, CSC);
    csc *C;// = initCsc(A->rowS,B->colS, 2*(A->nnz + B->nnz));
    
    print_version(A,B,C);
    t1 = tic();
    C = bmm_four_russians(A,B,C); t2 = toc(); time_elapsed(t2,t1);
    std::cout << std::endl;

    write_times(A->rowS, t1,t2, SEQUENTIAL,1);
    write_mtx_csc(C, fname3);
    destroyCsc(A); destroyCsc(B); // destroyCsc(C);
}

void test_csc2csr(){
    #define fname4 "../datasets/test/C_test.mtx"

    csc *A = (csc*)parse_data(fname4, CSC);

    // csc* A = initCsc(5,2,5);
    // free(A->col_ptr);
    // free(A->row);

    // int col_ptr[] = {0, 4, 5, 5, 5, 5};
    // int row[] = {0, 2, 3, 4, 0};
    // A->col_ptr = col_ptr;
    // A->row = row;

    printCsc(A);

    csr *B = csc2csr(A);
    printf("dims=(%d,%d)\n",B->rowS,B->colS);
    printCsr(B);

    write_mtx_csr(B, fname3);
    destroyCsr(B);
}

void test_vertical_chunking(){
#define fname4 "../datasets/test/C_test.mtx"

    csc *A = (csc*)parse_data(fname4, CSC);
    printCsc(A);

    csc *Ai = get_vertical_chunk(A,0); printCsc(Ai); destroyCsc(Ai);
    Ai = get_vertical_chunk(A,3);      printCsc(Ai); destroyCsc(Ai);
    Ai = get_vertical_chunk(A,6);     printCsc(Ai); destroyCsc(Ai);
    Ai = get_vertical_chunk(A,9);     printCsc(Ai); destroyCsc(Ai);
}

void test_horizontal_chunking(){
#define fname4 "../datasets/test/C_test.mtx"

    csc *A = (csc*)parse_data(fname4, CSC);
    csr *B = csc2csr(A);
    printCsr(B);

    csr *Bi = get_horizontal_chunk(B,0); printCsr(Bi); destroyCsr(Bi);
    Bi = get_horizontal_chunk(B,4);      printCsr(Bi); destroyCsr(Bi);
    Bi = get_horizontal_chunk(B,8);     printCsr(Bi); destroyCsr(Bi);

    destroyCsc(A); destroyCsr(B);
}

int main(){
    test_bmm_four_russians();

//    test_horizontal_chunking();
//    test_vertical_chunking();
//    test_csc2csr();

    cout << "Main is exiting successfully.." << endl;
}