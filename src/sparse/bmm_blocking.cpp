#include "bmm_blocking.hpp"
#include <omp.h>

#define A(p,q) Abl[( (p) * (nb) + (q) )]
#define B(p,q) Bbl[( (p) * (nb) + (q) )]
#define C(p,q) Cbl[( (p) * (nb) + (q) )]
#define Aps A(p,s)
#define Bsq B(s,q)
#define Cpq C(p,q)

template <typename T>
T *my_malloc(int size){
    T* res = (T*)malloc(size * sizeof(T));

    if (res == NULL){
        printf("Error in malloc\n");
        exit(1);
    }

    return res;
}

int *v;

void bmm_blocking(csc *A, csc *B, csc *C){

    /* Blocking Variables */
    int nb = BLOCKING_FACTOR; //Blocking Factor
    int n = A->rowS;
    int b = n/nb;

    csc ** Abl = create_blocks(A,nb); // Matrix A cut in blocks
    csc ** Bbl = create_blocks(B,nb); // Matrix B cut in blocks
    csc ** Cbl = create_blocks(C,nb); // Empty matrix to store C

    csc ** temps = (csc**)malloc( nb * sizeof(csc*) );
    v = my_malloc<int>(A->nnz);

    for (int q = 0; q < nb; q++)          // Calculate C(:,q)
    {

        for(int i=0;i<nb;i++) temps[i] = initCsc(b,b,2*(A->nnz + B->nnz));

        for (int p = 0; p < nb; p++)      // Being at C(p,q) = 0;
        {
            for (int s = 0; s < nb; s++) // C(p,q) = A(p,:)*B(:,q)
            {
//                temps[s] = initCsc(b,b, (A->nnz + B->nnz)/nb/nb*3);
                bmm(Aps, Bsq, temps[s]);
            }
            merge_blocks(Cpq,temps,nb);
        }
    }

    unify_blocks(C, Cbl, nb);
    free(v);
}

void merge_blocks(csc *dest, csc **temps, int nb){

    csc *Cb = dest;  Cb->nnz = 0; Cb->col_ptr[0] = 0;
    int *v = my_malloc<int>(nb * temps[0]->nnz + 2);
    int colSize = dest->colS;

    // For each Cpq(:,j)
    for(int j = 0, nnz_col=0; j<colSize; j++,nnz_col=0 ){
        // For each subblock of total nb subblocks
        for(int s=0; s<nb;s++){
            csc *temp = temps[s];
            int start = temps[s]->col_ptr[j];
            int end = temps[s]->col_ptr[j+1];

            for(int i=0;i<end-start;i++)
                v[nnz_col++] = temps[s]->row[start+i];

        }

        mergesort<int>(v,0,nnz_col-1);

        // Remove Duplicates
        nnz_col = unique<int>(Cb->row + Cb->col_ptr[j], v, nnz_col);

        // printV<int>(Cb->row + Cb->col_ptr[j],nnz_col);

        // Correct col_ptr
        Cb->col_ptr[j+1] = Cb->col_ptr[j] + nnz_col;
        Cb->nnz +=nnz_col;
    }

    // Clear temp csc sub blocks
    // for(int i=0;i<nb;i++) destroyCsc(temps[i]);
}

csc ** create_blocks(csc *A,int nb){

    int n = A->colS;
    int nnz = A->nnz;
    int b = n/nb;


    //printf("Block Creation, nb = %d, b = %d\n",nb,b);
    csc **Abl = (csc**)malloc( nb * nb * sizeof(csc*));

    //Initialize all matrix sub-blocks
    for (int q = 0; q < nb; q++)
        for (int p = 0; p < nb; p++){
            A(p,q) = initCsc(b,b,3*nnz/nb);
            A(p,q)->nnz = 0;
        }


    for (int j=0; j<n;j++) // Scan A(:,j)
    {
        int pend = A->col_ptr[j+1];
        int q = j/b; int p;
        //printf("Column %d\n",j);
        for(int pA = A->col_ptr[j]; pA < pend; pA++ )
        {
            int k = A->row[pA];

            //printf("A(%d,%d)->col_ptr[%d]++\n",k/b,q,j - (q*b));
            A(k/b,q)->col_ptr[j - (q*b)+1]++;
            A(k/b,q)->row[(A(k/b,q)->nnz)++] = k - (k/b)*b;
        }
    }

    /*Normalize col_ptr values and parameters for each matrix subblock*/
    for (int p = 0; p < nb; p++)
    {
        for (int q = 0; q < nb; q++)
        {
            for (int j = 1; j < b; j++)
            {
                A(p,q)->col_ptr[j] += A(p,q)->col_ptr[j-1];
            }
            A(p,q)->col_ptr[b] = A(p,q)->nnz;
            A(p,q)->rowS = b;
            A(p,q)->colS = b;
        }
    }

    // Debug Printfs
//    for (int p = 0; p < nb; p++)
//        for (int q = 0; q < nb; q++)
//            n < 10 && (PRINT_2) ? print_block(A(p, q), p, q) : (void)0;

    return Abl;
}

csc *unify_blocks(csc *C, csc **Cbl, int nb){

    // Find total nnz to make proper memory allocation
    int nnz = 0;
    for (int i = 0; i < nb*nb; i++)
        nnz += Cbl[i]->nnz;

    int rowS = Cbl[0]->rowS * nb;
    int colS = Cbl[0]->colS * nb;

    int b = rowS / nb;

    C->col_ptr[0] = 0;
    C->nnz = nnz;

    int j = 0;
    int count = 0;

    for (int q = 0; q < nb; q++)
    {
        for (int j_loc = 0; j_loc < rowS / nb; j_loc++,j++)
        {
            for (int p = 0; p < nb; p++)
            {
                Cpq->col_ptr[0] = 0;
                int start = Cpq->col_ptr[j_loc];
                int end = Cpq->col_ptr[j_loc+1];

                for (int i = 0; i < end - start; i++){
                    C->row[count++] = Cpq->row[start + i] + b*p;
                }
            }
            C->col_ptr[j+1] = count;
        }
    }

    // printf("dim(C)=(%d,%d), nnz = %d\n",C->rowS,C->colS,C->nnz);

    return C;
}
