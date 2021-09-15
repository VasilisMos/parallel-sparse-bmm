#include "sparse.hpp"
#include <iostream>

#define ind(x,y,n) ( (x) * (n) + (y) )
#define isOneBased 0

void write_mtx_csr(csr *A_csr, char *path){
    std::cout << "Write Mtx to disc..";

    FILE *f = fopen(path,"w");
    if (f == NULL) { printf("Could not open mtx file to write\n"); exit(1); }

    coo *A = csr2coo(A_csr);

    fprintf(f,"%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(f,"%d %d %d\n",A->rowS,A->colS,A->nnz);
    for(int i=0;i<A->nnz;i++)
        fprintf(f,"%d %d 1\n",A->row_coo[i]+1,A->col_coo[i]+1);

    destroyCoo(A);

    std::cout << "Done successfully" << std::endl;
}

void write_mtx_csc(csc *A_csc, char *path){
    std::cout << "Write Mtx to disc..";

    FILE *f = fopen(path,"w");
    if (f == NULL) { printf("Could not open mtx file to write\n"); exit(1); }

    coo *A = csc2coo(A_csc);

    fprintf(f,"%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(f,"%d %d %d\n",A->rowS,A->colS,A->nnz);
    for(int i=0;i<A->nnz;i++)
        fprintf(f,"%d %d 1\n",A->row_coo[i]+1,A->col_coo[i]+1);

    destroyCoo(A);

    fclose(f);
    
    std::cout << "Done successfully" << std::endl;
}

/*
 * --CSR to COO conversion--
 * Works with any CSR input matrix.
*/
coo *csr2coo(csr *A_csr){

    int rowCounter = 0, elem=0;
    coo *out = initCoo(A_csr->rowS,A_csr->colS,A_csr->nnz);

    while( rowCounter < A_csr->rowS ){
        int dif = A_csr->r_p[rowCounter+1] - A_csr->r_p[rowCounter];
        
        for(int i=0;i<dif;i++,elem++){
            out->row_coo[elem] = rowCounter;
            out->col_coo[elem] = A_csr->col[elem];
        }
        rowCounter++;
    }

    return out;
}

coo *csc2coo(csc *A_csc){

    int colCounter = 0, elem=0;
    coo *out = initCoo(A_csc->rowS,A_csc->colS,A_csc->nnz);

    while( colCounter < A_csc->rowS ){
        int dif = A_csc->col_ptr[colCounter+1] - A_csc->col_ptr[colCounter];
        
        for(int i=0;i<dif;i++,elem++){
            out->row_coo[elem] = A_csc->row[elem]; //colCounter;
            out->col_coo[elem] = colCounter;       //A_csc->row[elem];
        }
        colCounter++;
    }

    return out;
}

/*
 * --COO to CSR conversion--
 * Works with any COO input matrix.
*/
csr *coo2csr(coo *in){
    int n = in->rowS;
    int nnz = in->nnz;

    csr *out = (csr*)malloc(sizeof(csr));
    out->r_p = (int*)malloc( n * sizeof(int) );
    out->col = (int*)malloc( nnz * sizeof(int) );
    out->nnz = nnz;
    out->rowS = n; out->colS = n;


    //Coo to Csr
    for (int l = 0; l < n+1; l++)
        out->r_p[l] = 0;

    // ----- find the correct column sizes
    for (int l = 0; l < nnz; l++)
        out->r_p[in->row_coo[l] - isOneBased]++;

    // ----- cumulative sum
    for (int i = 0, cumsum = 0; i < n; i++) {
        int temp = out->r_p[i];
        out->r_p[i] = cumsum;
        cumsum += temp;
    }
    
    out->r_p[n] = in->nnz;
    // ----- copy the row indices to the correct place
    for (int l = 0; l < nnz; l++) {
        int row_l;
        row_l = in->row_coo[l] - isOneBased;

        int dst = out->r_p[row_l];
        out->col[dst] = in->col_coo[l] - isOneBased;

        out->r_p[row_l]++;
    }

    // ----- revert the column pointers
    for (int i = 0, last = 0; i < n; i++) {
        int temp = out->r_p[i];
        out->r_p[i] = last;
        last = temp;
    }

    return out;
}

/*
 * COO to CSC conversion
 * Working at least on mtx COO that is sorted by column index 
 * and for each distinct column index the entries are 
 * sorted by row index (The default MATLAB save matrix
 * function meets this requirement)
*/
csc *coo2csc(coo *in){
    int n = in->colS;
    int nnz = in->nnz;

    csc *out = (csc*)malloc(sizeof(csc));
    out->col_ptr = (int*)malloc( n * sizeof(int) );
    out->row = (int*)malloc( nnz * sizeof(int) );
    out->nnz = nnz;
    out->rowS = in->rowS; out->colS = in->colS;


    /*Coo to Csc*/
    for (int l = 0; l < n+1; l++)
        out->col_ptr[l] = 0;

    // ----- find the correct column sizes
    for (int l = 0; l < nnz; l++)
        out->col_ptr[in->col_coo[l] - isOneBased]++;

    // ----- cumulative sum
    for (int i = 0, cumsum = 0; i < n; i++) {
        int temp = out->col_ptr[i];
        out->col_ptr[i] = cumsum;
        cumsum += temp;
    }
    
    out->col_ptr[n] = in->nnz;
    // ----- copy the row indices to the correct place
    for (int l = 0; l < nnz; l++) {
        int col_l;
        col_l = in->col_coo[l] - isOneBased;

        int dst = out->col_ptr[col_l];
        out->row[dst] = in->row_coo[l] - isOneBased;

        out->col_ptr[col_l]++;
    }

    // ----- revert the column pointers
    for (int i = 0, last = 0; i < n; i++) {
        int temp = out->col_ptr[i];
        out->col_ptr[i] = last;
        last = temp;
    }

    return out;
}

/*Transform a 2D n*n symmetric dense matrix A
 * into an similar one "a" in csr format,
 * with columns for each node in ascending order
 * */
csr *Dense2Csr(int *A, int n/*n = nodes*/){
    
    csr *ret = (csr*)malloc(sizeof(csr));
    int nnz=0;
    /*Find non zero values*/
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(A[ind(i,j,n)] == 1)
                nnz++;
        }
    }

    ret->r_p = (int *)malloc((n+1)*sizeof(int));
    ret->col = (int *)malloc(nnz*sizeof(int));
    ret->rowS = n; ret->colS = n; ret->nnz = nnz;

    nnz = 0;
    ret->r_p[0]=0;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(A[ind(i,j,n)] == 1){
                ret->col[nnz] = j;
                nnz++;
            }

        }
        ret->r_p[i+1] = nnz;
    }
    ret->r_p[n] = nnz;

    return ret;
}

void *parse_data(char * fname, int type) {
    coo *M = readCoo(fname);

    if( type == COO )
        return (void*)M;
    
    if (type == CSR){
        csr *Mcsr = coo2csr(M);
        destroyCoo(M);
        return (void*)Mcsr;
    }
    
    if (type == CSC){
        csc *Mcsc = coo2csc(M);
        destroyCoo(M);
        return (void*)Mcsc;
    }
}

coo *readCoo(char filename[]){
  printf("Reading COO...");
    MM_typecode matcode;
    FILE *f;
    int M, N, nz ,ret_code;
    int i, *I, *J;
    coo *ret = (coo*)malloc(sizeof(coo));


    if ((f = fopen(filename, "r")) == NULL)
            exit(1);

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    /* reserve memory for matrices */

    I = (int *) malloc( nz * sizeof(int) );
    J = (int *) malloc( nz * sizeof(int) );


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    /* Replace missing val column with 1s and change the fscanf to match patter matrices*/

    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d 1\n", &I[i], &J[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
        
    }



    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    ret->row_coo = I;
    ret->col_coo = J;
    ret->rowS = N;
    ret->colS = N;
    ret->nnz = nz;


    printf(" finished\n");
	return ret;
}

void destroyCoo(coo *a){
    free(a->col_coo);
    free(a->row_coo);
    free(a);
}
void destroyCsr(csr *a){
    free(a->r_p);
    free(a->col);
    free(a);
}
void destroyCsc(csc *a){
    free(a->row);
    free(a->col_ptr);
    free(a);
}

coo *initCoo(int rowS,int colS, int nnz){
    coo *out = (coo*)malloc(sizeof(coo));
    out->row_coo = (int*)malloc( nnz * sizeof(int) );
    out->col_coo = (int*)malloc( nnz * sizeof(int) );
    out->rowS = rowS;
    out->colS = colS;
    out->nnz = nnz;

    return out;
}

csr *initCsr(int rowS,int colS, int nnz){
  csr *t = (csr*)malloc(sizeof(csr));
  t->col = (int *)malloc(nnz * sizeof(int));
  t->r_p = (int *)malloc((rowS+1) * sizeof(int));
  t->rowS = rowS;
  t->colS = colS;
  t->nnz = nnz;

  return t;
}

csc *initCsc(int rowS,int colS, int nnz){
  csc *t = (csc*)malloc(sizeof(csc));
  t->col_ptr = (int *)calloc( (colS+1), sizeof(int) );
  t->row = (int *)malloc(nnz * sizeof(int));
  t->rowS = rowS;
  t->colS = colS;
  t->nnz = nnz;

  return t;
}

/*
 * Prints a boolean CSR matrix A.
 * Working
 * */
void printCsr(csr *a){
  printf("col = ["); int i; int r_n = 1; int n=a->rowS;
  for(i=0;i<a->nnz;i++){
      if(i == a->r_p[r_n]){
          printf(",  ");
          r_n++;
      }
    printf("%d ",a->col[i]);
  }
  printf("]\n");

  printf("row = [");
  for(i=0;i<n+1;i++){
    printf("%d ",a->r_p[i]);
  }
  printf("]\n");
}

/*
 * Prints a boolean CSC matrix A.
 * Working
 * */
void printCsc(csc *a){
  printf("row = ["); int i; int r_n = 1; int n=a->rowS;
  for(i=0;i<a->nnz;i++){
      if(i == a->col_ptr[r_n]){
          printf(",  ");
          r_n++;
      }
    printf("%d ",a->row[i]);
  }
  printf("]\n");

  printf("col = [");
  for(i=0;i<n+1;i++){
    printf("%d ",a->col_ptr[i]);
  }
  printf("]\n");
}

/*
 * Prints a boolean COO matrix A.
 * Working
 * */
void printCoo(coo *a){
  printf("col = ["); int i;
  for(i=0;i<a->nnz;i++)
    printf("%d ",a->col_coo[i]);
  printf("]\n");

  printf("row = [");
  for(i=0;i<a->nnz;i++)
    printf("%d ",a->row_coo[i]);
  printf("]\n");
}

//End Of File