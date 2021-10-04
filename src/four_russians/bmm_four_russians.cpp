#include "bmm_four_russians.hpp"
#define chunk_start ((i-1)*block_size)
#define chunk_end (i*block_size-1)

csr *csc2csr(csc *A){
    int rowS = A->rowS, colS = A->colS;
    int nnz = A->nnz;
    int cumsum=0;
    csr *B = initCsr(rowS,colS, nnz);
    B->nnz = nnz;

    for(int i=0;i<rowS+1;i++)
        B->r_p[i] = 0;

    for(int i=0;i<nnz;i++)
        B->r_p[A->row[i]]++;

    for(int i=0;i<rowS;i++){
        int temp  = B->r_p[i];
        B->r_p[i] = cumsum;
        cumsum   += temp;
    }

    for (int col = 0; col < colS; col++)
    {
        for(int j=A->col_ptr[col]; j<A->col_ptr[col+1]; j++){
            int row = A->row[j];
            int dest = B->r_p[row];

            B->col[dest] = col;
            B->r_p[row]++;
        }
    }

    B->r_p[rowS] = nnz;

    for(int row = 0, last = 0; row <= rowS; row++){
        int temp    = B->r_p[row];
        B->r_p[row] = last;
        last        = temp;
    }
    
    return B;
}

csc *csr2csc(csr *A){
    const int rowS = A->rowS;
    const int colS = A->colS;
    const int nnz = A->r_p[rowS];
    csc *B = initCsc(rowS, colS, nnz);

    for(int i=0;i<colS+1;i++)
        B->col_ptr[i] = 0;

    for (int n = 0; n < nnz; n++)       
        B->col_ptr[A->col[n]]++;

    //cumsum the nnz per column to get Bp[]
    for(int col = 0, cumsum = 0; col < colS; col++){     
        int temp  = B->col_ptr[col];
        B->col_ptr[col] = cumsum;
        cumsum += temp;
    }
    B->col_ptr[colS] = nnz; 

    for(int row = 0; row < rowS; row++){
        for(int jj = A->r_p[row]; jj < A->r_p[row+1]; jj++){
            int col  = A->col[jj];
            int dest = B->col_ptr[col];

            B->row[dest] = row;

            B->col_ptr[col]++;
        }
    }  

    for(int col = 0, last = 0; col <= colS; col++){
        int temp  = B->col_ptr[col];
        B->col_ptr[col] = last;
        last    = temp;
    }

    return B;
}

csc *bmm_four_russians(csc *A, csc *B,csc *C){
    int n = A->rowS;
    int m = (int) floor(log2(n));
    int i_max = (int)ceil(n/m);

    csr *B_csr = csc2csr(B);
    csr *Rs = initCsr((int)pow(2,m), n, 3*(nnz(A) + nnz(B)));

    csr *Ci = initCsr(A->rowS, B->colS, 3*(nnz(A) + nnz(B)));
    csr *C_csr = initCsr(A->rowS, B->colS, 3*(nnz(A) + nnz(B)));
    C_csr->nnz=0; for(int i=0;i<C_csr->rowS+1;i++) C_csr->r_p[i]=0;

    csc *Ai[i_max];
    csr *Bi[i_max],*A_csr[i_max];
    csr *temp[i_max+1];  temp[0] = C_csr;

    zero_padding(A,n);
    zero_padding(B_csr,n);

    for(int i=1;i<=i_max;i++){
        int bp=1, k=0;
        int block_size = m;
//        cout << "Chunk start: " << chunk_start << ", Chunk end:" << chunk_end << endl;

        Ai[i-1] = get_vertical_chunk(A,chunk_start);
        Bi[i-1] = get_horizontal_chunk(B_csr, chunk_start);
        A_csr[i-1] = csc2csr(Ai[i-1]);
        
        // Rs[0] = ZERO_ROW
        Rs->r_p[0]=0; Rs->r_p[1]=0;
        Rs->col[0]=0;

        for(int j=1;j<pow(2,m);j++){
            // Rs(j,:) = or(Rs(j-2^k,:),Bi(m+1 - (k+1),:));
            or_rows(Bi[i-1],m-(k+1),    Rs,j-pow(2,k),   j,0);
            loop_logic(&bp,&j,&k);
        }

        Rs->nnz = Rs->r_p[1<<m];
//        printCsr(Rs);

        int nnz_c = 0; Ci->r_p[0] = 0;
        for(int j=0;j<n;j++){
            
            // ind = num(Ai(j,:))
            int s = A_csr[i-1]->r_p[j];
            int f = A_csr[i-1]->r_p[j+1];
            int ind = num(A_csr[i-1]->col + s,f-s,0);
//            cout << ind << endl;

            //Ci(j,:) = Rs(ind,:);
            s = Rs->r_p[ind];
            f = Rs->r_p[ind+1];

            for(int k=s;k<f;k++){
                Ci->col[nnz_c++] = Rs->col[k]; //cout << Rs->col[k] << endl;
            }
            Ci->r_p[j+1] = nnz_c;
            Ci->nnz = nnz_c;
        }

        cout << endl << endl << "Ci:" << endl; printCsr(Ci);
        temp[i] = or_csr(temp[i-1],Ci);
        cout << "temp:" << endl;                printCsr(temp[i]);

//        destroyCsc(Ai[i-1]); destroyCsr(Bi[i-1]); destroyCsr(A_csr[i-1]);
    }

    for (int i = 0; i < i_max; ++i) {
        destroyCsr(Bi[i]);
        destroyCsr(A_csr[i]);
        destroyCsc(Ai[i]);
    }

    C = csr2csc(C_csr);
//    printCsr(C_csr);

    destroyCsr(Ci); destroyCsr(B_csr); destroyCsr(C_csr); destroyCsr(Rs);
    return C;
}

void or_rows(csr *A, int i1, csr *Rs, int j1, int i_targ, int offset){
    // Assumption length(A(i,:)) == length(Rs(j,:))
    int n = A->colS;
    int i=0,j=0,un=0;

    int size1 = A->r_p[i1+1] - A->r_p[i1];
    int size2 = Rs->r_p[j1+1] - Rs->r_p[j1];

    int *r1 = A->col + (A->r_p[i1]);
    int *r2 = Rs->col + (Rs->r_p[j1]);
    int *targ = Rs->col + (Rs->r_p[i_targ]);

    while(i<size1 && j<size2){
        if(r1[i]-offset < r2[j])
            targ[un++] = r1[i++] - offset;
        else if(r1[i]-offset > r2[j])
            targ[un++] = r2[j++];
        else {
            targ[un++] = r1[i++] - offset;
            j++;
        }
    }

    while(i<size1){
        targ[un++] = r1[i++] - offset;
    }

    while(j<size2){
        targ[un++] = r2[j++];
    }

    Rs->r_p[i_targ+1] = Rs->r_p[i_targ] + un;

//    printf("Bi(s=%d):",size1);
//    printV<int>(r1,size1);
//
//    printf("Rs(j,:)(s=%d):",size2);
//    printV<int>(r2,size2);
//
//    cout << "targ";
//    printV<int>(targ,un); cout << endl << endl;
}

/*
 * num(): Sum Of Powers of 2, Specified by Input Vector Ai
 * Example Ai={0,3,5}
 * result: 2^0+2^3+2^5=41
 */
int num(int *Ai, int n, int offset){
    if (n==0) return 0;
    int logic = Ai[0]-offset;
    int res = (logic == 0) ? 1 : 0;
    int i_0 = (logic == 0) ? 1 : 0;
    for(int i=i_0;i<n;i++){
        res+= 2<<(Ai[i]-1-offset);
    }

//    printf("v(s=%d)=",n);
//    printV<int>(Ai,n);

    return res;
}

csc *get_vertical_chunk(csc *A, int start){
    int n = A->rowS;
    int m = (int) floor(log2(n));

    int nnz_chnk = A->col_ptr[start+m] - A->col_ptr[start];
    int offs = A->col_ptr[start];

    csc *Ai = initCsc(n,m, nnz_chnk);

    for(int c=0;c<m+1;c++)
        Ai->col_ptr[c] = A->col_ptr[start+c] - offs;

    for(int r=0;r<nnz_chnk;r++)
        Ai->row[r] = A->row[r+offs];

    Ai->col_ptr[m] = nnz_chnk;

    return Ai;
}

csr *get_horizontal_chunk(csr *B, int start){
    int n = B->colS;
    int m = (int) floor(log2(n));

    int nnz_chnk = B->r_p[start+m] - B->r_p[start];
    int offs = B->r_p[start];

    csr *Bi = initCsr(m,n,nnz_chnk);

    for(int c=0;c<m+1;c++)
        Bi->r_p[c] = B->r_p[start+c] - offs;

    for(int r=0;r<nnz_chnk;r++)
        Bi->col[r] = B->col[r+offs];

    Bi->r_p[m] = nnz_chnk;

    return Bi;
}

/*
 * Calculates the hadamard OR result of 
 * sparse matrices A and B (in CSR format)
 */
csr *or_csr(csr *A, csr *B){
    if(!((A->colS == B->colS) && (A->rowS == B->rowS))) {
        printf("Or function ERROR: Matrice Dims Mismatch\n");
        exit(1);
    } 

    csr *res = (csr*) initCsr(A->rowS,A->colS,nnz(A) + nnz(B));
    int *v = my_malloc<int>(A->rowS);

    res->r_p[0] = 0;
    int nnz = 0;

    for(int r=0;r<A->rowS;r++){
        int s1 = A->r_p[r];
        int s2 = B->r_p[r];
        int n1 = A->r_p[r+1] - s1;
        int n2 = B->r_p[r+1] - s2;

        int *v1 = A->col + s1;
        int *v2 = B->col + s2;

        int nnz_row = vector_or(v1,n1,v2,n2,v+nnz);
        nnz+=nnz_row;
        res->r_p[r+1] = nnz;
    }

    free(res->col);
    res->col = v;

    return res;
}

int vector_or(int *r1, int size1, int *r2, int size2, int *dest){

    int i=0,j=0,un=0;

    while(i<size1 && j<size2){
        if(r1[i] < r2[j])
            dest[un++] = r1[i++];
        else if(r1[i] > r2[j])
            dest[un++] = r2[j++];
        else {
            dest[un++] = r1[i++];
            j++;
        }
    }

    while(i<size1){
        dest[un++] = r1[i++];
    }

    while(j<size2){
        dest[un++] = r2[j++];
    }

    return un;
}