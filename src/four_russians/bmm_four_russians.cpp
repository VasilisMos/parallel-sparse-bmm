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
    csr *Ci = initCsr(A->rowS, B->colS, 3*(nnz(A) + nnz(B)));
    csr *Rs = initCsr((int)pow(2,m), n, 3*(nnz(A) + nnz(B)));

    zero_padding(A,n);
    zero_padding(B_csr,n);

    for(int i=1;i<=i_max;i++){
        int bp=1;
        int k=0;
        int block_size = m;

        csc *Ai = get_vertical_chunk(A,chunk_start, chunk_end);
        csr *A_csr = csc2csr(Ai);

        csr *Bi = (csr*)malloc(sizeof(csr));
        Bi->rowS = m; Bi->colS = B_csr->colS;
        Bi->r_p= B_csr->r_p + chunk_start;
        Bi->col = B_csr->col + B_csr->r_p[chunk_start];
        Bi->nnz = B_csr->r_p[chunk_end] - B_csr->r_p[chunk_start];

        printCsr(Bi);
        
        // Rs[0] = ZERO_ROW
        Rs->r_p[0]=0; Rs->r_p[1]=0;
        Rs->col[0]=0;

        for(int j=1;j<pow(2,m);j++){
            // Rs(j,:) = or(Rs(j-2^k,:),Bi(m+1 - (k+1),:));
            or_rows(Bi,m-(k+1),    Rs,j-pow(2,k),   j,chunk_start);

            if(bp==1){
                bp=j+1;
                k++;
            }
            else 
                bp--;
        }

        free(Bi);

        Ci->r_p[0] = 0;
        int nnz_c = 0;
        for(int j=0;j<n;j++){
            
            // ind = num(Ai(j,:))
            int s = A_csr->r_p[j];
            int f = A_csr->r_p[j+1];
            int ind = num(A_csr->col + (j+chunk_start),f-s,chunk_start);
            ind = 1;

            //Ci(j,:) = Rs(ind,:);
            s = Rs->r_p[ind];
            f = Rs->r_p[ind+1];

            for(int k=s;k<f;k++){
                Ci->col[nnz_c++] = Rs->col[k];
//                cout << Rs->col[k] << endl;
            }
//            cout << endl;
            Ci->r_p[j+1] = nnz_c;
        }

        destroyCsr(A_csr);
    }

//    printCsr(Ci);

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
}

/*
 * num(): Sum Of Powers of 2, Specified by Input Vector Ai
 * Example Ai={0,3,5}
 * result: 2^0+2^3+2^5=41
 */
int num(int *Ai, int n, int offset){
    int logic = Ai[0]-offset;
    int res = (logic == 0) ? 1 : 0;
    int i_0 = (logic == 0) ? 1 : 0;
    for(int i=i_0;i<n;i++){
        res+= 2<<(Ai[i]-1-offset);
    }

    return res;
}

csc *get_vertical_chunk(csc *A, int start, int end){
    int n = A->rowS;
    int m = (int) floor(log2(n));

    int nnz_chnk = A->col_ptr[end] - A->col_ptr[start];
    int offs = A->col_ptr[start];

    csc *Ai = initCsc(n,m, nnz_chnk);

    for(int c=0;c<m+1;c++)
        Ai->col_ptr[c] = A->col_ptr[start+c] - offs;

    for(int r=0;r<nnz_chnk;r++)
        Ai->row[r] = A->row[r+offs];


    return Ai;
}