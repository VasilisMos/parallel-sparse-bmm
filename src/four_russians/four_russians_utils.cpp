#include "four_russians_utils.hpp"


/*
 * num(): Sum Of Powers of 2, Specified by Input Vector Ai
 * Example Ai={0,3,5}
 * result: 2^0+2^3+2^5=41
 */
int num(int *Ai, int n, int offset,int dim){
    if (n==0) return 0;

    int res = 0; // = (logic == 0) ? 1 : 0;
    for(int i=0;i<n;i++){
        if(Ai[i] < dim-1)
            res+= 2<<(dim - (Ai[i])-2);
        else
            res++;
    }

    return res;
}