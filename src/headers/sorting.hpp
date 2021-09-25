#ifndef SORTING_HPP
#define SORTING_HPP

template <typename T>
inline void swap(T *a, T *b){
    T temp = *a;
    *a = *b;
    *b = temp;
}


template <typename T>
inline void merge(T *v, int l, int m, int r)
{
    int n1 = m - l + 1;
    int n2 = r - m;
 
    // Create temp arrays
    T v1[n1], v2[n2];
 
    // Copy data to temp arrays L[] and R[]
    for (int i = 0; i < n1; i++)
        v1[i] = v[l + i];
    for (int j = 0; j < n2; j++)
        v2[j] = v[m + 1 + j];
 
    // Merge the temp arrays back into arr[l..r]
 
    // Initial index of first subarray
    int i = 0;
 
    // Initial index of second subarray
    int j = 0;
 
    // Initial index of merged subarray
    int k = l;
 
    while (i < n1 && j < n2) {
        if (v1[i] <= v2[j]) 
            v[k++] = v1[i++];
        else
            v[k++] = v2[j++];
    }
 
    // Copy the remaining elements of
    // L[], if there are any
    while (i < n1) {
        v[k++] = v1[i++];
    }
 
    // Copy the remaining elements of
    // R[], if there are any
    while (j < n2) {
        v[k++] = v2[j++];
    }
}
 
// l is for left index and r is
// right index of the sub-array
// of arr to be sorted */
template <typename T>
inline void mergesort(T arr[],int l,int r){
    if(l>=r) return;//returns recursively

    int m =l+ (r-l)/2;
    mergesort(arr,l,m);
    mergesort(arr,m+1,r);
    merge(arr,l,m,r);
}

/*
 * -------------------------------------
 * # Bubble Sort Implemention, source : https://www.geeksforgeeks.org/bubble-sort/
 * -------------------------------------
 */
template <typename T>
void bubblesort(T a[], int l, int r)
{
   T* arr = a + l;
   int n = r-l+1;

   int i, j;
   bool swapped;
   for (i = 0; i < n-1; i++)
   {
     swapped = false;
     for (j = 0; j < n-i-1; j++)
     {
        if (arr[j] > arr[j+1])
        {
           swap<T>(&arr[j], &arr[j+1]);
           swapped = true;
        }
     }
  
     // IF no two elements were swapped by inner loop, then break
     if (swapped == false)
        break;
   }
}

/* Takes as input a sorted int array (src) of size and fills dest array
 * with unique values of src (also sorted).
 * Return value is the number of unique elements (<=n) in src 
 */
template <typename T>
inline int unique(T *dest, T *src, int n){
    if(n==0) return 0;
    
    if(n==1) {
        dest[0] = src[0];
        return 1;
    }
    int un = 1;
    dest[0] = src[0];
    for(int i=1;i<n;i++){
        if(src[i] == src[i-1]) continue;
        
        dest[un++] = src[i];
    }
    
    return un;
}

#endif