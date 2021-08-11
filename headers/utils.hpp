#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include "file_io.h"
#include <string.h>

#define SEQUENTIAL 1
#define DISTRIBUTED 2

template <typename T>
inline void printV(T *vec, int n){
  std::cout << "[";
  for(int i=0;i<n-1;i++)
    std::cout << vec[i] << " ";
  
  std::cout << vec[n-1] << "]" << std::endl;
}

inline void do_some_work(){
    
    int n = 1000000;
    float A[n];

    for(int i=1;i<n;i++) A[i] = 0.5*A[i-1];
    printf("Job - Done - Last is %f\n",A[n-1]);
}

inline void print_tuple(int a, int b) { printf("(%d,%d)",a,b); }
inline void println_tuple(int a, int b) { printf("(%d,%d)\n",a,b); }

inline void write_times(int N,struct timespec t1, struct timespec t2, int TYPE, int proc_num){
  char type[20];
  char msg[70];
  char fname[] = "./logs/times.csv";

  if(TYPE == SEQUENTIAL)
    strcpy(type,"Sequential");
  if(TYPE == DISTRIBUTED)
    strcpy(type,"Distributed");

  long seconds = t2.tv_sec - t1.tv_sec;
  long nanoseconds = t2.tv_nsec - t1.tv_nsec;
  if(t2.tv_nsec<t1.tv_nsec){
    seconds--;
    nanoseconds +=1000000000;
  }

  double total = ((double)seconds)*(1e9) + ((double)nanoseconds);
  total /= (double) 1e9;

  // printf("Time elapsed %ld seconds %ld nanoSeconds (%.7f sec)\n", seconds, nanoseconds,total);

  sprintf(msg, "%d,%.7f,%s,%d", N, total, type, proc_num );

  file_writer(fname, msg, APPEND);
}
inline void clear_times(){
  char fname[] = "./logs/times.csv";
  file_writer(fname, "N,time_sec,type, procs", OVERWRITE);
}


#endif
