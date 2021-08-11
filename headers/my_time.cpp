#include "my_time.hpp"

struct timespec t1,t2;

void time_elapsed(struct timespec t2, struct timespec t1){
  long seconds = t2.tv_sec - t1.tv_sec;
  long nanoseconds = t2.tv_nsec - t1.tv_nsec;
  if(t2.tv_nsec<t1.tv_nsec){
    seconds--;
    nanoseconds +=1000000000;
  }

  double total = ((double)seconds)*(1e9) + ((double)nanoseconds);
  total /= (double) 1e9;

  printf("Time elapsed %ld seconds %ld nanoSeconds (%.7f sec)\n", seconds, nanoseconds,total);
}

struct timespec tic() { clock_gettime(CLOCK_MONOTONIC,&t1); return t1; }
struct timespec toc() { clock_gettime(CLOCK_MONOTONIC,&t2); /*time_elapsed(t2,t1);*/ return t2; }
