#ifndef MY_TIME_HPP
#define MY_TIME_HPP

#include <time.h>
#include <stdio.h>

void time_elapsed(struct timespec t2, struct timespec t1);

struct timespec tic();
struct timespec toc();

#endif