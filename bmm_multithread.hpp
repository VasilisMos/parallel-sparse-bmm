#ifndef BMM_MULTITHREADED_HPP
#define BMM_MULTITHREADED_HPP

#include "bmm_blocking.hpp"
#include <omp.h>

#define TIC tic();
#define TOC toc();

void run_bmm_multithread(int total_procs);
csc *bmm_multithreads(csc **Abl_total, csc **Bbl_total, int total_procs);

#endif