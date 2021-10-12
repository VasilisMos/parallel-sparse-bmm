#ifndef MASKED_BMM_HPP
#define MASKED_BMM_HPP

#include "matrix.hpp"
#include "../headers/parameters.hpp"
#include "bmm_blocking.hpp"


void masked_bmm(csc *A, csc *B, csc *C,csc *F);
void masked_bmm_blocking(csc *A, csc *B, csc *C,csc *F);

#endif