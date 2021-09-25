#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#define TEST
// #define SMALL
#define BLOCKING_FACTOR 4

#ifdef SMALL
#define fname1 "./datasets/small/small.mtx"
#define fname2 "./datasets/small/small.mtx"
#define fname3 "./datasets/small/C_result.mtx"
#endif

#ifdef TEST
#define fname1 "../datasets/test/A_test.mtx"
#define fname2 "../datasets/test/B_test.mtx"
#define fname3 "../datasets/test/C_result.mtx"
#endif


#endif
