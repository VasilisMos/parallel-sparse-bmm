cmake_minimum_required(VERSION 3.10)
project(Exercise-4)

set(headers headers/mmio.cpp headers/my_time.cpp headers/file_io.cpp)
set(sparse_headers ./sparse/sparse.cpp ./sparse/matrix.cpp)

add_executable(seq ${headers} ${sparse_headers} sequential1.cpp bmm_blocking.cpp )
add_executable(multithreaded ${headers} ${sparse_headers} omp.cpp bmm_blocking.cpp bmm_multithread.cpp)