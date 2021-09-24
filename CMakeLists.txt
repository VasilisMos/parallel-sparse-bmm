cmake_minimum_required(VERSION 3.10)
project(Exercise-4)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-w")

# General Header Files needed for all Versions
set(headers headers/mmio.cpp headers/my_time.cpp headers/file_io.cpp)
set(sparse_headers ./sparse/sparse.cpp ./sparse/matrix.cpp)

# Build Sequential Version
add_executable(seq ${headers} ${sparse_headers} sequential1.cpp bmm_blocking.cpp )

# Build Multithreaded Version
find_package(OpenMP REQUIRED)
add_executable(multithreaded ${headers} ${sparse_headers} omp.cpp bmm_blocking.cpp bmm_multithread.cpp)
target_link_libraries(multithreaded PRIVATE "${OpenMP_CXX_FLAGS}")
target_compile_options(multithreaded PRIVATE "${OpenMP_CXX_FLAGS}")
