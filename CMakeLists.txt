cmake_minimum_required(VERSION 3.10)
project(Exercise-4)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-w")

# General Header Files needed for all Versions
set(src/headers src/headers/mmio.cpp src/headers/my_time.cpp src/headers/file_io.cpp)
set(sparse_headers src/sparse/sparse.cpp src/sparse/matrix.cpp)

# Build Sequential Version
add_executable(seq ${headers} ${sparse_headers} src/sequential/sequential1.cpp src/sparse/bmm_blocking.cpp)

# Build Multithreaded Version
find_package(OpenMP REQUIRED)
add_executable(src/multithreaded ${headers} ${sparse_headers} src/multithreaded/omp.cpp src/sparse/bmm_blocking.cpp src/multithreaded/bmm_multithread.cpp)
target_link_libraries(src/multithreaded PRIVATE "${OpenMP_CXX_FLAGS}")
target_compile_options(src/multithreaded PRIVATE "${OpenMP_CXX_FLAGS}")
