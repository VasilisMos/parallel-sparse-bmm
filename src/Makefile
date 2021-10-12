MPI_PROCS=4

# Define Appropriate Compilers
CC=g++
MPICC=mpiCC

# Define Various Flags used on builds
OPTIMIZATION=-O3
FLAGS=-w
DEBUG=-g

HEADERS=headers/mmio.cpp headers/my_time.cpp headers/file_io.cpp
SPARSE_HEADERS=./sparse/sparse.cpp ./sparse/matrix.cpp ./sparse/bmm_blocking.cpp ./sparse/masked_bmm.cpp
SEQ_FILES=./sequential/sequential1.cpp
MULTITHREADED_FILES=./multithreaded/bmm_multithread.cpp ./multithreaded/omp.cpp ./multithreaded/masked_bmm_multithreaded.cpp
DISTRIBUTED_FILES=./distributed/distributed.cpp ./distributed/bmm_distributed.cpp ./headers/my_mpi_comms.cpp
DISTRIBUTED_FILES2=./distributed2/distributed2.cpp ./distributed2/bmm_distributed2.cpp ./headers/my_mpi_comms.cpp
HYBRID_FILES=./hybrid/bmm_hybrid.cpp ./headers/my_mpi_comms.cpp

all: sequential distributed multithreaded 

.PHONY: sequential distributed multithreaded four_russians masked distributed2 hybrid

# ----- Matlab/Octave I/O
data:
		bash ./bash_src/set_parameters.sh $(N)
		date "+%H:%M:%S   %d/%m/%y" > ../logs/generate_data.txt
		echo "----Generating Data----" >> ../logs/generate_data.txt
		cd ./matlab && octave generateDatasets.m >> ../../logs/generate_data.txt
		echo "----Data Generated Successfully" >> ../logs/generate_data.txt

filter:
		bash ./bash_src/set_filter_parameters.sh $(N) $(D)
		date "+%H:%M:%S   %d/%m/%y" > ../logs/generate_data.txt
		echo "----Generating Filter Data----" >> ../logs/generate_data.txt
		cd ./matlab && octave generate_filter.m >> ../../logs/generate_data.txt
		echo "----Filter Data Generated Successfully" >> ../logs/generate_data.txt

validate_results:
		cd ./matlab && octave getResults.m



# ----- Execute a Full Application Test
test: data sequential
		clear && ./seq.out && rm ./seq.out

test_distributed: data distributed
		clear && mpirun -np $(MPI_PROCS) ./distributed.out && rm ./distributed.out
test_multithreaded: data multithreaded
		clear && ./multithreaded.out ${THREAD_NUM} && rm ./multithreaded.out

test_masked: filter masked
		clear && ./masked.out && rm ./masked.out

fast_test: sequential
		./seq.out && rm ./seq.out

fast_test_distributed: distributed
		mpirun -np $(MPI_PROCS) ./distributed.out &&	rm ./distributed.out

fast_test_four_russians: four_russians
		./four_russ.out && rm ./four_russ.out

ultra_fast_test_sequential:
		./seq.out 

ultra_fast_test_distributed:
		mpirun -np $(MPI_PROCS) ./distributed.out

ultra_fast_test_distributed2:
		mpirun -np $(MPI_PROCS) ./distributed2.out

ultra_fast_test_hybrid:
		export OMP_NUM_THREADS=1
		mpirun -np $(MPI_PROCS) ./hybrid.out

ultra_fast_test_multithreaded:
		./multithreaded.out ${THREAD_NUM}


# ----- Build Project
sequential:
		$(CC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) -o seq.out ./sequential/sequential1.cpp $(FLAGS)

multithreaded:
		$(CC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) -o multithreaded.out $(MULTITHREADED_FILES) -fopenmp $(FLAGS)

distributed:
		$(MPICC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) -o distributed.out $(DISTRIBUTED_FILES) $(FLAGS)

distributed2:
		$(MPICC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) -o distributed2.out $(DISTRIBUTED_FILES2) $(FLAGS)

four_russians:
		$(CC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) -o four_russ.out ./four_russians/bmm_four_russians.cpp ./four_russians/main.cpp ./four_russians/four_russians_utils.cpp $(FLAGS)

masked:
		$(CC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) -o masked.out ./sequential/sequential_masked.cpp $(FLAGS)

hybrid:
		$(MPICC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) -o hybrid.out ./hybrid/hybrid.cpp $(HYBRID_FILES) $(FLAGS) -fopenmp



# ----- Measure Performance
performance:
		cd ./matlab && python main.py && cd ..




# ----- Clean/Purge
clear_times:
		rm ../logs/times.csv && echo "dim,time,type,procs" >> "../logs/times.csv"
clean:
		echo "Removing executable.."
		rm ./*.out
		echo "Done Successfully"
		echo "Removing Datasets.."
		rm ../datasets/test/*_test.mtx ../datasets/test/*result.mtx
		echo "Done Successfully"
		