MPI_PROCS=4

# Define Appropriate Compilers
CC=g++
MPICC=mpiCC

# Define Various Flags used on builds
OPTIMIZATION=-O3
FLAGS=-w
DEBUG=-g

HEADERS=headers/mmio.cpp headers/my_time.cpp headers/file_io.cpp
SPARSE_HEADERS=./sparse/sparse.cpp ./sparse/matrix.cpp ./sparse/bmm_blocking.cpp
SEQ_FILES=./sequential/sequential1.cpp
MULTITHREADED_FILES=./multithreaded/bmm_multithread.cpp ./multithreaded/omp.cpp
DISTRIBUTED_FILES=./distributed/distributed.cpp ./distributed/bmm_distributed.cpp ./headers/my_mpi_comms.cpp

all: sequential distributed multithreaded

.PHONY: sequential distributed multithreaded four_russians

# ----- Matlab/Octave I/O
data:
		bash ./bash_src/set_parameters.sh $(N)
		date "+%H:%M:%S   %d/%m/%y" > ../logs/generate_data.txt
		echo "----Generating Data----" >> ../logs/generate_data.txt
		cd ./matlab && octave generateDatasets.m >> ../../logs/generate_data.txt
		echo "----Data Generated Successfully" >> ../logs/generate_data.txt

validate_results:
		cd ./matlab && octave getResults.m



# ----- Execute a Full Application Test
test: data sequential
		clear && ./seq.out && rm ./seq.out

test_distributed: data distributed
		clear && mpirun -np $(MPI_PROCS) ./distributed.out && rm ./distributed.out
test_multithreaded: data multithreaded
		clear && ./multithreaded.out ${THREAD_NUM} && rm ./multithreaded.out
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

ultra_fast_test_multithreaded:
		./multithreaded.out ${THREAD_NUM}


# ----- Build Project
sequential:
		$(CC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) -o seq.out ./sequential/sequential1.cpp $(FLAGS)

multithreaded:
		$(CC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) -o multithreaded.out $(MULTITHREADED_FILES) -fopenmp $(FLAGS)

distributed:
		$(MPICC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) -o distributed.out $(DISTRIBUTED_FILES) $(FLAGS)

four_russians:
		$(CC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) -o four_russ.out ./four_russians/bmm_four_russians.cpp ./four_russians/main.cpp $(FLAGS)

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
		