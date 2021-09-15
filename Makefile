MPI_PROCS=4

# Define Appropriate Compilers
CC=g++
MPICC=mpiCC

# Define Various Flags used on builds
OPTIMIZATION=-O3
FLAGS=-w
DEBUG=-g
HEADERS=headers/mmio.c headers/my_time.cpp headers/file_io.c
SPARSE_HEADERS=./sparse/sparse.cpp ./sparse/matrix.cpp

all: sequential distributed multithreaded


# ----- Matlab/Octave I/O
data:
		bash ./bash_src/set_parameters.sh $(N)
		date "+%H:%M:%S   %d/%m/%y" > ./logs/generate_data.txt
		echo "----Generating Data----" >> ./logs/generate_data.txt
		cd ./matlab && octave generateDatasets.m >> ../logs/generate_data.txt
		echo "----Data Generated Successfully" >> ./logs/generate_data.txt

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

ultra_fast_test_distributed: 
		mpirun -np $(MPI_PROCS) ./distributed.out

ultra_fast_test_multithreaded: 
		./multithreaded.out ${THREAD_NUM}


# ----- Build Project
sequential:
		$(CC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) -o seq.out sequential1.cpp bmm_blocking.cpp $(FLAGS)

distributed:
		$(MPICC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) ./headers/my_mpi_comms.cpp -o distributed.out distributed.cpp bmm_distributed.cpp bmm_blocking.cpp $(FLAGS)
multithreaded:
		$(CC) $(OPTIMIZATION) $(HEADERS) $(SPARSE_HEADERS) -o multithreaded.out omp.cpp bmm_blocking.cpp bmm_multithread.cpp -fopenmp $(FLAGS)


# ----- Measure Performance
performance:
		cd ./matlab && python main.py && cd ..




# ----- Clean/Purge
clear_times:
		rm ./logs/times.csv && echo "dim,time,type,procs" >> "./logs/times.csv"
clean: 
		echo "Removing executable.."
		rm ./*.out
		echo "Done Successfully"
		echo "Removing Datasets.."
		rm ./datasets/test/*_test.mtx ./datasets/test/*result.mtx
		echo "Done Successfully"
