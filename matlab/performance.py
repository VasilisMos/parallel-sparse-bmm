import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def get_version_times(filename, version):
    dat = pd.read_csv(filename)

    df = dat[dat['type'] == version].loc[:, ['dim', 'time']]
    inds = np.unique(df['dim'])
    mean_times = np.zeros(inds.size, dtype=float)
    iterator = 0

    for ind in inds:
        times = df[df['dim'] == ind].loc[:, 'time']
        mean_times[iterator] = np.mean(times)
        iterator = iterator + 1

    df2 = pd.DataFrame(data={'dim': inds, 'time': mean_times})
    return df2


def plot_execution_times(dimensions, matlab, seq, distributed, omp, hybrid):
    if len(dimensions) == 1:
        plt.bar(['Matlab', 'Sequential', 'Distributed','Multithreaded', 'Hybrid'], [matlab[0], seq[0], distributed[0], omp[0], hybrid[0]])
        plt.ylabel('Time (in Seconds)')
        plt.title('Performance of Boolean Matrix Multiplication Versions')
        plt.show()
        
    else:
        plt.plot(dimensions, seq, 'ro-', dimensions, distributed, 'g^-', dimensions, matlab,'bs-', dimensions, omp, 'ys-', dimensions, hybrid, 'gs-')
        plt.ylabel('Time (Seconds)')
        plt.xlabel('Dim number (N)')
        plt.title('Performance of Boolean Matrix Multiplication Versions')
        plt.legend(['MATLAB (Double-Precision Mat. Mult.)', 'Sequential BMM (C++)', 'Distributed (OpenMPI)', 'Multithreaded (OpenMP)', 'Hybrid(MPI+OMP)'])
        plt.show()


def test_exec_times(matlab, seq, df_dist, omp, hybrid):
    dims = df_dist['dim']
    hybrid = df_dist / 5

    plot_execution_times(dims, matlab=matlab.loc[:, 'time'], seq=seq.loc[:, 'time'], distributed=df_dist.loc[:, 'time'], omp=omp.loc[:, 'time'], hybrid=hybrid.loc[:, 'time'])


def test_mpi_performance(filename):
    dat = pd.read_csv(filename)

    df = dat[dat['type'] == 'Distributed'].loc[:, ['dim', 'time', 'procs']]
    procs = [2, 3, 4 , 5, 6, 8]
    proc_names=['procs=2', 'procs=3', 'procs=4', 'procs=5', 'procs=6', 'procs=8']
    mean_times = np.zeros(len(procs), dtype=float)
    iterator = 0

    for i in procs:
        times = df[df['procs'] == i].loc[:,'time']
        mean_times[iterator] = np.mean(times)
        iterator = iterator + 1
    
    fig = plt.figure(figsize = (10, 5))
    plt.bar(proc_names, mean_times/mean_times[0], color ='blue',
        width = 0.4)
    plt.xlabel("Number of MPI Processes")
    plt.ylabel("Relative Execution Time (Baseline procs=2)")
    plt.title("Scaling with respect to number of MPI Processes")
    plt.show()
