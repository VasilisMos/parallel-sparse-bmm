import performance

# filename = './times.csv'
filename = '../../logs/times.csv'
df_dist = performance.get_version_times(filename, 'Distributed')
# df_seq = performance.get_version_times(filename, 'Sequential').loc[:, 'time']
# df_matlab = performance.get_version_times(filename, 'MATLAB').loc[:, 'time']
# df_omp = performance.get_version_times(filename, 'Multithreaded').loc[:, 'time']
# df_hybrid = performance.get_version_times(filename, 'Hybrid').loc[:, 'time']
# df_matlab = performance.get_version_times(filename, 'MATLAB').loc[:, 'time']

df_seq = performance.get_version_times(filename, 'Sequential')
df_matlab = performance.get_version_times(filename, 'MATLAB')
df_omp = performance.get_version_times(filename, 'Multithreaded')
df_hybrid = performance.get_version_times(filename, 'Hybrid')
df_matlab = performance.get_version_times(filename, 'MATLAB')

# performance.test_exec_times(df_matlab)
# performance.test_mpi_performance(filename=filename)
performance.test_exec_times(df_matlab, df_seq, df_dist, df_omp, df_hybrid)
