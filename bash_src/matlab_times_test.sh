#!/bin/bash
#SBATCH -J Matlab-Job
#SBATCH -p batch
#SBATCH -t 02:00

#### load module ####
module load matlab/R2021a

### Go to matlab-src directory ###
cd ../matlab

#### RUN MATLAB JOB ####
echo "Starting Matlab"
matlab -nodisplay << EOF &> matlab.out

disp("Sanity Check");
a = 15;
disp(a)

timing_performance;

exit
EOF