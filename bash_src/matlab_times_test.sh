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
matlab -nodisplay << EOF &> ../logs/matlab.out

disp("Sanity Check");
a = 15;
disp(a)

n = (1:4) * 1.2e4;

for i=1:size(n)
    generateDatasets(1.2e4);
end
exit
EOF