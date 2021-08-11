#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=batch
#SBATCH --time=00:07:00

module purge
module load gcc/7.3.0 cuda/10.0.130 libpng
make sequential
