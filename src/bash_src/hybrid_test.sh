#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=batch
#SBATCH --time=00:01:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4

echo "Srun with 4/4"
srun -np 4 ./hybrid.out 4