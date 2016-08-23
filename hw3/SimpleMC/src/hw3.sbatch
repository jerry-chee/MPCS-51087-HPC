#!/bin/bash
#SBATCH --job-name=mpcs51087_hw3
#SBATCH --output=hw3.out
#SBATCH --error=hw3.err
#SBATCH --time=00:10:00
#SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9
#SBATCH --open-mode=append|truncate

module load openmpi

mpirun -n 9 ./transport -particles 10000000 -bins 16
