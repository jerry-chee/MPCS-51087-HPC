#!/bin/bash
#SBATCH --job-name=mpcs51087_hw2
#SBATCH --output=hw2_ltbd.out
#SBATCH --error=hw2_ltbd.err
#SBATCH --time=0:10:00
#SBATCH --partition=sandyb
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1

module load openmpi
#ompi_info --param mpi all

mpirun -n 2 ../src/latency_bandwidth
