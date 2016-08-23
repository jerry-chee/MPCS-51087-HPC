#!/bin/bash
#SBATCH --job-name=mpcs51087_hw2
#SBATCH --output=hw2.out
#SBATCH --error=hw2.err
#SBATCH --time=00:10:00
#SBATCH --partition=sandyb
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2

module load openmpi

N=10000
K=100
T=1
MAX=5

NPROCS=1
while [	$NPROCS -lt $MAX ]; do
	mpirun -n $NPROCS ../src/mandelbrot $N $K $T 
	let NPROCS=NPROCS*2
done
