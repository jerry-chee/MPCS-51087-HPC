#!/bin/bash
#SBATCH --job-name=mpcs51087_hw5_openmp
#SBATCH --output=hw5_openmp.out
#SBATCH --error=hw5_openmp.err
#SBATCH --time=00:10:00
#SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --open-mode=append|truncate

N=1 
MAX=10 #3^5+1

export OMP_NUM_THREADS=16
while [	$N -lt $MAX ]; do
	../src/simpson_openmp $N
	let N=N+1
done
