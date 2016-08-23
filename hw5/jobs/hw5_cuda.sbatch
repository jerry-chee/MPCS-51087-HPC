#!/bin/bash
#SBATCH --account=mpcs51087
#SBATCH --job-name=mpcs51087_hw5_cuda
#SBATCH --output=hw5_cuda.out
#SBATCH --error=hw5_cuda.err
#SBATCH --time=00:10:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --exclusive
#SBATCCH --open-mode=append|truncate

module load cuda

N=1 
MAX=10 

while [ $N -lt $MAX ]; do
	../src/simpson_cuda $N
	let N=N+1
done

