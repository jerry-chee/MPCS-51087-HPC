#!/bin/bash
#SBATCH --job-name=mpcs51087_hw4
#SBATCH --output=hw4.out
#SBATCH --error=hw4.err
#SBATCH --time=00:5:00
#SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --open-mode=append|truncate

export OMP_NUM_THREADS=8
./transport
