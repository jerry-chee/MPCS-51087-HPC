#!/bin/bash
#SBATCH --job-name=mpcs51087_hw1
#SBATCH --output=hw1.out
#SBATCH --error=hw1.err
#SBATCH --time=00:10:00
#SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load intel
module load mkl

./hw1_p1
./hw1_p2
