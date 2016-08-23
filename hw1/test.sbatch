#!/bin/bash
#SBATCH --job-name=mpcs_test
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --time=00:05:00
#SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load intel
module load mkl

./test
