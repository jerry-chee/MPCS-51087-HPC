#!/bin/bash
#SBATCH --job-name=mpcs51087_test
#SBATCH --account=mpcs51087
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --time=00:1:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

module load cuda
./test 10 10
