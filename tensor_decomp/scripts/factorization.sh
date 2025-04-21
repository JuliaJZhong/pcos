#!/bin/bash

# Job Flags
#SBATCH --job-name=factorization_%j
#SBATCH --output=outfiles/factorization_%j.out
#SBATCH --error=stderr/factorization_%j.err
#SBATCH --partition=sched_mit_hill
#SBATCH --time=2:00:00
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jjzhong@mit.edu

# Set up environment
# module load cuda/12.4.0
source /home/jjzhong/miniforge3/bin/activate tensorz

# Run your application
python src/factorization.py