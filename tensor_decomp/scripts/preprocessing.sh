#!/bin/bash

# Job Flags
#SBATCH --job-name=preprocessing
#SBATCH --output=outfiles/preprocessing.out
#SBATCH --error=stderr/preprocessing.err
#SBATCH --partition=mit_normal
#SBATCH -t 0-00:30
#SBATCH --mem-per-cpu=32GB
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jjzhong@mit.edu

# Set up environment
source /home/jjzhong/miniforge3/bin/activate tensorz

nvidia-smi
# Run your application
python src/preprocessing.py