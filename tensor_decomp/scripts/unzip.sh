#!/bin/bash

# Job Flags
#SBATCH --job-name=unzip
#SBATCH --output=outfiles/unzip.out
#SBATCH --partition=mit_normal
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jjzhong@mit.edu

# Set up environment
source /home/jjzhong/miniforge3/bin/activate tensorz

# Run your application
python src/unzip.py