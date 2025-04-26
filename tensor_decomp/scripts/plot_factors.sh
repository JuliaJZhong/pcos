#!/bin/bash

# Job Flags
#SBATCH --job-name=factors-04-22_%j
#SBATCH --output=outfiles/factors-04-22_%j.out
#SBATCH --error=stderr/factors-04-22_%j.err
#SBATCH --partition=mit_normal
#SBATCH -t 0-00:30
#SBATCH --mem-per-cpu=32GB
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jjzhong@mit.edu

# Set up environment
source /home/jjzhong/miniforge3/bin/activate tensorz

# Run your application
python src/plot_factors.py