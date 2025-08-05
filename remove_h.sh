#!/bin/bash
#SBATCH --job-name=RemoveH_rna   # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=96gb                   # Job memory request
#SBATCH --time=96:00:00             # Time limit hrs:min:sec
#SBATCH --output=RemoveH_rnas.log   # Standard output and error log
pwd; hostname; date

module load conda
conda activate pymol

python remove_h_pymol.py

date