#!/bin/bash
#SBATCH --job-name=AddH   # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=96gb                   # Job memory request
#SBATCH --time=96:00:00             # Time limit hrs:min:sec
#SBATCH --output=AddH_ligands.log   # Standard output and error log
pwd; hostname; date

module load conda
conda activate pymol

python AddH_pymol.py

date