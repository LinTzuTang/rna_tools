#!/bin/bash
#SBATCH --job-name=cif2pdb   # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=96gb                   # Job memory request
#SBATCH --time=96:00:00             # Time limit hrs:min:sec
#SBATCH --output=c2p.log   # Standard output and error log
pwd; hostname; date

module load conda
conda activate hariboss

python convert_cif2pdb.py

date