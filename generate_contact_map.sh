#!/bin/bash
#SBATCH --job-name=generate_contact_map_removed_h    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=64
#SBATCH --mem=96gb                   # Job memory request
#SBATCH --partition=bigmem
#SBATCH --time=96:00:00             # Time limit hrs:min:sec
#SBATCH --output=generate_contact_map_removed_h.log   # Standard output and error log
pwd; hostname; date

module load conda
conda activate hariboss

python generate_contact_map.py

date
