#!/bin/bash

#SBATCH --job-name=docking
#SBATCH -n 100
##SBATCH --partition=gpu
#SBATCH --nodes=5
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=20
##SBATCH --gres=gpu:1
#SBATCH --output=RNA-error.out
##SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00
##SBATCH --nodelist=c1007a-s11
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept-b
#SBATCH --mem=20G

# Load the conda environment
#source /path/to/conda.sh  # Adjust the path to where your conda.sh is located
#conda activate pymol_env  # Replace with the name of your conda environment

# Define directories
INPUT_DIR="/blue/mdisney/a.taghavi/HARIBOSS-DOCKING/convert-all-pdb-to-pdbqt/docking_benchmark-docking_file-rna_chain_remove"
OUTPUT_DIR="/blue/mdisney/a.taghavi/HARIBOSS-DOCKING/convert-all-pdb-to-pdbqt/docking_benchmark-docking_file-rna_chain_remove_pdbqt"
NUM_WORKERS=$SLURM_CPUS_PER_TASK

# Run the master Python script
python process_pdb_files.py $INPUT_DIR $OUTPUT_DIR $NUM_WORKERS

# Deactivate the conda environment
conda deactivate

