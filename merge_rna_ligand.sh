#!/bin/bash
#SBATCH --job-name=merge_rna_ligand    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=64
#SBATCH --mem=96gb                   # Job memory request
#SBATCH --partition=bigmem
#SBATCH --qos=yanjun.li-b
#SBATCH --time=96:00:00             # Time limit hrs:min:sec
#SBATCH --output=merge_rna_ligand.log   # Standard output and error log
pwd; hostname; date

module load conda
conda activate hariboss

python merge_rna_ligand.py ../hariboss_process/hariboss_clean_20240621_error_revised_.csv ../hariboss_process/hariboss_20240621_error_revised_result/parsed_rna3db_rnas_pdb_add_H/ ../hariboss_process/hariboss_20240621_error_revised_result/parsed_ligands_pdb_add_H/ ../hariboss_process/hariboss_20240621_error_revised_result/merged_rna_ligand_add_H/
date