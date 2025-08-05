import pandas as pd
import os
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import numpy as np
import pickle
import json
import warnings
import logging

# Configure logging
logging.basicConfig(filename='contact_map_processing_removed_h.log', level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s:%(message)s')

# Suppress PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

# Define constants and paths
MAX_BLOCKS = 50  # Define the max blocks constant
result_path = "../hariboss_process/hariboss_20240621_error_revised_result/contact_map_results"
dirs = {
    'ligand_row_dir': "../hariboss_process/hariboss_20240621_error_revised_result/parsed_ligands_pdb/",
    'ligand_with_h_dir': "../hariboss_process/hariboss_20240621_error_revised_result/parsed_ligands_pdb_add_H/",
    'ligand_removed_h_dir': "../hariboss_process/hariboss_20240621_error_revised_result/parsed_ligands_pdb_remove_H/",
    'rna_row_dir': "../hariboss_process/hariboss_20240621_error_revised_result/parsed_rna3db_rnas_pdb/",
    'rna_with_h_dir': "../hariboss_process/hariboss_20240621_error_revised_result/parsed_rna3db_rnas_pdb_add_H/",
    'rna_removed_h_dir': "../hariboss_process/hariboss_20240621_error_revised_result/parsed_rna3db_rnas_pdb_remove_H/"
}

# Read the CSV file into a DataFrame
df = pd.read_csv("../hariboss_process/hariboss_clean_20240621_error_revised.csv", keep_default_na=False)

# Check if DataFrame is loaded correctly
logging.info(f"DataFrame loaded with {len(df)} rows")

def generate_contact_map(distance_matrix, threshold=10):
    """
    Generate a contact map given a distance matrix and a threshold.
    """
    return (distance_matrix < threshold).astype(int)

def calculate_block_size(rna_length, max_blocks=MAX_BLOCKS):
    """
    Calculate the block size for RNA given its length and the maximum number of blocks.
    """
    if rna_length <= max_blocks:
        return 1
    return (rna_length + max_blocks - 1) // max_blocks

def get_file_paths(row):
    """
    Get the file paths for the ligand and RNA files from the dataframe row.
    """
    ligand_row_file = os.path.join(dirs['ligand_row_dir'], f"{row['pdbid']}_{row['ligand_chain']}_{row['ligand_id']}.pdb")
    ligand_with_h_file = os.path.join(dirs['ligand_with_h_dir'], f"{row['pdbid']}_{row['ligand_chain']}_{row['ligand_id']}_with_hydrogens.pdb")
    ligand_removed_h_file = os.path.join(dirs['ligand_removed_h_dir'], f"{row['pdbid']}_{row['ligand_chain']}_{row['ligand_id']}_without_hydrogens.pdb")
    rna_row_file = os.path.join(dirs['rna_row_dir'], f"{row['pdbid']}_{row['rna_chain']}.pdb")
    rna_with_h_file = os.path.join(dirs['rna_with_h_dir'], f"{row['pdbid']}_{row['rna_chain']}_with_hydrogens.pdb")
    rna_removed_h_file = os.path.join(dirs['rna_removed_h_dir'], f"{row['pdbid']}_{row['rna_chain']}_without_hydrogens.pdb")
    return pd.Series([ligand_row_file, ligand_with_h_file, ligand_removed_h_file, rna_row_file, rna_with_h_file, rna_removed_h_file])

def main():
    os.makedirs(result_path, exist_ok=True)

    # Apply the function to get file paths
    df[["ligand_pdb_row", "ligand_pdb_with_h", "ligand_pdb_removed_h", "rna_pdb_row", "rna_pdb_with_h", "rna_pdb_removed_h"]] = df.apply(get_file_paths, axis=1)

    # PDB parser
    pdb_parser = PDBParser()

    # Initialize dictionaries to store results
    # contact_map_dict_with_h = {}
    contact_map_dict_removed_h = {}
    # contact_map_dict_row = {}

    # Define the contact threshold
    contact_threshold = 10

    # Loop through each row in the DataFrame and process the data
    for index, row in df.iterrows():
        logging.info(f"Processing row {index + 1}/{len(df)}: {row['pdbid']}")
        use_h = "removed_h"
        # for use_h in ["with_h", "removed_h", "row"]:
        try:
            ligand_file = row[f"ligand_pdb_{use_h}"]
            rna_file = row[f"rna_pdb_{use_h}"]

            # Check if files exist before processing
            if not os.path.exists(ligand_file):
                logging.warning(f"File not found: {ligand_file}")
                continue
            if not os.path.exists(rna_file):
                logging.warning(f"File not found: {rna_file}")
                continue

            # Parse the structures
            ligand_structure = pdb_parser.get_structure("ligand", ligand_file)
            ligand = list(ligand_structure.get_atoms())[0].parent

            rna_structure = pdb_parser.get_structure("rna", rna_file)
            rna_chain = rna_structure[0][row["rna_chain"][0]]  # Ensure chain name fits PDB format

            # Calculate block size
            block_size = calculate_block_size(len(rna_chain))

            # Initialize matrices
            num_rna_blocks = (len(rna_chain) + block_size - 1) // block_size
            num_ligand_atoms = len(ligand)

            block_distance_matrix = np.full((num_rna_blocks, num_ligand_atoms), np.inf)
            original_distance_matrix = np.full((len(rna_chain), num_ligand_atoms), np.inf)

            rna_residues = list(rna_chain)

            # Calculate distances
            for i, rna_residue in enumerate(rna_residues):
                for j, ligand_atom in enumerate(ligand):
                    min_distance_in_block = np.inf
                    atoms = rna_residue.get_atoms()
                    for rna_atom in atoms:
                        dist = np.linalg.norm(rna_atom.coord - ligand_atom.coord)
                        original_distance_matrix[i, j] = min(original_distance_matrix[i, j], dist)
                        min_distance_in_block = min(min_distance_in_block, dist)
                    block_idx = i // block_size
                    block_distance_matrix[block_idx, j] = min(block_distance_matrix[block_idx, j], min_distance_in_block)

            # Generate contact maps
            min_distance = np.min(original_distance_matrix)
            contact_map_org = generate_contact_map(original_distance_matrix, contact_threshold)

            contact_map = np.where(np.sum(contact_map_org, axis=1) == 0, 0, 1)
            one_index = np.where(contact_map == 1)[0]

            # Create contact map dictionary
            contact_map_entry = {
                "pdbid": row["pdbid"],
                "ligand_id": row["ligand_id"],
                "ligand_resnum": row["ligand_resnum"],
                "rna_chain": row["rna_chain"],
                "ligand_chain": row["ligand_chain"],
                "min_contact_distance": min_distance,
                "contact_percentage": np.sum(contact_map == 1) / len(contact_map),
                "contact_map_2d": contact_map_org.tolist(),
                "contact_map_1d": contact_map.tolist(),
                "contact_residue_index": one_index.tolist(),
            }

            contact_map_dict_removed_h[row["pdbid"]] = contact_map_entry

            logging.info(f"Processed: {row['pdbid']} with hydrogens: {use_h}")

        except Exception as e:
            logging.error(f"Failed to process {row['pdbid']} with hydrogens: {use_h}: {e}")

    # Save the contact map dictionary
    with open(f"{result_path}/contact_map_dict_removed_h.pkl", "wb") as f:
        pickle.dump(contact_map_dict_removed_h, f)
    with open(f"{result_path}/contact_map_dict_removed_h.json", "w") as f:
        json.dump(contact_map_dict_removed_h, f)

    # Uncomment and save the contact map dictionaries for row and with_h if needed
    # with open(f"{result_path}/contact_map_dict_with_h.pkl", "wb") as f:
    #     pickle.dump(contact_map_dict_with_h, f)
    # with open(f"{result_path}/contact_map_dict_with_h.json", "w") as f:
    #     json.dump(contact_map_dict_with_h, f)

    # with open(f"{result_path}/contact_map_dict_row.pkl", "wb") as f:
    #     pickle.dump(contact_map_dict_row, f)
    # with open(f"{result_path}/contact_map_dict_row.json", "w") as f:
    #     json.dump(contact_map_dict_row, f)

if __name__ == "__main__":
    main()
