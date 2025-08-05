import os
import argparse
import pandas as pd
from Bio import PDB

def merge_pdb_files(rna_file, ligand_file, output_file):
    """
    Merges an RNA PDB file and a ligand PDB file into a single PDB file.

    :param rna_file: Path to the RNA PDB file
    :param ligand_file: Path to the ligand PDB file
    :param output_file: Path where the merged PDB file will be saved
    """
    parser = PDB.PDBParser(QUIET=True)
    io = PDB.PDBIO()

    try:
        rna_structure = parser.get_structure("RNA", rna_file)
        ligand_structure = parser.get_structure("LIGAND", ligand_file)
    except Exception as e:
        print(f"Error parsing PDB files: {e}")
        return

    # Collect existing chain IDs to avoid conflicts
    existing_chain_ids = {chain.id for model in rna_structure for chain in model}

    # Renaming ligand chains to avoid conflicts
    for model in ligand_structure:
        for chain in model:
            new_chain_id = chain.id
            while new_chain_id in existing_chain_ids:
                new_chain_id = (
                    chr(ord(new_chain_id) + 1) if new_chain_id.isalpha() else "A"
                )
            chain.id = new_chain_id
            existing_chain_ids.add(new_chain_id)
            rna_structure[0].add(chain.copy())

    try:
        io.set_structure(rna_structure)
        io.save(output_file)
        print(f'Successfully merged: {rna_file} and {ligand_file} into {output_file}')
    except Exception as e:
        print(f"Error saving merged PDB file: {e}")

def get_file_starting_with_base_id(directory, base_id):
    """
    Retrieves a file from the directory that starts with the given base_id.

    :param directory: Directory to search for the file
    :param base_id: Base ID that the file name should start with
    :return: File name that starts with the base ID or None if not found
    """
    for file_name in os.listdir(directory):
        if file_name.startswith(base_id):
            return file_name
    return None

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Merge RNA and ligand PDB files.")
    parser.add_argument("input_csv", help="Path to the input index CSV file.")
    parser.add_argument("rna_dir", help="Directory containing RNA PDB files.")
    parser.add_argument("ligand_dir", help="Directory containing ligand PDB files.")
    parser.add_argument("output_dir", help="Directory to save the merged PDB files.")

    # Parse arguments
    args = parser.parse_args()

    # Load the input CSV file
    rrna_shorter_than_621 = pd.read_csv(args.input_csv, keep_default_na=False)

    # Create the output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Iterate through the DataFrame and merge the corresponding RNA and ligand files
    for index, row in rrna_shorter_than_621.iterrows():
        rna_base_id = row["pdbid"] + "_" + row["rna_chain"]
        rna_file = get_file_starting_with_base_id(args.rna_dir, rna_base_id)
        rna_path = os.path.join(args.rna_dir, rna_file)
        
        ligand_base_id = row["pdbid"] + "_" + row["ligand_chain"] + "_" + row["ligand_id"]
        ligand_file = get_file_starting_with_base_id(args.ligand_dir, ligand_base_id)
        ligand_path = os.path.join(args.ligand_dir, ligand_file) 
        
        output_filename = f"{row['pdbid']}_{row['rna_chain']}_{row['ligand_chain']}_{row['ligand_id']}.pdb"
        output_path = os.path.join(args.output_dir, output_filename)
        
        print(f"Merging files: {rna_path} and {ligand_path}")
        merge_pdb_files(rna_path, ligand_path, output_path)
        print(f'Merged: {rna_path} and {ligand_path} into {output_path}')

if __name__ == "__main__":
    main()
