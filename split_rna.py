import sys

sys.path.append("..")
import pandas as pd
import os

# Append the 'rlclip/preprocessor' subdirectory to the system path
sys.path.append(os.path.join("..", "rlclip", "preprocessor"))

# Import classes from the 'rna3db' module in the 'rlclip.preprocessor' package
from rlclip.preprocessor.rna3db import StructureFile, ModificationHandler, mmCIFParser


def write_rna_rna3db(cif_path, rna_chain_id, pdbid, output_dir):
    modification_handler = ModificationHandler()
    complex_structure_rna3db = StructureFile(
        cif_path, modification_handler, include_atoms=True
    )

    complex_structure_rna3db.write_mmcif_chain(
        f"{output_dir}/{pdbid}_{rna_chain_id}.cif",
        author_id=rna_chain_id,
    )


def main(input_csv, output_dir, cif_base_path):
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read the input CSV file
    df = pd.read_csv(input_csv)

    # Process each row in the DataFrame
    for index, row in df.iterrows():
        pdb_id = row["PDB_ID"]
        rna_chain_id = row["RNA_Chain"]
        cif_path = os.path.join(cif_base_path, f"{pdb_id}.cif")
        write_rna_rna3db(cif_path, rna_chain_id, pdb_id, output_dir)


if __name__ == "__main__":
    import argparse

    # Setting up argument parser
    parser = argparse.ArgumentParser(
        description="Process RNA chains from a CSV file and write to mmCIF files."
    )
    parser.add_argument(
        "input_csv",
        type=str,
        help="Path to the input CSV file containing PDB_ID and RNA_Chain columns",
    )
    parser.add_argument(
        "output_dir", type=str, help="Directory where the output files will be saved"
    )
    parser.add_argument(
        "cif_base_path",
        type=str,
        help="Base path to the directory containing CIF files",
    )

    args = parser.parse_args()

    main(args.input_csv, args.output_dir, args.cif_base_path)

# python process_rna_chains.py input.csv ./output_dir ../hariboss_process/complex_database/complexes_mmcif


# # Create a list of tuples with PDB IDs and RNA chains
# data = [
#     ('4tua', 'XA'),
#     ('4tub', 'XA'),
#     ('4tud', 'XA'),
#     ('4tue', 'XA'),
#     ('5doy', '2A'),
#     ('6az3', '2'),
#     ('6v39', 'AN1'),
#     ('6ydp', 'AA'),
#     ('6ydw', 'AA'),
#     ('6yl5', 'K'),
#     ('6ymi', 'O'),
#     ('6ymj', 'C'),
#     ('6ymj', 'I'),
#     ('6ymj', 'O'),
#     ('7tql', '3'),
#     ('8ceu', 'a')
# ]

# # Create a DataFrame
# df = pd.DataFrame(data, columns=['PDB_ID', 'RNA_Chain'])
