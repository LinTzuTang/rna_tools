import pandas as pd
import logging
import os
import pickle
import json
from rdkit import Chem
import openbabel.pybel as pybel
import requests

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define file paths
ligands_dir = "../hariboss_process/hariboss_20240621_error_revised_result/parsed_ligands_pdb/"
index_path = "../hariboss_process/hariboss_clean_20240621_error_revised_rrna_removed_non_rnasm_removed.csv"
result_path = "../hariboss_process/hariboss_20240621_error_revised_result_rrna_removed/ligands_smiles_fingerprint"


def extract_smiles_from_pdb(pdb_path):
    """
    Extract SMILES from a PDB file using RDKit.
    """
    try:
        mol = Chem.MolFromPDBFile(pdb_path)
        if mol:
            smiles = Chem.MolToSmiles(mol)
            return smiles
        else:
            logging.error(f"RDKit could not parse the PDB file: {pdb_path}")
            return None
    except Exception as e:
        logging.error(f"Error extracting SMILES from PDB file {pdb_path}: {e}")
        return None

def smiles2fp2(smiles):
    """
    Convert a SMILES string to a 2D fingerprint array.
    """
    mol = pybel.readstring("smi", smiles)
    fp2_bits = mol.calcfp(fptype="FP2").bits
    fp2_bit_one_hot = [0] * 1024
    for bit in fp2_bits:
        fp2_bit_one_hot[bit] = 1
    return fp2_bit_one_hot

def standardize_smiles(smiles_string):
    """
    Standardize the SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles_string)
    return Chem.MolToSmiles(mol)

def download_ligand_smiles(ligand_id):
    """
    Download ligand SMILES from PDB Bank.
    """
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data["rcsb_chem_comp_descriptor"]["smilesstereo"]
    return None

def process_ligands(df):
    contact_map_entries = []

    for index, row in df.iterrows():
        pdbid = row["pdbid"]
        ligand_chain_id = row["ligand_chain"]
        ligand_id = row["ligand_id"]
        
        ligand_path = f"{ligands_dir}/{pdbid}_{ligand_chain_id}_{ligand_id}.pdb"
        
        if os.path.exists(ligand_path):
            try:
                ligand_smiles = extract_smiles_from_pdb(ligand_path)
                canonical_ligand_smiles = standardize_smiles(ligand_smiles)
            except Exception as e:
                logging.error(f"Error processing PDB file {ligand_path}: {e}")
                continue

            try:
                download_ligand_smiles_str = download_ligand_smiles(ligand_id)
                canonical_download_ligand_smiles = standardize_smiles(download_ligand_smiles_str)
            except Exception as e:
                logging.error(f"Error downloading SMILES for ligand {ligand_id}: {e}")
                canonical_download_ligand_smiles = None
            
            parsed_fp2 = smiles2fp2(canonical_ligand_smiles)
            download_fp2 = smiles2fp2(canonical_download_ligand_smiles) if canonical_download_ligand_smiles else None

            contact_map_entry = {
                "pdbid": pdbid,
                "ligand_id": ligand_id,
                "ligand_resnum": row["ligand_resnum"],
                "rna_chain": row["rna_chain"],
                "ligand_chain": ligand_chain_id,
                "parsed_ligand_smiles": canonical_ligand_smiles if canonical_ligand_smiles else ligand_smiles,
                "downloaded_ligand_smiles": canonical_download_ligand_smiles if canonical_download_ligand_smiles else download_ligand_smiles_str,
                "parsed_fp2": parsed_fp2,
                "download_fp2": download_fp2,
            }
            
            contact_map_entries[f"{pdbid}"] = contact_map_entry
        else:
            logging.warning(f"Ligand file {ligand_path} does not exist.")
    
    return contact_map_entries

def save_results(contact_map_entries):
    os.makedirs(result_path, exist_ok=True)
    with open(f"{result_path}/ligands_smiles_fingerprint.pkl", "wb") as f:
        pickle.dump(contact_map_entries, f)
    with open(f"{result_path}/ligands_smiles_fingerprint.json", "w") as f:
        json.dump(contact_map_entries, f)
    logging.info("Processing complete. Results saved.")

def main():
    # Read the CSV file into a DataFrame
    df = pd.read_csv(index_path, keep_default_na=False)
    
    # Process ligands and generate contact map entries
    contact_map_entries = process_ligands(df)
    
    # Save the results
    save_results(contact_map_entries)

if __name__ == "__main__":
    main()
