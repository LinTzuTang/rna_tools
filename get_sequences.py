from Bio.PDB.MMCIFParser import MMCIFParser
import os

def extract_sequences_from_cif(cif_file):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('structure', cif_file)
    sequences = {}

    # Extract PDB ID from the CIF file metadata
    cif_dict = parser._mmcif_dict
    pdb_id = cif_dict["_entry.id"][0]
    
    for model in structure:
        for chain in model:
            sequence = []
            for residue in chain:
                if residue.id[0] == ' ':  # Skip hetero residues
                    residue_id = residue.resname
                    sequence.append(residue_id)
            sequences[f'{pdb_id}'] = ''.join(sequence)
    
    return sequences

def process_all_cif_files(directory):
    all_sequences = {}
    for filename in os.listdir(directory):
        if filename.endswith(".cif"):
            cif_file = os.path.join(directory, filename)
            sequences = extract_sequences_from_cif(cif_file)
            all_sequences.update(sequences)
    return all_sequences

def save_to_fasta(sequences, fasta_file):
    with open(fasta_file, 'w') as f:
        for chain, sequence in sequences.items():
            f.write(f">{chain}\n")
            f.write(f"{sequence}\n")

def main():
    # # Replace 'your_directory_path' with the path to your directory containing CIF files
    directory = '../hariboss_process/hariboss_20240621_error_revised_result/parsed_rna3db_rnas'
    # # Replace 'output.fasta' with the desired path for the output FASTA file
    fasta_file = 'hariboss_20240621_error_revised_rna.fasta'
    
    all_sequences = process_all_cif_files(directory)
    save_to_fasta(all_sequences, fasta_file)
    print(f"Sequences have been saved to {fasta_file}")

if __name__ == "__main__":
    main()
