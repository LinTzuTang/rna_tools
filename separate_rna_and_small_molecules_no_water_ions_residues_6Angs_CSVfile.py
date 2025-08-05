from Bio.PDB import MMCIFParser, MMCIFIO, NeighborSearch, Selection, Structure, Model, Chain
import pandas as pd
import os

class RNAandSmallMoleculeSelect:
    def __init__(self, small_molecule_ids):
        self.rna_residues = {'A', 'U', 'G', 'C'}  # RNA residues
        self.small_molecule_ids = small_molecule_ids  # List of small molecule IDs

    def accept_residue(self, residue):
        # Check for RNA residues
        if residue.resname.strip() in self.rna_residues:
            return False
        # Check for small molecules
        is_protein = residue.id[0] == " "  # Protein residues have a blank hetflag
        is_small_molecule = not is_protein and residue.resname.strip() in self.small_molecule_ids
        return is_small_molecule

def extract_first_model(structure):
    """Extract the first model if there are multiple models in the structure."""
    for model in structure:
        return model

def select_nearby_residues(structure, small_molecule_atoms, cutoff=6.0):
    """Select residues within the cutoff distance from small molecule atoms."""
    ns = NeighborSearch(Selection.unfold_entities(structure, 'A'))  # Search across all atoms
    nearby_residues = set()

    for atom in small_molecule_atoms:
        neighbors = ns.search(atom.coord, cutoff)
        print(f"Atom {atom} has {len(neighbors)} neighbors within {cutoff} Å")  # Debug print
        for neighbor in neighbors:
            residue = neighbor.get_parent()  # Get the residue of the nearby atom
            nearby_residues.add(residue)  # Add residue to the set

    return nearby_residues

def read_small_molecule_ids(csv_file):
    """Read small molecule IDs from a CSV file."""
    df = pd.read_csv(csv_file)
    return set(df['id'].str.strip())  # Return a set of IDs

def process_mmCIF(input_file, output_file, csv_file, cutoff=6.0):
    print(f"Processing file: {input_file} with cutoff: {cutoff}")  # Debug print
    small_molecule_ids = read_small_molecule_ids(csv_file)
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("", input_file)

    model = extract_first_model(structure)
    
    selected_residues = set()
    small_molecule_atoms = []

    # Initialize the selector with small molecule IDs
    selector = RNAandSmallMoleculeSelect(small_molecule_ids)

    # Step 1: Select small molecules based on CSV file
    for chain in model:
        for residue in chain:
            if selector.accept_residue(residue):
                if residue.id[0] != " ":  # Small molecules
                    small_molecule_atoms.extend(residue.get_atoms())
                    selected_residues.add(residue)  # Include small molecules

    print(f"Selected {len(selected_residues)} residues and {len(small_molecule_atoms)} small molecule atoms.")  # Debug print

    # Step 2: Select residues within the cutoff distance from small molecules
    if small_molecule_atoms:
        nearby_residues = select_nearby_residues(model, small_molecule_atoms, cutoff)
        print(f"Found {len(nearby_residues)} residues within cutoff distance.")  # Debug print
        selected_residues.update(nearby_residues)

    # Create a new structure to store the selected residues
    new_structure = Structure.Structure("Selected_Structure")
    new_model = Model.Model(0)
    new_structure.add(new_model)

    for residue in selected_residues:
        chain_id = residue.parent.id
        if chain_id not in new_model.child_dict:
            new_chain = Chain.Chain(chain_id)
            new_model.add(new_chain)
        new_model[chain_id].add(residue.copy())  # Safely add the residue

    # Save to output file
    io = MMCIFIO()
    io.set_structure(new_structure)
    io.save(output_file)
    print(f"Saved to {output_file}")

def main(input_dir, output_dir, csv_file, cutoff=6.0):
    os.makedirs(output_dir, exist_ok=True)
    
    for filename in os.listdir(input_dir):
        if filename.endswith(".cif"):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, filename.replace(".cif", "_rna_small_nearby.cif"))
            print(f"Processing {filename}")
            process_mmCIF(input_file, output_file, csv_file, cutoff)
            print(f"Finished processing {filename}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Process mmCIF files to extract small molecules from a CSV list and residues within a given distance from these small molecules.")
    parser.add_argument("input_dir", help="Directory containing mmCIF files.")
    parser.add_argument("output_dir", help="Directory to save output files.")
    parser.add_argument("csv_file", help="CSV file containing small molecule IDs.")
    parser.add_argument("--cutoff", type=float, default=6.0, help="Distance cutoff in angstroms for selecting nearby residues.")
    
    args = parser.parse_args()
    main(args.input_dir, args.output_dir, args.csv_file, args.cutoff)

