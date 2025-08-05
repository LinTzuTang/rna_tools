
# 🧬 RNA Structure & Ligand Analysis Toolkit

This toolkit provides a collection of automated Python scripts for RNA structure processing, sequence/ligand feature extraction, and batch analysis.

---

## ⚙️ Requirements

Install the core dependencies with:

```bash
pip install -r requirements.txt
```


Note:

For 3D visualization and file conversion, you may also need to install `PyMOL` and `OpenBabel` as system tools (conda install -c conda-forge openbabel pymol-open-source).

For parallel processing, `joblib` is used.
`RDKit` is best installed via conda if possible.

## 📁 Tool Overview

| File Name | Brief Description |
| :--- | :--- |
| `remove_h_pymol.py` | Remove hydrogen atoms from PDB files using PyMOL. |
| `merge_rna_ligand.py` | Merge separate RNA and ligand PDB files into a single structure. |
| `generate_contact_map.py` | Calculate RNA–ligand contact maps based on atomic distances. |
| `split_rna.py` | Split PDB/CIF structures into RNA, ligand, and protein components. |
| `convert_cif2pdb.py` | Convert mmCIF to PDB format, keeping chain info. |
| `separate_rna_and_small_molecules_no_water_ions_residues_6Angs_CSVfile.py` | Extract RNA and small molecules (exclude water/ions) within 6Å. Batch processing via CSV. |
| `AddH_pymol_save_pdbqt_obabel_parallel_HPC.py` | Add hydrogens & convert to PDBQT using PyMOL + OpenBabel, parallel HPC support. |
| `get_sequences.py` | Extract RNA sequences (FASTA) or ligand SMILES from structures. |
| `AddH_pymol.py` | Add hydrogens to PDB structure via PyMOL. |
| `get_ligands_smiles_fingerprint.py` | Extract ligand SMILES & generate RDKit fingerprints. |


---

## ⛓️ Pipeline Overview

These tools can be combined for a full RNA-ligand structural workflow:
1. **Split complex structures:** `split_rna.py`
2. **Standardize structures:** `remove_h_pymol.py` → `AddH_pymol.py`
3. **Merge components:** `merge_rna_ligand.py`
4. **Generate features:** `generate_contact_map.py`, `get_sequences.py`, `get_ligands_smiles_fingerprint.py`
5. **Batch/high-throughput:** Use parallel scripts for scaling (`AddH_pymol_save_pdbqt_obabel_parallel_HPC.py`, CSV-driven tools).

