import argparse
import os
from Bio.PDB import PDBParser, PDBList

# Command-line arguments configuration
parser = argparse.ArgumentParser(description="Download and parse a protein structure directly from the PDB.")
parser.add_argument("-id", "--pdbid", required=True, help="4-letter PDB code (e.g., 1ubq or 1crn)")

args = parser.parse_args()
pdb_id = args.pdbid.lower()

# PDB download module configuration
pdbl = PDBList()
expected_file = f"pdb{pdb_id}.ent"

# Check for local file to avoid redundant downloads
if not os.path.exists(expected_file):
    print(f"Connecting to PDB servers to download: {pdb_id.upper()}")
    file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
else:
    print(f"Local file found for structure {pdb_id.upper()}. Skipping download.")
    file_path = expected_file

# Load PDB structure
pdb_parser = PDBParser(QUIET=True)
structure = pdb_parser.get_structure(pdb_id.upper(), file_path)

print("")
print(f"Structure {structure.get_id()} successfully loaded.")

# Calculate total atoms as a validation step
total_atoms = len(list(structure.get_atoms()))
print(f"Total atoms in the structure: {total_atoms}")