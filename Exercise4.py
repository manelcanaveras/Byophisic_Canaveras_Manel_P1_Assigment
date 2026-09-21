import argparse
import sys
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_1to3

# Command-line arguments configuration
parser = argparse.ArgumentParser(description="List CA coordinates for a specific amino acid type.")
parser.add_argument("-i", "--input", required=True, help="Input PDB file")
parser.add_argument("-t", "--type", required=True, help="Residue type")

args = parser.parse_args()

# Convert 1-letter amino acid code to 3-letter code if necessary
residue_type = args.type.upper()

if len(residue_type) == 1:
    try:
        residue_type = protein_letters_1to3[residue_type].upper()
    except KeyError:
        print(f"Error: '{residue_type}' is not a standard amino acid letter code.")
        sys.exit(1)

# Load PDB structure
pdb_parser = PDBParser(QUIET=True)
structure = pdb_parser.get_structure("Protein", args.input)

# Use only the first model to prevent duplicate entries
model = structure[0]

# Collect matching residues
found_residues = []

for chain in model:
    for residue in chain:
        if residue.get_resname() == residue_type and 'CA' in residue:
            found_residues.append((chain.get_id(), residue, residue['CA']))

# Sort by chain ID and residue sequence number
found_residues.sort(key=lambda x: (x[0], x[1].get_id()[1]))

print(f"Searching CA atoms for residue type: {residue_type}")
print("")

# Formatted output
for chain_id, residue, atom_ca in found_residues:
    res_num = residue.get_id()[1]
    coords = atom_ca.get_coord()
    
    print(f"Chain {chain_id} | Residue: {residue_type} {res_num:<3} | CA Coordinates: {coords}")
print("")
print(f"Total {residue_type} residues found: {len(found_residues)}")