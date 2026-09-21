import argparse
from Bio.PDB import PDBParser

# Command-line arguments configuration
parser = argparse.ArgumentParser(description="Maps backbone peptide bonds (C-N) based on distance threshold.")
parser.add_argument("-i", "--input", required=True, help="Input PDB file (e.g., 1ubq.pdb)")
parser.add_argument("-d", "--distance", type=float, default=2.5, help="Distance cutoff in Angstroms (default: 2.5)")

args = parser.parse_args()

# Load PDB structure
pdb_parser = PDBParser(QUIET=True)
structure = pdb_parser.get_structure("Protein", args.input)

# Use the first model to avoid NMR duplicate artifacts
model = structure[0]

# Separate backbone Carbon and Nitrogen atoms
c_atoms = []
n_atoms = []

for chain in model:
    for residue in chain:
        if 'C' in residue:
            c_atoms.append(residue['C'])
        if 'N' in residue:
            n_atoms.append(residue['N'])

peptide_bonds = []

# Measure distances to identify bonded atoms
for c in c_atoms:
    for n in n_atoms:
        # A peptide bond connects two different residues
        if c.get_parent() != n.get_parent():
            distance = c - n
            
            if distance < args.distance:
                peptide_bonds.append((c, n, distance))

# Sort by chain ID and the sequence number of the Carbon residue
peptide_bonds.sort(key=lambda x: (
    x[0].get_parent().get_parent().get_id(), 
    x[0].get_parent().get_id()[1]
))

# Print results
for c, n, dist in peptide_bonds:
    res_c = c.get_parent()
    res_n = n.get_parent()
    
    print(f"{res_c.get_resname()} {res_c.get_id()[1]:<3} (C) --> "
          f"{res_n.get_resname()} {res_n.get_id()[1]:<3} (N) : {dist:.2f} A")

print(f"\nTotal peptide bonds found: {len(peptide_bonds)}")