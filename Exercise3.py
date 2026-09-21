import argparse
from Bio.PDB import PDBParser

# Command-line arguments configuration
parser = argparse.ArgumentParser(description="Find potential hydrogen bonds (O, N, S atoms) within a given distance.")
parser.add_argument("-i", "--input", required=True, help="Input PDB file (e.g., 1ubq.pdb)")
parser.add_argument("-d", "--distance", type=float, default=3.5, help="Maximum distance in Angstroms (default: 3.5)")

args = parser.parse_args()

# Load PDB structure
pdb_parser = PDBParser(QUIET=True)
structure = pdb_parser.get_structure("Protein", args.input)

# Extract first model to avoid multi-model redundancy
model = structure[0]

# Filter polar atoms (O, N, S)
polar_elements = ('O', 'N', 'S')
polar_atoms = [atom for atom in model.get_atoms() if atom.get_name()[0] in polar_elements]

h_bonds = []

# Calculate distances between polar atoms
for i in range(len(polar_atoms)):
    for j in range(i + 1, len(polar_atoms)):
        atom1 = polar_atoms[i]
        atom2 = polar_atoms[j]
        
        # Skip atoms belonging to the same residue
        if atom1.get_parent() != atom2.get_parent():
            distance = atom1 - atom2
            
            if distance < args.distance:
                h_bonds.append((atom1, atom2, distance))

# Sort results by residue sequence number, then by atom serial number
h_bonds.sort(key=lambda x: (
    x[0].get_parent().get_id()[1], 
    x[0].get_serial_number(),
    x[1].get_parent().get_id()[1],
    x[1].get_serial_number()
))

# Formatted output
for a1, a2, dist in h_bonds:
    res1 = a1.get_parent()
    res2 = a2.get_parent()
    
    print(f"{res1.get_resname()} {res1.get_id()[1]:<3} ({a1.get_name():<3}) -- "
          f"{res2.get_resname()} {res2.get_id()[1]:<3} ({a2.get_name():<3}) : {dist:.2f} A")

print(f"\nTotal potential hydrogen bonds: {len(h_bonds)}")