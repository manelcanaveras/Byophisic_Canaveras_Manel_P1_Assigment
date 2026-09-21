import argparse
from Bio.PDB import PDBParser

# Command-line arguments configuration
parser = argparse.ArgumentParser(description="Finds disulfide bridges (S-S) between CYS residues.")
parser.add_argument("-i", "--input", required=True, help="Input PDB file (e.g., 1ubq.pdb)")
parser.add_argument("-d", "--distance", type=float, default=3.0, help="Maximum distance in Angstroms (default: 3.0)")

args = parser.parse_args()

# Load PDB structure
pdb_parser = PDBParser(QUIET=True)
structure = pdb_parser.get_structure("Protein", args.input)

# Use the first model to prevent multi-model redundancy
model = structure[0]

# Collect SG atoms specifically from Cysteine residues
sulfur_atoms = []

for chain in model:
    for residue in chain:
        if residue.get_resname() == 'CYS' and 'SG' in residue:
            sulfur_atoms.append(residue['SG'])

disulfide_bridges = []

# Calculate distances to identify S-S bridges
for i in range(len(sulfur_atoms)):
    for j in range(i + 1, len(sulfur_atoms)):
        s1 = sulfur_atoms[i]
        s2 = sulfur_atoms[j]
        
        # Ensure atoms belong to different residues
        if s1.get_parent() != s2.get_parent():
            distance = s1 - s2
            
            if distance < args.distance:
                disulfide_bridges.append((s1, s2, distance))

# Sort by chain and residue sequence numbers to meet formatting requirements
disulfide_bridges.sort(key=lambda x: (
    x[0].get_parent().get_parent().get_id(),
    x[0].get_parent().get_id()[1],
    x[1].get_parent().get_parent().get_id(),
    x[1].get_parent().get_id()[1]
))

# Print formatted results
for s1, s2, dist in disulfide_bridges:
    res1 = s1.get_parent()
    res2 = s2.get_parent()
    
    print(f"S-S Bridge: {res1.get_resname()} {res1.get_id()[1]:<3} <--> "
          f"{res2.get_resname()} {res2.get_id()[1]:<3} : {dist:.2f} A")

print(f"Total disulfide bridges found: {len(disulfide_bridges)}")