import argparse
import sys
from Bio.PDB import PDBParser

# Command-line arguments configuration
parser = argparse.ArgumentParser(description="Calculate distances between all atoms of two specified residues.")
parser.add_argument("-i", "--input", required=True, help="Input PDB file")
parser.add_argument("-c", "--chain", default="A", help="Chain ID")
parser.add_argument("-r1", "--res1", required=True, type=int, help="Sequence number of the first residue")
parser.add_argument("-r2", "--res2", required=True, type=int, help="Sequence number of the second residue")

args = parser.parse_args()

# Load PDB structure
pdb_parser = PDBParser(QUIET=True)
structure = pdb_parser.get_structure("Protein", args.input)

try:
    # Access hierarchy: Structure -> Model 0 -> Chain -> Residue
    residue1 = structure[0][args.chain][args.res1]
    residue2 = structure[0][args.chain][args.res2]
    
    # Extract atoms and sort them by serial number to ensure ordered output
    atoms_r1 = sorted(list(residue1.get_atoms()), key=lambda a: a.get_serial_number())
    atoms_r2 = sorted(list(residue2.get_atoms()), key=lambda a: a.get_serial_number())
    
    print(f"Calculating distances between {residue1.get_resname()} {args.res1} and {residue2.get_resname()} {args.res2}:")
    print("")
    
    # Calculate distances between all atom pairs
    for a1 in atoms_r1:
        for a2 in atoms_r2:
            distance = a1 - a2
            print(f"Atom {a1.get_name():<4} <--> Atom {a2.get_name():<4} : {distance:.2f} A")
            
    print("")
    print(f"Total distances calculated: {len(atoms_r1) * len(atoms_r2)}")

except KeyError:
    print(f"Error: Residue {args.res1} or {args.res2} not found in chain {args.chain}.")
    sys.exit(1)