import argparse
from Bio.PDB import PDBParser

# Command-line argument configuration
cmd_parser = argparse.ArgumentParser(description="List atoms of a specific residue.")
cmd_parser.add_argument("-i", "--input", required=True, help="Input PDB file (e.g., 1ubq.pdb)")
cmd_parser.add_argument("-c", "--chain", default="A", help="Chain ID (default: A)")
cmd_parser.add_argument("-r", "--resnum", required=True, type=int, help="Residue sequence number")

args = cmd_parser.parse_args()

# PDB structure initialization
parser = PDBParser(QUIET=True)
structure = parser.get_structure("Protein", args.input)

try:
    # Access hierarchy: Structure -> Model 0 -> Chain -> Residue
    residue = structure[0][args.chain][args.resnum]
    
    # Extract atoms and sort them by atom serial number as required
    atoms = sorted(residue.get_atoms(), key=lambda a: a.get_serial_number())
    
    print(f"Atoms for {residue.get_resname()} {args.resnum} (Chain {args.chain}):")
    
    for atom in atoms:
        atom_name = atom.get_name()
        coords = atom.get_coord()
        
        # Formatted output for readability
        print(f"Atom: {atom_name:<4} | Coordinates (x, y, z): {coords}")

except KeyError:
    print(f"Error: Residue {args.resnum} not found in chain {args.chain}.")