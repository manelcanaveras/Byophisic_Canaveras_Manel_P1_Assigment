import argparse
from Bio.PDB import PDBParser

def find_ca_contacts(pdb_path, max_distance):
    """
    Parses a PDB file and finds all CA-CA contacts within a specified distance.
    Returns a list of tuples containing the two residues and their distance, 
    sorted by residue numbers.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    
    # Extract all CA atoms from the structure
    ca_atoms = [atom for atom in structure.get_atoms() if atom.get_name() == 'CA']
    
    contacts = []
    
    # Calculate distances between all unique CA pairs
    for i in range(len(ca_atoms)):
        for j in range(i + 1, len(ca_atoms)):
            distance = ca_atoms[i] - ca_atoms[j]
            
            if distance < max_distance:
                res1 = ca_atoms[i].get_parent()
                res2 = ca_atoms[j].get_parent()
                contacts.append((res1, res2, distance))
                
    # Ensure the output is sorted by the first residue number, then the second
    contacts.sort(key=lambda x: (x[0].get_id()[1], x[1].get_id()[1]))
    
    return contacts

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find CA atom pairs within a specified distance in a PDB file.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input PDB file (e.g., 1ubq.pdb)")
    parser.add_argument("-d", "--distance", type=float, default=5.0, help="Maximum distance threshold in Angstroms")
    
    args = parser.parse_args()
    
    contacts = find_ca_contacts(args.input, args.distance)
    
    # Format and print the sorted results
    for res1, res2, dist in contacts:
        res1_name = res1.get_resname()
        res1_id = res1.get_id()[1]
        res2_name = res2.get_resname()
        res2_id = res2.get_id()[1]
        
        print(f"{res1_name} {res1_id} -- {res2_name} {res2_id} : {dist:.2f} A")
        
    print(f"\nTotal contacts found: {len(contacts)}")