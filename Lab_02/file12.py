from Bio.PDB import PDBParser

def compute_distance(coord1, coord2):
    x1, y1, z1 = coord1
    x2, y2, z2 = coord2
    distance = ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
    return distance

def find_farthest_atom(structure):
    atoms = structure.get_atoms()
    max_distance = 0
    farthest_atom = None

    for atom1 in atoms:
        current_distance = 0
        for atom2 in atoms:
            current_distance += compute_distance(atom1.get_coord(), atom2.get_coord())

        if current_distance > max_distance:
            max_distance = current_distance
            farthest_atom = atom1

    return farthest_atom

if __name__ == "__main__":
    pdb_file = "4DFR.pdb"
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein_structure", pdb_file)

    farthest_atom = find_farthest_atom(structure)

    print(f"Information for the Farthest Atom:")
    print(f"Atom Name: {farthest_atom.get_name()}")
    print(f"Residue Name: {farthest_atom.get_parent().get_resname()}")
    print(f"Chain ID: {farthest_atom.get_full_id()[2]}")
    print(f"Coordinates: {farthest_atom.get_coord()}")

