from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser, is_aa
import conditions_biopython
from timeit import default_timer as timer
import os
from sys import argv
from psutil import Process


def main():
    """
    Main function for processing PDB or CIF files to detect and analyze atomic interactions.

    This function performs the following steps:
    1. Sets CPU affinity to a specific core if a core number is provided as an argument.
    2. Parses command-line arguments to determine the folder containing files to process and optional output settings.
    3. Scans the specified folder for files or uses the provided file directly if only one file is given.
    4. Processes each file to extract protein structures and compute atomic interactions.
    5. Filters interactions based on predefined conditions and outputs the results.
    6. Reports the total processing time.

    Command-line Arguments:
    - folder (str): The path to the folder containing PDB or CIF files. If a single file is specified, it is processed directly.
    - core_id (optional int): The CPU core number to which processing should be confined. If not provided, the function will run without setting CPU affinity.
    """
    
    if len(argv) == 3:
        core_id = argv[2]
        p = Process(os.getpid())
        p.cpu_affinity([int(core_id)])
        print(f"Running on core {core_id}")

    global_start = timer()

    folder = argv[1]
    files = [folder] if folder.endswith(".cif") else [entry.name for entry in os.scandir(folder)]
    files = sorted(files)

    for file in files:
        file_path = os.path.join(folder, file) if len(files) > 1 else file
        
        result = process_file(file_path)
        if result:
            protein_size, interactions_count, process_time = result
            print(file.split(".")[0], protein_size, interactions_count, f"{process_time:.4f}")
    
    print(f"Total time: {timer() - global_start}")


def process_file(file_path):
    start = timer()
    """
    Processes a single PDB or CIF file to extract and analyze atomic interactions.
    
    Parameters:
    - file_path: Path to the PDB or CIF file.

    Returns:
    - tuple: (protein_size, len(interactions), processing_time) or None if skipped.
    """
    try:
        parser = PDBParser(QUIET=True) if file_path.endswith(".pdb") else MMCIFParser(QUIET=True)
        structure = next(parser.get_structure('structure_id', file_path).get_models())
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return None

    protein_size = sum(1 for _ in structure.get_residues() if is_aa(_))
    if protein_size > 25000:
        print(f"Skipping ID '{file_path}'. Size: {protein_size} residues")
        return None

    # Filter out waters and heteroatoms
    atoms = [atom for atom in structure.get_atoms() if not atom.get_parent().get_resname() == 'HOH' and not atom.get_parent().id[0] == 'W']
    
    ns = NeighborSearch(atoms)
    radius = 6.0
    interactions = []

    for atom in atoms:
        residue_number = atom.get_parent().get_id()[1]
        chain1 = atom.get_parent().get_parent().get_id()
        atom1_name = f"{atom.get_parent().get_resname()}:{atom.get_name()}"
        
        if atom1_name not in conditions_biopython.contact_types:
            continue

        for neighbor in ns.search(atom.coord, radius):
            neighbor_number = neighbor.get_parent().get_id()[1]
            chain2 = neighbor.get_parent().get_parent().get_id()
            
            if atom == neighbor or residue_number >= neighbor_number:
                continue

            atom2_name = f"{neighbor.get_parent().get_resname()}:{neighbor.get_name()}"
            if atom2_name not in conditions_biopython.contact_types:
                continue

            distance = atom - neighbor
            
            if atom1_name in conditions_biopython.contact_types and atom2_name in conditions_biopython.contact_types:
                for interaction, condition in conditions_biopython.contact_conditions.items():
                    if interaction == 'hydrogen_bond' and (abs(neighbor_number - residue_number) <= 3):
                        continue
                    min_distance, max_distance = conditions_biopython.categories[interaction]
                    if min_distance <= distance <= max_distance and condition(atom1_name, atom2_name):
                        interactions.append((atom1_name, residue_number, atom2_name, neighbor_number, interaction, distance, chain1, chain2))
        
    interactions = sorted(interactions, key=lambda x:x[4])
    
    return protein_size, len(interactions), timer() - start


if __name__ == "__main__":
    main()