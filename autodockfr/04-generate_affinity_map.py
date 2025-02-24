## 2. Generate the Affinity Maps and Gridbox Parameters
### Generate affinity maps using AGFR

import os
import subprocess
import multiprocessing
        
        
# Set the paths to the input directories and the output directory
input_protein_pdbqt = "input_protein_pdbqt"
input_ligand_pdbqt = "input_ligand_pdbqt"
output_directory = "affinity_maps"


import os
import subprocess
import multiprocessing

def process_file(receptor_file, receptor_dir, ligand_dir, output_dir):
    """
    Process a single receptor file and generate the affinity map using AutodockFR.

    Parameters:
        receptor_file (str): The name of the receptor file to process.
        receptor_dir (str): Path to the directory containing receptor protein files in pdbqt format.
        ligand_dir (str): Path to the directory containing ligand files in pdbqt format.
        output_dir (str): Path to the directory where affinity map files (.trg) and log files (.log) will be saved.
    """
    # Extract the protein name (without the _protein suffix and .pdbqt extension)
    protein_name = os.path.splitext(receptor_file)[0].replace("_protein", "")

    # Construct the full paths to the receptor and ligand files
    receptor_path = os.path.join(receptor_dir, receptor_file)
    ligand_path = os.path.join(ligand_dir, f"{protein_name}_ligand.pdbqt")

    # Construct the output path for the .trg file
    target_file = os.path.join(output_dir, f"{protein_name}.trg")

    # Construct the command to run the agfr tool with the output path
    command = f"agfr -r {receptor_path} -l {ligand_path} -o {target_file}"

    # Run the command using subprocess
    subprocess.run(command, shell=True)

def generate_affinity_maps(receptor_dir="input_protein_pdbqt", ligand_dir="input_ligand_pdbqt", output_dir="input_affinity_maps"):
    """
    Generate affinity maps and log files for molecular docking using AutodockFR.

    Parameters:
        receptor_dir (str, optional): Path to the directory containing receptor protein files in pdbqt format.
                                      Default is "input_protein_pdbqt".
        ligand_dir (str, optional): Path to the directory containing ligand files in pdbqt format.
                                    Default is "input_ligand_pdbqt".
        output_dir (str, optional): Path to the directory where affinity map files (.trg) and log files (.log) will be saved.
                                    Default is "input_affinity_maps".
    """
    # Get a list of receptor files in the input_protein_pdbqt directory
    receptor_files = [f for f in os.listdir(receptor_dir) if f.endswith(".pdbqt")]

    # Make sure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Use multiprocessing.Pool to process receptor files in parallel
    num_processes = int(multiprocessing.cpu_count() * 0.9)  # Get the number of CPU cores
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Use starmap to pass multiple arguments to process_file function
        pool.starmap(process_file, [(receptor_file, receptor_dir, ligand_dir, output_dir) for receptor_file in receptor_files])   


# Call the function to generate affinity maps in parallel
generate_affinity_maps(input_protein_pdbqt, input_ligand_pdbqt, output_directory)
