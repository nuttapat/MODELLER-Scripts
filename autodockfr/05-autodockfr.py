import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
import shutil

#### Molecular docking and virtual screening
def perform_molecular_docking_parallel(ligand_file, affinity_map_file, output_dir, nb_runs, max_evals, no_improve_stop, max_gens, seed_value, adfr_options='--overwriteFiles'):
    ligand_path = os.path.join(ligand_dir, ligand_file)
    ligand_name = os.path.splitext(ligand_file)[0]
    ligand_short = os.path.splitext(ligand_file)[0].replace("_dock", "")

    affinity_map_path = os.path.join(affinity_map_dir, affinity_map_file)
    affinity_map_name = os.path.splitext(affinity_map_file)[0]
    output_prefix = os.path.join(output_dir, f"{affinity_map_name}-{ligand_name}")

    # Construct the adfr command
    adfr_command = f"adfr -l {ligand_path} -t {affinity_map_path} -o {output_prefix} --nbRuns {nb_runs} --maxEvals {max_evals} --noImproveStop {no_improve_stop} --maxGens {max_gens} --seed {seed_value} {adfr_options}"

    # Print information and perform molecular docking using adfr
    print(f"##############################")
    print(f"Performing molecular docking for {ligand_name} with {affinity_map_name}")
    print(f"  - ligand_file       = {ligand_file}")
    print(f"  - ligand_path       = {ligand_path}")
    print(f"  - ligand_name       = {ligand_name}")
    print(f"  - ligand_short      = {ligand_short}")
    print(f"  - affinity_map_file = {affinity_map_file}")
    print(f"  - affinity_map_path = {affinity_map_path}")
    print(f"  - affinity_map_name = {affinity_map_name}")
    print(f"  - output_prefix     = {output_prefix}")
    subprocess.run(adfr_command, shell=True)

    # Rename the output files and move them to their respective directories
    output_pdbqt = f"{output_prefix}_out.pdbqt"
    output_dlg = f"{output_prefix}_summary.dlg"
    output_dro = f"{output_prefix}.dro"

    os.rename(f"{output_prefix}_out.pdbqt", output_pdbqt)
    os.rename(f"{output_prefix}_summary.dlg", output_dlg)
    os.rename(f"{output_prefix}.dro", output_dro)

    output_pdbqt_dir = "output_dock_pdbqt"
    os.makedirs(output_pdbqt_dir, exist_ok=True)
    shutil.move(output_pdbqt, os.path.join(output_pdbqt_dir, f"{affinity_map_name}-{ligand_short}_out.pdbqt"))

    docking_objects_dir = os.path.join(output_dir, "docking_objects")
    os.makedirs(docking_objects_dir, exist_ok=True)
    shutil.move(output_dro, os.path.join(docking_objects_dir, f"{affinity_map_name}-{ligand_short}.dro"))

def perform_molecular_docking(ligand_dir, affinity_map_dir, output_dir, nb_runs, max_evals, no_improve_stop, max_gens, seed_value, adfr_options='--overwriteFiles'):
    # Get the list of ligand files and affinity map files
    ligand_files = [file for file in os.listdir(ligand_dir) if file.endswith('_dock.pdbqt')]
    affinity_map_files = [file for file in os.listdir(affinity_map_dir) if file.endswith('.trg')]

    # Check if ligand and affinity map files are available
    if not ligand_files:
        raise ValueError(f"No ligand files found in the ligand directory ({ligand_dir}).")
    if not affinity_map_files:
        raise ValueError(f"No affinity map files found in the affinity map directory ({affinity_map_dir}).")

    # Make sure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Specify the number of parallel jobs
    #num_workers = int(os.cpu_count() * 0.9)  # use 90% of CPU cores
    num_workers =  int((os.cpu_count() * 0.9)/2)
    
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        for ligand_file in ligand_files:
            for affinity_map_file in affinity_map_files:
                executor.submit(perform_molecular_docking_parallel, ligand_file, affinity_map_file, output_dir, nb_runs, max_evals, no_improve_stop, max_gens, seed_value, adfr_options)

# Specify the input/output files directory
ligand_dir = "input_dock_pdbqt"
affinity_map_dir = "affinity_maps"
output_dir = "docking_results"

# Define the docking parameters
nb_runs = 50
max_evals = 2500000
no_improve_stop = 5
max_gens = 10000000
seed_value = 8
adfr_options = "--maxCores 2 --overwriteFiles"

perform_molecular_docking(ligand_dir, affinity_map_dir, output_dir, nb_runs, max_evals, no_improve_stop, max_gens, seed_value, adfr_options)
