import os
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor

def process_file(file_path):
    print(f"Processing file: {file_path}")
    
    # Generate the new output file name by changing the suffix
    new_output_file = file_path.replace("_protein_Repair.pdb", "_protein.pdbqt")

    # Add your second shell command with the new output file name as an argument
    prepare_receptor = f"prepare_receptor -r {file_path} -o {new_output_file} -A bonds_hydrogens"
    additional_args = shlex.split(prepare_receptor)
    subprocess.run(additional_args)

def get_cpu_count():
    return os.cpu_count() or 1

def main():
    cpu_count = get_cpu_count()
    available_cpus = int(cpu_count * 0.9)

    input_file_suffix = "_protein_Repair.pdb"
    input_files = [file for file in os.listdir() if file.endswith(input_file_suffix)]

    with ThreadPoolExecutor(max_workers=available_cpus) as executor:
        executor.map(process_file, input_files)

if __name__ == "__main__":
    main()
