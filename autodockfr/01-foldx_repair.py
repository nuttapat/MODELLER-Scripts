import os
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor

def process_file(file_path):
    print(f"Processing file: {file_path}")
    repair_command = f"foldx --command=RepairPDB --pdb={file_path}"
    repair_args = shlex.split(repair_command)
    subprocess.run(repair_args)

def get_cpu_count():
    return os.cpu_count() or 1

def main():
    cpu_count = get_cpu_count()
    available_cpus = int(cpu_count * 0.9)

    input_file_suffix = "_protein.pdb"
    input_files = [file for file in os.listdir() if file.endswith(input_file_suffix)]

    with ThreadPoolExecutor(max_workers=available_cpus) as executor:
        executor.map(process_file, input_files)

if __name__ == "__main__":
    main()
