### Convert the output pdbqt to mol and mol2

import os
import subprocess
from openbabel import openbabel

# Check if Open Babel is available and import the module
try:
    from openbabel import openbabel as ob
    OPENBABEL_AVAILABLE = True
except ImportError:
    OPENBABEL_AVAILABLE = False


# Input and output directories
INPUT_DIRECTORY = "output_best_pdbqt"
OUTPUT_DIRECTORY_MOL = "output_best_mol"
OUTPUT_DIRECTORY_MOL2 = "output_best_mol2"


def convert_pdbqt_to_mol(input_dir, output_dir):
    """Convert .pdbqt files in the input directory to .mol format and save them in the output directory."""

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    pdbqt_files = [f for f in os.listdir(input_dir) if f.endswith(".pdbqt")]
    converted_files = []

    if OPENBABEL_AVAILABLE:
        for pdbqt_file in pdbqt_files:
            input_file_path = os.path.join(input_dir, pdbqt_file)
            output_file_path = os.path.join(output_dir, os.path.splitext(pdbqt_file)[0] + ".mol")

            convert_file_using_openbabel(input_file_path, output_file_path, "pdbqt", "mol")
            converted_files.append(output_file_path)
    else:
        for pdbqt_file in pdbqt_files:
            input_file_path = os.path.join(input_dir, pdbqt_file)
            output_file_path = os.path.join(output_dir, os.path.splitext(pdbqt_file)[0] + ".mol")

            convert_file_using_subprocess(input_file_path, output_file_path, "mol")
            converted_files.append(output_file_path)

    return converted_files


def convert_pdbqt_to_mol2(input_dir, output_dir):
    """Convert .pdbqt files in the input directory to .mol2 format and save them in the output directory."""

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    pdbqt_files = [f for f in os.listdir(input_dir) if f.endswith(".pdbqt")]
    converted_files = []

    if OPENBABEL_AVAILABLE:
        for pdbqt_file in pdbqt_files:
            input_file_path = os.path.join(input_dir, pdbqt_file)
            output_file_path = os.path.join(output_dir, os.path.splitext(pdbqt_file)[0] + ".mol2")

            convert_file_using_openbabel(input_file_path, output_file_path, "pdbqt", "mol2")
            converted_files.append(output_file_path)
    else:
        for pdbqt_file in pdbqt_files:
            input_file_path = os.path.join(input_dir, pdbqt_file)
            output_file_path = os.path.join(output_dir, os.path.splitext(pdbqt_file)[0] + ".mol2")

            convert_file_using_subprocess(input_file_path, output_file_path, "mol2")
            converted_files.append(output_file_path)

    return converted_files


def convert_file_using_openbabel(input_file_path, output_file_path, in_format, out_format):
    """Convert a file using Open Babel."""

    conv = ob.OBConversion()
    conv.SetInFormat(in_format)
    conv.SetOutFormat(out_format)
    mol = ob.OBMol()
    conv.ReadFile(mol, input_file_path)
    conv.WriteFile(mol, output_file_path)


def convert_file_using_subprocess(input_file_path, output_file_path, format_type):
    """Convert a file using the subprocess module."""

    command = f"obabel {input_file_path} -O {output_file_path} -d {format_type}"
    subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


# Convert pdbqt to mol and mol2
convert_pdbqt_to_mol(INPUT_DIRECTORY, OUTPUT_DIRECTORY_MOL)
convert_pdbqt_to_mol2(INPUT_DIRECTORY, OUTPUT_DIRECTORY_MOL2)
