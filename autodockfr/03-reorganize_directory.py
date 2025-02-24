### Reorganize the working space
import os
import shutil

def move_files_to_folders():
    """
    Moves output files to their respective folders.

    Returns:
        None
    """
    # Move files from split_pdb_structure() output to extracted_pdb folder
    split_output_directory = 'extracted_pdb'
    os.makedirs(split_output_directory, exist_ok=True)

    split_files = [file for file in os.listdir() if file.endswith('_protein.pdb') or file.endswith('_ligand.pdb')]
    split_duplicates = []
    for file in split_files:
        destination_file = os.path.join(split_output_directory, file)
        if os.path.exists(destination_file):
            split_duplicates.append(file)

    if split_duplicates:
        print("The following files already exist in the 'extracted_pdb' folder:")
        for file in split_duplicates:
            print(f"- {file}")
        confirm = input("Do you want to replace these files? (y/n) ")
        if confirm.lower() != 'y':
            print("Skipping moving files to 'extracted_pdb' folder.")
        else:
            for file in split_files:
                if file in split_duplicates:
                    shutil.move(file, os.path.join(split_output_directory, file))
    else:
        for file in split_files:
            shutil.move(file, os.path.join(split_output_directory, file))

    # Move files from prepare_protein_with_reduce() output to reduced_pdb folder
    reduce_output_directory = 'reduced_pdb'
    os.makedirs(reduce_output_directory, exist_ok=True)

    reduce_files = [file for file in os.listdir() if file.endswith('_proteinFH.pdb')]
    reduce_duplicates = []
    for file in reduce_files:
        destination_file = os.path.join(reduce_output_directory, file)
        if os.path.exists(destination_file):
            reduce_duplicates.append(file)

    if reduce_duplicates:
        print("The following files already exist in the 'reduced_pdb' folder:")
        for file in reduce_duplicates:
            print(f"- {file}")
        confirm = input("Do you want to replace these files? (y/n) ")
        if confirm.lower() != 'y':
            print("Skipping moving files to 'reduced_pdb' folder.")
        else:
            for file in reduce_files:
                if file in reduce_duplicates:
                    shutil.move(file, os.path.join(reduce_output_directory, file))
    else:
        for file in reduce_files:
            shutil.move(file, os.path.join(reduce_output_directory, file))

    # Move files from convert_protein_to_pdbqt() output to input_protein_pdbqt folder
    protein_pdbqt_output_directory = 'input_protein_pdbqt'
    os.makedirs(protein_pdbqt_output_directory, exist_ok=True)

    protein_pdbqt_files = [file for file in os.listdir() if file.endswith('_proteinFH.pdbqt')]
    protein_pdbqt_duplicates = []
    for file in protein_pdbqt_files:
        destination_file = os.path.join(protein_pdbqt_output_directory, file)
        if os.path.exists(destination_file):
            protein_pdbqt_duplicates.append(file)

    if protein_pdbqt_duplicates:
        print("The following files already exist in the 'input_protein_pdbqt' folder:")
        for file in protein_pdbqt_duplicates:
            print(f"- {file}")
        confirm = input("Do you want to replace these files? (y/n) ")
        if confirm.lower() != 'y':
            print("Skipping moving files to 'input_protein_pdbqt' folder.")
        else:
            for file in protein_pdbqt_files:
                if file in protein_pdbqt_duplicates:
                    shutil.move(file, os.path.join(protein_pdbqt_output_directory, file))
    else:
        for file in protein_pdbqt_files:
            shutil.move(file, os.path.join(protein_pdbqt_output_directory, file))

    # Move files with '_ligand' suffix to input_ligand_pdbqt folder
    ligand_pdbqt_output_directory = 'input_ligand_pdbqt'
    os.makedirs(ligand_pdbqt_output_directory, exist_ok=True)

    ligand_pdbqt_files = [file for file in os.listdir() if file.endswith('_ligand.pdbqt')]
    ligand_pdbqt_duplicates = []
    for file in ligand_pdbqt_files:
        destination_file = os.path.join(ligand_pdbqt_output_directory, file)
        if os.path.exists(destination_file):
            ligand_pdbqt_duplicates.append(file)

    if ligand_pdbqt_duplicates:
        print("The following files already exist in the 'input_ligand_pdbqt' folder:")
        for file in ligand_pdbqt_duplicates:
            print(f"- {file}")
        confirm = input("Do you want to replace these files? (y/n) ")
        if confirm.lower() != 'y':
            print("Skipping moving files to 'input_ligand_pdbqt' folder.")
        else:
            for file in ligand_pdbqt_files:
                if file in ligand_pdbqt_duplicates:
                    shutil.move(file, os.path.join(ligand_pdbqt_output_directory, file))
    else:
        for file in ligand_pdbqt_files:
            shutil.move(file, os.path.join(ligand_pdbqt_output_directory, file))
            
    # Move files with '_dock' suffix to input_dock_pdbqt folder
    dock_pdbqt_output_directory = 'input_dock_pdbqt'
    os.makedirs(dock_pdbqt_output_directory, exist_ok=True)

    dock_pdbqt_files = [file for file in os.listdir() if file.endswith('_dock.pdbqt')]
    dock_pdbqt_duplicates = []
    for file in dock_pdbqt_files:
        destination_file = os.path.join(dock_pdbqt_output_directory, file)
        if os.path.exists(destination_file):
            dock_pdbqt_duplicates.append(file)

    if dock_pdbqt_duplicates:
        print("The following files already exist in the 'input_dock_pdbqt' folder:")
        for file in dock_pdbqt_duplicates:
            print(f"- {file}")
        confirm = input("Do you want to replace these files? (y/n) ")
        if confirm.lower() != 'y':
            print("Skipping moving files to 'input_dock_pdbqt' folder.")
        else:
            for file in dock_pdbqt_files:
                if file in dock_pdbqt_duplicates:
                    shutil.move(file, os.path.join(dock_pdbqt_output_directory, file))
    else:
        for file in dock_pdbqt_files:
            shutil.move(file, os.path.join(dock_pdbqt_output_directory, file))

    print("Files have been moved successfully.")

move_files_to_folders()