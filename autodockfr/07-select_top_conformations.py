import os
import logging
import shutil

# Constants
BEST_LABEL = "best"
WORST_LABEL = "worst"
TOP_N_MODELS = 10  # Number of top models to consider

def read_summary_file(summary_file_path):
    """
    Reads the summary file and returns a list of lines excluding the header.

    Args:
        summary_file_path (str): Path to the summary file.

    Returns:
        list: List of lines from the summary file.
    """
    with open(summary_file_path, 'r') as summary_file:
        return summary_file.readlines()[1:]


def extract_top_n_models(lines, label, top_n=1):
    """
    Extracts the top N models based on binding affinity.

    Args:
        lines (list): List of lines from the summary file.
        label (str): Label indicating whether to extract the best or worst models.
        top_n (int): Number of top models to extract.

    Returns:
        list: List containing tuples (receptor, ligand, affinity, model_index) for the top N models.
    """
    models = []

    for line in lines:
        columns = line.strip().split()
        receptor, ligand, model_index, affinity = columns[0], columns[1], int(columns[2]), float(columns[3])
        models.append((receptor, ligand, affinity, model_index))

    # Sort models based on affinity
    if label == BEST_LABEL:
        models = sorted(models, key=lambda x: x[2])
    elif label == WORST_LABEL:
        models = sorted(models, key=lambda x: x[2], reverse=True)

    # Extract top N models based on the label
    if label == BEST_LABEL:
        top_n_models = models[:top_n]
    elif label == WORST_LABEL:
        top_n_models = models[:top_n]
    else:
        top_n_models = []

    return top_n_models


def extract_model_data(pdbqt_file_path, model_index):
    """
    Extracts the model data from the PDBQT file based on the given model_index.

    Args:
        pdbqt_file_path (str): Path to the PDBQT file.
        model_index (int): Index of the model to extract.

    Returns:
        list: List containing the extracted model data lines.
    """
    with open(pdbqt_file_path, 'r') as pdbqt_file:
        pdbqt_data = pdbqt_file.readlines()

    start_idx, end_idx = None, None

    for idx, line in enumerate(pdbqt_data):
        if line.startswith("MODEL") and int(line.split()[1]) == model_index:
            start_idx = idx
        elif line.startswith("ENDMDL") and start_idx is not None:
            end_idx = idx
            break

    if start_idx is None or end_idx is None:
        raise ValueError(f"Invalid PDBQT file format or model_index for ligand {receptor}-{ligand_rename}.")

    return pdbqt_data[start_idx:end_idx + 1]


def save_model_data(output_file_path, mode_data):
    """
    Saves the extracted model data to the specified output file.

    Args:
        output_file_path (str): Path to the output file.
        mode_data (list): List containing the extracted model data lines.

    Returns:
        None
    """
    with open(output_file_path, 'w') as output_file:
        output_file.writelines(mode_data)


def copy_and_rename_files(source_dir, dest_dir, file_list, label):
    """
    Copies files from source directory to destination directory and renames them.

    Args:
        source_dir (str): Source directory path.
        dest_dir (str): Destination directory path.
        file_list (list): List of files to copy and rename.
        label (str): Label indicating the type of files (e.g., best or worst).

    Returns:
        None
    """
    os.makedirs(dest_dir, exist_ok=True)

    for idx, (receptor, ligand, affinity, model_index) in enumerate(file_list):
        ligand_rename = ligand.replace('_dock', '')

        # Copy and rename PDBQT file for the model
        source_file_path = os.path.join(source_dir, f"{receptor}-{ligand_rename}_out.pdbqt")
        dest_file_path = os.path.join(dest_dir, f"{receptor}-{ligand_rename}_{label}_{str(idx + 1).zfill(2)}.pdbqt")
        
        try:
            shutil.copy(source_file_path, dest_file_path)
            logging.info(f"File copied and renamed: {dest_file_path}")
        except FileNotFoundError:
            logging.error(f"Error: Source file not found: {source_file_path}")

        # Copy and rename receptor file for ligand (assuming receptor files are in "input_protein_pdbqt")
        receptor_file_path = os.path.join("input_protein_pdbqt", f"{receptor}_protein.pdbqt")
        dest_receptor_file_path = os.path.join(dest_dir, f"{receptor}_{label}_{str(idx + 1).zfill(2)}.pdbqt")
        
        try:
            shutil.copy(receptor_file_path, dest_receptor_file_path)
            logging.info(f"Receptor file copied and renamed: {dest_receptor_file_path}")
        except FileNotFoundError:
            logging.error(f"Error: Receptor file not found: {receptor_file_path}")


def main():
    # Specify input arguments
    summary_file_path = 'summary_binding_score.txt'
    pdbqt_dir = 'output_dock_pdbqt'

    # Output directories
    output_best_dir = "output_best_pdbqt"
    output_worst_dir = "output_worst_pdbqt"

    # Create the output directories if they do not exist
    os.makedirs(output_best_dir, exist_ok=True)
    os.makedirs(output_worst_dir, exist_ok=True)

    # Read the summary file
    lines = read_summary_file(summary_file_path)

    # Extract the top 10 best models based on affinity
    top_10_best_models = extract_top_n_models(lines, BEST_LABEL, TOP_N_MODELS)

    # Extract and save top 10 best models
    for idx, (receptor, ligand, affinity, model_index) in enumerate(top_10_best_models):
        ligand_rename = ligand.replace('_dock', '')
        pdbqt_file_path = os.path.join(pdbqt_dir, f"{receptor}-{ligand_rename}_out.pdbqt")
        output_file_path = os.path.join(output_best_dir, f"{receptor}-{ligand_rename}_{BEST_LABEL}_{str(idx + 1).zfill(2)}.pdbqt")

        mode_data = extract_model_data(pdbqt_file_path, model_index)
        save_model_data(output_file_path, mode_data)

        logging.info(f"{BEST_LABEL.capitalize()} binding ligand for receptor {receptor}: {ligand}, "
                     f"Affinity: {affinity:.2f}, Model: {model_index}")
        logging.info(f"Model extracted and saved to: {output_file_path}")

    # Extract the top 10 worst models based on affinity
    top_10_worst_models = extract_top_n_models(lines, WORST_LABEL, TOP_N_MODELS)

    # Extract and save top 10 worst models
    for idx, (receptor, ligand, affinity, model_index) in enumerate(top_10_worst_models):
        ligand_rename = ligand.replace('_dock', '')
        pdbqt_file_path = os.path.join(pdbqt_dir, f"{receptor}-{ligand_rename}_out.pdbqt")
        output_file_path = os.path.join(output_worst_dir, f"{receptor}-{ligand_rename}_{WORST_LABEL}_{str(idx + 1).zfill(2)}.pdbqt")

        mode_data = extract_model_data(pdbqt_file_path, model_index)
        save_model_data(output_file_path, mode_data)

        logging.info(f"{WORST_LABEL.capitalize()} binding ligand for receptor {receptor}: {ligand}, "
                     f"Affinity: {affinity:.2f}, Model: {model_index}")
        logging.info(f"Model extracted and saved to: {output_file_path}")

    # Copy and rename files for the top 10 best models
    copy_and_rename_files(pdbqt_dir, "top_receptor", top_10_best_models, BEST_LABEL)

    # Copy and rename files for the top 10 worst models
    copy_and_rename_files(pdbqt_dir, "top_receptor", top_10_worst_models, WORST_LABEL)


if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()