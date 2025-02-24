## Summary of Results

import os
import pandas as pd
import re


# Set the path to the docking_results directory
docking_results_dir = "docking_results"


def collect_binding_scores(docking_results_dir, output_file="summary_binding_score.txt"):
    """
    Collect the best binding energy for each ligand in the docking_results directory.

    Parameters:
        docking_results_dir (str): Path to the directory containing docking results (.dlg files).
        output_file (str, optional): Path to the output file to save the summary. Default is 'summary_binding_score.txt'.
    """

    binding_scores = []
    dlg_files = [f for f in os.listdir(docking_results_dir) if f.endswith(".dlg")]

    for dlg_file in dlg_files:
        ligand_name = os.path.splitext(dlg_file)[0]

        with open(os.path.join(docking_results_dir, dlg_file), "r") as f:
            lines = f.readlines()

        start_index = None
        receptor = None
        ligand = None

        for i, line in enumerate(lines):
            if "Unpacking maps" in line:
                receptor = line.split("/")[-1].split(".")[0]
            elif "reading ligand" in line:
                ligand = line.split("/")[-1].split(".")[0]
            elif "mode |  affinity" in line:
                start_index = i + 2
                break

        if start_index is None or receptor is None or ligand is None:
            continue

        for line in lines[start_index:]:
            if re.match(r'^\s*\d', line): # only append rows that start with a number
                binding_info = line.split()
                if len(binding_info) >= 8:
                    mode = binding_info[0]
                    affinity = binding_info[1]
                    clust_rmsd = binding_info[2]
                    ref_rmsd = binding_info[3]
                    clust_size = binding_info[4]
                    rmsd_stdv = binding_info[5]
                    energy_stdv = binding_info[6]
                    best_run = binding_info[7]

                    binding_scores.append({
                        "receptor": receptor,
                        "ligand": ligand,
                        "mode": mode,
                        "affinity_(kcal/mol)": affinity,
                        "clust_rmsd": clust_rmsd,
                        "ref_rmsd": ref_rmsd,
                        "clust_size": clust_size,
                        "rmsd_stdv": rmsd_stdv,
                        "energy_stdv": energy_stdv,
                        "best_run": best_run,
                    })

    binding_scores_df = pd.DataFrame(binding_scores)
    binding_scores_df = binding_scores_df[['receptor', 'ligand', 'mode', 'affinity_(kcal/mol)', 'clust_rmsd',
                                           'ref_rmsd', 'clust_size', 'rmsd_stdv', 'energy_stdv', 'best_run']]

    binding_scores_df.to_csv(output_file, sep="\t", index=False)

def save_top_scores(summary_file="summary_binding_score.txt"):
    binding_scores_df = pd.read_csv(summary_file, sep="\t")

    top_10_best_scores = binding_scores_df.sort_values(by='affinity_(kcal/mol)', ascending=True).head(10)
    top_10_best_scores.to_csv("top-10-best-score.txt", sep="\t", index=False)

    top_10_worst_scores = binding_scores_df.sort_values(by='affinity_(kcal/mol)', ascending=False).head(10)
    top_10_worst_scores.to_csv("top-10-worst-score.txt", sep="\t", index=False)


# Collect the binding scores and save to 'summary_binding_score.txt'
collect_binding_scores(docking_results_dir)

# Save top 10 best and worst scores to separate files
save_top_scores()
