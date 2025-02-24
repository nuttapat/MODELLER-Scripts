import os
from modeller import *
from modeller.scripts import complete_pdb

# Set up the Modeller environment
env = Environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# Function to calculate the DOPE score for a given .pdb file
def calculate_dope_score(pdb_file):
    mdl = complete_pdb(env, pdb_file)
    atmsel = Selection(mdl.chains[0])
    score = atmsel.assess_dope()
    return score

# Open the summary file for writing
with open("summary_dope_score.txt", "w") as summary_file:
    # Write header
    summary_file.write("Filename\tDOPE Score\n")

    # Iterate over each .pdb file in the working directory and calculate the DOPE score
    for filename in os.listdir("."):
        if filename.endswith("_protein.pdb"):
            print(f"Processing {filename}...")
            dope_score = calculate_dope_score(filename)
            print(f"DOPE score: {dope_score}")

            # Extract filename without extension
            file_name_without_extension = os.path.splitext(filename)[0]

            # Write data to the summary file in real-time
            summary_file.write(f"{file_name_without_extension}\t{dope_score}\n")

print("Summary data saved to 'summary_dope_score.txt'.")


#### old script ####
'''
import os
from modeller import *
from modeller.scripts import complete_pdb

# Set up the Modeller environment
env = Environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# Function to calculate the DOPE score for a given .pdb file
def calculate_dope_score(pdb_file):
    mdl = complete_pdb(env, pdb_file)
    atmsel = Selection(mdl.chains[0])
    score = atmsel.assess_dope()
    return score

# Create a list to store the summary data
summary_data = []

# Iterate over each .pdb file in the working directory and calculate the DOPE score
for filename in os.listdir("."):
    if filename.endswith("_protein.pdb"):
        print(f"Processing {filename}...")
        dope_score = calculate_dope_score(filename)
        print(f"DOPE score: {dope_score}")

        # Extract filename without extension
        file_name_without_extension = os.path.splitext(filename)[0]

        # Append data to the summary list
        summary_data.append((file_name_without_extension, dope_score))

# Save the summary data to a text file in tabular format
with open("summary_dope_score.txt", "w") as summary_file:
    # Write header
    summary_file.write("Filename\tDOPE Score\n")

    # Write data
    for entry in summary_data:
        summary_file.write(f"{entry[0]}\t{entry[1]}\n")

print("Summary data saved to 'summary_dope_score.txt'.")

'''