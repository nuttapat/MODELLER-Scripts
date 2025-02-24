# Import MODELLER package and call required functions
import os
import re
import pickle
import shutil
import pandas as pd
from modeller import *
from modeller.automodel import *
from modeller.parallel import Job, LocalWorker


# Calculate number of CPUs for parallel computing
def calculate_num_cpus(percentage=0.95):
    """
    Calculates the number of CPUs to use for parallel computing, using a
    percentage of the total available CPUs.

    Args:
        percentage (float, optional): Percentage of total CPUs to use. Defaults to 0.95.

    Returns:
        int: The number of CPUs to use for parallel computing.
    """

    num_cpus = int(os.cpu_count() * percentage)

    # Ensure at least one CPU is used
    if num_cpus < 1:
        num_cpus = 1

    return num_cpus


# Generates the MODELLER command with the specified number of workers
def parallel_setup(num_cpus):
    """
    This function generates a MODELLER command that utilizes multiple CPUs for parallel processing.

    Args:
        num_cpus (int): The desired number of CPUs to use in the parallel setup.

    Returns:
        str: The formatted MODELLER command string.
    """

    command = "# Use {} CPUs in a parallel job on this machine\n".format(num_cpus)
    command += "j = Job()\n"

    for i in range(num_cpus):
        command += "j.append(LocalWorker())\n"

    return command


# Create a script check point
def confirm_continue():
    """Asks the user for continuation, returns True if 'yes', False otherwise."""
    while True:
        answer = input("Do you want to continue? (yes/no): ").lower()
        if answer == 'yes':
            return True
        elif answer == 'no':
            return False
        else:
            print("Invalid input. Please enter 'yes' or 'no'.")


# Remove suffix of alignment header in multiple sequence alignment
def replace_fit_pdb(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith("structureX"):
                line = line.replace("_fit.pdb", ".pdb")
            f_out.write(line)


# Adding ligand information into the alignment file
def process_pir_file(input_file, output_file):
    # Read the input file
    with open(input_file, 'r', encoding="utf8", errors='ignore') as file:
        lines = file.readlines()

    # Find the line starting with "structureX"
    pattern = r'^structureX.*$'
    for i, line in enumerate(lines):
        if re.match(pattern, line):
            # Find the fifth column and increment its value by 1
            columns = line.split(':')
            if len(columns) >= 5:
                try:
                    number = columns[4].lstrip('+')  # Remove the leading "+"
                    number = int(number)
                    columns[4] = '+' + str(number + 1)  # Add the "+" back
                    lines[i] = ':'.join(columns)

                except ValueError:
                    print('Invalid number in the fifth column.')

    # Add "." before "*" in each line
    for i, line in enumerate(lines):
        lines[i] = line.replace('*', '.*')

    # Check if the output file already exists
    if os.path.exists(output_file):
        # Prompt the user for confirmation to replace the file
        user_input = input(f"File '{output_file}' already exists. Do you want to replace it? (y/n): ")
        if user_input.lower() != 'y':
            print('Operation cancelled. Output file was not replaced.')
            return

    # Write the modified lines to the output file
    with open(output_file, 'w', encoding="utf8", errors='ignore') as file:
        file.writelines(lines)

    print(f"Output file '{output_file}' has been created.")


# Homology modeling of single template model
def single_auto_model(alignment_file, template_code, target_seq_code, start_index, end_index, include_ligand, num_cpus):
    """
    This function performs homology modeling using Modeller based-on single template.
    It takes input parameters, generates models, ranks them based on DOPE score, and saves the top model.

    Args:
        alignment_file (str): Alignment input filename.
        template_code (str): Code of the templates in alignment the file.
        target_seq_code (str): Code of target in the alignment file.
        start_index (int): Index of the first model to generate.
        end_index (int): Index of the last model to generate.
        include_ligand (bool): Whether to include HETATM records.
        num_cpus (int): Number of CPUs to use for parallel processing.
    """

    # Step 1: Parallel Configuration Setup
    j = Job()
    for _ in range(num_cpus):
        j.append(LocalWorker())

    # Step 2: Modeller Environment Setup
    env = Environ()  # Create a new Modeller environment
    env.io.atom_files_directory = ['.', '../atom_files/']  # Set input atom file directories

    if include_ligand:
        env.io.hetatm = True  # Read HETATM records from template PDBs

    # Step 3: Model Generation
    a = AutoModel(env,
                  alnfile=alignment_file,
                  knowns=template_code,
                  sequence=target_seq_code,
                  assess_methods=(assess.DOPE, assess.GA341))

    a.starting_model = start_index
    a.ending_model = end_index

    a.use_parallel_job(j)  # Enable parallel processing

    a.md_level = refine.slow  # Set a slow MD optimization level

    a.library_schedule = autosched.slow  # Set a slow optimization schedule
    a.max_var_iterations = 300  # Increase maximum variable iterations
    a.repeat_optimization = 3  # Repeat optimization 3 times
    a.max_molpdf = 1e6  # Set a maximum objective function value

    a.make()  # Build the models

    # Step 4: Summarize the outputs
    ok_models_single = [x for x in a.outputs]  # List all generated models

    # Save all modeller outputs to pickle format
    with open('output_models_single.pickle', "wb") as f:
        pickle.dump(ok_models_single, f)

    # Export summary output to dataframe
    data_single = a.outputs

    # Convert the list of dictionaries to a dataframe
    df_single = pd.DataFrame(data_single)

    # Write the summary of model evaluating score
    df_single.to_csv('summary_single_model.csv', sep=",", header=True, index=True, index_label='model')

    # Step 5: Rank the models by DOPE score
    key = 'DOPE score'
    ok_models_single.sort(key=lambda a: a[key])

    # Get top model
    mtop_single = ok_models_single[0]
    print("Top single model: %s (DOPE score %.3f)" % (mtop_single['name'], mtop_single[key]))

    # Save the top single model to pickle format
    with open('mtop_single.pickle', "wb") as f:
        pickle.dump(mtop_single, f)

    # Rename the best model's pdb file
    shutil.copy(mtop_single['name'], 'best_single_model.pdb')


# Homology modeling of multiple template model
def mult_auto_model(alignment_file, template_tuple, target_seq_code, start_index, end_index, include_ligand,
                    num_cpus):
    """
    This function performs homology modeling using Modeller based-on multiple template.
    It takes input parameters, generates models, ranks them based on DOPE score, and saves the top model.

    Args:
        alignment_file (str): Alignment input filename.
        template_tuple (tuple): Code of the templates in alignment the file.
        target_seq_code (str): Code of target in the alignment file.
        start_index (int): Index of the first model to generate.
        end_index (int): Index of the last model to generate.
        include_ligand (bool): Whether to include HETATM records.
        num_cpus (int): Number of CPUs to use for parallel processing.
    """

    # Step 1: Parallel Configuration Setup
    j = Job()
    for _ in range(num_cpus):
        j.append(LocalWorker())

    # Step 2: Modeller Environment Setup
    env = Environ()  # Create a new Modeller environment
    env.io.atom_files_directory = ['.', '../atom_files/']  # Set input atom file directories

    if include_ligand:
        env.io.hetatm = True  # Read HETATM records from template PDBs

    if os.path.isfile('template_multiple_tuple.pickle'):
        with open('template_multiple_tuple.pickle', "rb") as f:
            template_tuple = pickle.load(f)

    # Step 3: Model Generation
    a = AutoModel(env,
                  alnfile=alignment_file,
                  knowns=template_tuple,
                  sequence=target_seq_code,
                  assess_methods=(assess.DOPE, assess.GA341))

    a.starting_model = start_index
    a.ending_model = end_index

    a.use_parallel_job(j)  # Enable parallel processing

    a.md_level = refine.slow  # Set a slow MD optimization level

    a.library_schedule = autosched.slow  # Set a slow optimization schedule
    a.max_var_iterations = 300  # Increase maximum variable iterations
    a.repeat_optimization = 3  # Repeat optimization 3 times
    a.max_molpdf = 1e6  # Set a maximum objective function value

    a.make()  # Build the models

    # Step 4: Summarize the outputs
    ok_models_mult = [x for x in a.outputs]  # List all generated models

    # Save all modeller outputs to pickle format
    with open('output_models_single.pickle', "wb") as f:
        pickle.dump(ok_models_mult, f)

    # Export summary output to dataframe
    data_mult = a.outputs

    # Convert the list of dictionaries to a dataframe
    df_mult = pd.DataFrame(data_mult)

    # Write the summary of model evaluating score
    df_mult.to_csv('summary_mult_model.csv', sep=",", header=True, index=True, index_label='model')

    # Step 5: Rank the models by DOPE score
    key = 'DOPE score'
    ok_models_mult.sort(key=lambda a: a[key])

    # Get top model
    mtop_mult = ok_models_mult[0]
    print("Top multi-template model: %s (DOPE score %.3f)" % (mtop_mult['name'], mtop_mult[key]))

    # Save the top single model to pickle format
    with open('mtop_mult.pickle', "wb") as f:
        pickle.dump(mtop_mult, f)

    # Rename the best model's pdb file
    shutil.copy(mtop_mult['name'], 'best_mult_model.pdb')


# Homology modeling of single-template model with AutoLoop Refinement
def single_loop_model(alignment_file, template_code, target_seq_code, start_index, end_index, start_loop_index,
                      end_loop_index, include_ligand, num_cpus):
    """
    This function performs homology modeling using Modeller based on single template and performs an automatic
    loop refinement. It takes input parameters, generates models, ranks them based on DOPE score, and saves
    the top model.

    Args:
        alignment_file (str): Alignment input filename.
        template_code (str): Code of the templates in alignment the file.
        target_seq_code (str): Code of target in the alignment file.
        start_index (int): Index of the first core model to generate.
        end_index (int): Index of the last core model to generate.
        start_loop_index (int): Index of the first loop model to generate.
        end_loop_index (int): Index of the last loop model to generate.
        include_ligand (bool): Whether to include HETATM records.
        num_cpus (int): Number of CPUs to use for parallel processing.
    """

    # Step 1: Parallel Configuration Setup
    j = Job()
    for _ in range(num_cpus):
        j.append(LocalWorker())

    # Step 2: Modeller Environment Setup
    env = Environ()  # Create a new Modeller environment
    env.io.atom_files_directory = ['.', '../atom_files/']  # Set input atom file directories

    if include_ligand:
        env.io.hetatm = True  # Read HETATM records from template PDBs

    # Step 3: Model Generation
    a = LoopModel(env,
                  alnfile=alignment_file,
                  knowns=template_code,
                  sequence=target_seq_code,
                  assess_methods=(assess.DOPE, assess.GA341),
                  loop_assess_methods=(assess.DOPE, assess.GA341))

    a.starting_model = start_index
    a.ending_model = end_index

    a.loop.starting_model = start_loop_index
    a.loop.ending_model = end_loop_index

    a.use_parallel_job(j)  # Enable parallel processing

    a.md_level = refine.slow  # Set a slow MD optimization level

    a.library_schedule = autosched.slow  # Set a slow optimization schedule
    a.max_var_iterations = 300  # Increase maximum variable iterations
    a.repeat_optimization = 3  # Repeat optimization 3 times
    a.max_molpdf = 1e6  # Set a maximum objective function value

    a.make()  # Build the models

    # Step 4: Summarize the outputs
    ok_models_single = [x for x in a.outputs]  # List all generated core models
    ok_models_loop = [x for x in a.loop.outputs]  # List all generated loop models
    ok_models = ok_models_single + ok_models_loop  # Combine all results

    # Save all modeller outputs to pickle format
    with open('output_loop_single.pickle', "wb") as f:
        pickle.dump(ok_models, f)

    # Export summary output to dataframe
    data_single = a.outputs
    data_loop = a.loop.outputs

    # Update the value of 'failure' key in each dictionary to string
    # for key in data_single:
    #    key['failure'] = str(key['failure'])
    # for key in data_loop:
    #    key['failure'] = str(key['failure'])

    # Convert the list of dictionaries to a dataframe
    df_single = pd.DataFrame(data_single)
    df_loop = pd.DataFrame(data_loop)

    # Write the summary of model evaluating score
    df_single.to_csv('summary_single_model.csv', sep=",", header=True, index=True, index_label='model')
    df_loop.to_csv('summary_loop_single_model.csv', sep=",", header=True, index=True, index_label='model')

    # Step 5: Rank the models by DOPE score
    key = 'DOPE score'
    ok_models.sort(key=lambda a: a[key])

    # Get top model
    mtop_single = ok_models[0]
    print("Top single model: %s (DOPE score %.3f)" % (mtop_single['name'], mtop_single[key]))

    # Save the top single model to pickle format
    with open('mtop_single.pickle', "wb") as f:
        pickle.dump(mtop_single, f)

    # Rename the best model's pdb file
    shutil.copy(mtop_single['name'], 'best_single_model.pdb')


# Homology modeling of multi-template model with AutoLoop Refinement
def mult_loop_model(alignment_file, template_tuple, target_seq_code, start_index, end_index, start_loop_index,
                    end_loop_index, include_ligand, num_cpus):
    """
    This function performs homology modeling using Modeller based on multiple templates and performs an automatic
    loop refinement. It takes input parameters, generates models, ranks them based on DOPE score, and saves
    the top model.

    Args:
        alignment_file (str): Alignment input filename.
        template_tuple (str): Code of the templates in alignment the file.
        target_seq_code (str): Code of target in the alignment file.
        start_index (int): Index of the first core model to generate.
        end_index (int): Index of the last core model to generate.
        start_loop_index (int): Index of the first loop model to generate.
        end_loop_index (int): Index of the last loop model to generate.
        include_ligand (bool): Whether to include HETATM records.
        num_cpus (int): Number of CPUs to use for parallel processing.
    """

    # Step 1: Parallel Configuration Setup
    j = Job()
    for _ in range(num_cpus):
        j.append(LocalWorker())

    # Step 2: Modeller Environment Setup
    env = Environ()  # Create a new Modeller environment
    env.io.atom_files_directory = ['.', '../atom_files/']  # Set input atom file directories

    if include_ligand:
        env.io.hetatm = True  # Read HETATM records from template PDBs

    if os.path.isfile('template_multiple_tuple.pickle'):
        with open('template_multiple_tuple.pickle', "rb") as f:
            template_tuple = pickle.load(f)

    # Step 3: Model Generation
    a = LoopModel(env,
                  alnfile=alignment_file,
                  knowns=template_tuple,
                  sequence=target_seq_code,
                  assess_methods=(assess.DOPE, assess.GA341),
                  loop_assess_methods=(assess.DOPE, assess.GA341))

    a.starting_model = start_index
    a.ending_model = end_index

    a.loop.starting_model = start_loop_index
    a.loop.ending_model = end_loop_index

    a.use_parallel_job(j)  # Enable parallel processing

    a.md_level = refine.slow  # Set a slow MD optimization level

    a.library_schedule = autosched.slow  # Set a slow optimization schedule
    a.max_var_iterations = 300  # Increase maximum variable iterations
    a.repeat_optimization = 3  # Repeat optimization 3 times
    a.max_molpdf = 1e6  # Set a maximum objective function value

    a.make()  # Build the models

    # Step 4: Summarize the outputs
    ok_models_mult = [x for x in a.outputs]  # List all generated core models
    ok_models_loop = [x for x in a.loop.outputs]  # List all generated loop models
    ok_models = ok_models_mult + ok_models_loop  # Combine all results

    # Save all modeller outputs to pickle format
    with open('output_loop_multiple.pickle', "wb") as f:
        pickle.dump(ok_models, f)

    # Export summary output to dataframe
    data_mult = a.outputs
    data_loop = a.loop.outputs

    # Update the value of 'failure' key in each dictionary to string
    # for key in data_mult:
    #    key['failure'] = str(key['failure'])
    # for key in data_loop:
    #    key['failure'] = str(key['failure'])

    # Convert the list of dictionaries to a dataframe
    df_mult = pd.DataFrame(data_mult)
    df_loop = pd.DataFrame(data_loop)

    # Write the summary of model evaluating score
    df_mult.to_csv('summary_mult_model.csv', sep=",", header=True, index=True, index_label='model')
    df_loop.to_csv('summary_loop_mult_model.csv', sep=",", header=True, index=True, index_label='model')

    # Step 5: Rank the models by DOPE score
    key = 'DOPE score'
    ok_models.sort(key=lambda a: a[key])

    # Get top model
    mtop_mult = ok_models[0]
    print("Top multiple model: %s (DOPE score %.3f)" % (mtop_mult['name'], mtop_mult[key]))

    # Save the top multiple model to pickle format
    with open('mtop_multiple.pickle', "wb") as f:
        pickle.dump(mtop_mult, f)

    # Rename the best model's pdb file
    shutil.copy(mtop_mult['name'], 'best_mult_model.pdb')
