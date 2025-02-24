# *** Homology Modeling with MODELLER ***

# Script version 3.0
# Last updated: 18 March 2024

# Created by Nuttapat Anuwongcharoen, Ph.D.
# Department of Community Medical Technology, Mahidol University.

#################################################################

# Import essential packages
import sys
import os
import pickle
from modeller import *

# Import core_func.py in script_main
sys.path.append('../script_main')
from core_func import *

# *** Step 1: Specify input variables *** #

# Define the scope of homology modeling
single_template_modeling = True
ligand_presence_single = True
loop_model_single = True

multi_template_modeling = False
ligand_presence_multiple = True
loop_model_multiple = True

# Identify template for homology modeling
template_auto_detect = True
template_auto_alignment = True

# Specify the index number of starting and ending models
single_start_index = 1  # index of the first model
single_end_index = 100  # index of the last model (number of single-template models)
single_loop_start_index = 1  # index of the first model
single_loop_end_index = 10  # index of the last model (number of loop models for single-template)

mult_start_index = 1  # index of the first model
mult_end_index = 100  # index of the last model (number of multi-template models)
mult_loop_start_index = 1  # index of the first model
mult_loop_end_index = 10  # index of the last model (number of loop models for multi-template)

# Sequence database filename
database_search = False
seq_database = 'pdball.pir'

# Target filename in PIR format
target_seq_code = 'hkkp'

# Template code and chains for single template modeling
template_single_code = '4biw'
template_single_chain = 'A'

# Template code and chains for multiple template modeling
template_manual = [('kp_alphafold', 'A'), ('4biw', 'A')]

# Specify the desired number of CPUs
num_cpus = calculate_num_cpus(0.95)  # use 95% of available CPU cores
print(num_cpus)

######################################################################

# *** Step 2: Set alias variables *** #

# Set alias related to the target
target_seq_file = target_seq_code + '.ali'
target_single_profile = target_seq_code + '.profile'
target_multiple_profile = target_seq_code + '_mult.profile'
target_loop_profile = target_seq_code + '_loop.profile'
target_ligand_profile = target_seq_code + '_ligand.profile'

# Set alias related to the template
template_single = template_single_code + template_single_chain
template_single_file = template_single_code + '.pdb'
template_single_profile = template_single_code + '.profile'
template_seq_code = template_single

# Set alias related to the output
ok_models_single_pickle = 'ok_models_single.pickle'
ok_models_mult_pickle = 'ok_models_mult.pickle'
ok_models_loop_pickle = 'ok_models_loop.pickle'
mtop_single_pickle = 'mtop_single.pickle'
mtop_multiple_pickle = 'mtop_multiple.pickle'
mtop_loop_pickle = 'mtop_loop.pickle'
template_multiple_tuple_pickle = "template_multiple_tuple.pickle"

# Specify the input filename as variables
# aln_single_ali = target_seq_code + '_' + template_single + '.ali'
# aln_single_pap = target_seq_code + '_' + template_single + '.pap'
aln_multiple_ali = target_seq_code + '_mult.ali'

######################################################################


# *** Step 3: Perform target-template auto-alignment *** #

if template_auto_alignment:
    # Alignment for single template modeling
    if single_template_modeling:

        env = Environ()  # Create a new Modeller environment to build this model
        aln = Alignment(env)  # Read in the target sequence/alignment

        # specify a template structure (in complex with ligand)
        mdl = Model(env,
                    file=template_single_code,
                    model_segment=('FIRST:A', 'LAST:A'))

        aln.append_model(mdl,
                         align_codes=template_single,
                         atom_files=template_single_file)

        # specify a target sequence
        aln.append(file=target_seq_file, align_codes=target_seq_code)

        # Perform sequence alignment
        aln.align2d(max_gap_length=50)

        # Specify the output filename as variables
        aln_single_ali = target_seq_code + '_' + template_single + '.ali'
        aln_single_pap = target_seq_code + '_' + template_single + '.pap'

        # write the output file
        aln.write(file=aln_single_ali, alignment_format='PIR')
        aln.write(file=aln_single_pap, alignment_format='PAP', alignment_features='INDICES HELIX BETA')

        # Process the PIR file and write the output to a new file
        if ligand_presence_single:
            print('Perform a homology modeling with ligand/HETATM information')

            aln_single_ligand_ali = target_seq_code + '_' + template_single + '_ligand' + '.ali'
            process_pir_file(aln_single_ali, aln_single_ligand_ali)

            aln_single_input = aln_single_ligand_ali

        else:
            print('Perform a homology modeling without ligand/HETATM information')

            aln_single_input = aln_single_ali

    # Alignment for multiple template modeling
    if multi_template_modeling:
        # Automatic detection of multi-template
        if not template_auto_detect:
            print('Template automatic detection is turn-off')
            print("Please specify the all the templates in ('pdb_id', 'chain_id') format")

            template_multiple = template_manual

            # Create a variable of all template IDs
            template_multiple_tuple = tuple([x + y for x, y in template_multiple])

        else:
            print('Fetching all templates in pdb format inside the current directory')
            # Automatic detect template in the working directory
            template_multiple = []

            # Get the files in the current working directory with .pdb extension
            pdb_files = [file for file in os.listdir() if file.endswith('.pdb')]

            # Extract the filename and chain ID from each pdb file
            for pdb_file in pdb_files:
                with open(pdb_file, 'r', encoding="utf8", errors='ignore') as file:
                    for line in file:
                        if line.startswith('ATOM'):
                            line_data = line.split()
                            chain_id = line_data[4]
                            template_multiple.append((pdb_file.split('.')[0], chain_id))
                            break

            print(template_multiple)

            # Create a variable of all template IDs
            template_multiple_tuple = tuple([x + y for x, y in template_multiple])

        # Exporting template_multiple_tuple
        with open("template_multiple_tuple.pickle", "wb") as f:
            pickle.dump(template_multiple_tuple, f)

    if multi_template_modeling:
        # Create a new Modeller environment to build this model
        env = Environ()

        # Directories for input atom files
        env.io.atom_files_directory = ['.', '../atom_files/']

        # Read in the target sequence/alignment
        aln = Alignment(env)

        for (code, chain) in template_multiple:
            m = Model(env, file=code, model_segment=('FIRST:' + chain, 'LAST:' + chain))

            aln.append_model(m, align_codes=code + chain, atom_files=code)

        for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                            ((1., 0.5, 1., 1., 1., 0.), False, True),
                                            ((1., 1., 1., 1., 1., 0.), True, False)):
            aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                       rr_file='$(LIB)/as1.sim.mat', overhang=30,
                       gap_penalties_1d=(-450, -50),
                       gap_penalties_3d=(0, 3),
                       gap_gap_score=0,
                       gap_residue_score=0,
                       dendrogram_file='fm00495.tree',
                       alignment_type='tree',  # If 'progresive', the tree is not
                       # computed and all structues will be
                       # aligned sequentially to the first
                       feature_weights=weights,  # For a multiple sequence alignment only
                       # the first feature needs to be non-zero
                       improve_alignment=True, fit=True, write_fit=write_fit,
                       write_whole_pdb=whole, output='ALIGNMENT QUALITY')

        # Alignment output files
        aln.write(file='fm00495.pap', alignment_format='PAP')
        aln.write(file='fm00495.ali', alignment_format='PIR')

        # Preprocess the output files
        replace_fit_pdb('fm00495.ali', 'fm00495_edited.ali')

        # Calculate quality score
        aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
                   rr_file='$(LIB)/as1.sim.mat', overhang=30,
                   gap_penalties_1d=(-450, -50),
                   gap_penalties_3d=(0, 3),
                   gap_gap_score=0,
                   gap_residue_score=0,
                   dendrogram_file='1is3A.tree',
                   alignment_type='progressive', feature_weights=[0] * 6,
                   improve_alignment=False, fit=False, write_fit=True,
                   write_whole_pdb=False, output='QUALITY')

    if multi_template_modeling:

        # Create a new Modeller environment to build this model
        env = Environ()

        # Read topology
        env.libs.topology.read(file='$(LIB)/top_heav.lib')

        # Read aligned structure(s):
        aln = Alignment(env)
        aln.append(file='fm00495_edited.ali', align_codes='all')
        aln_block = len(aln)

        # Read aligned sequence(s):
        aln.append(file=target_seq_file, align_codes=target_seq_code)

        # Structure sensitive variable gap penalty sequence-sequence alignment:
        aln.salign(output='', max_gap_length=20,
                   gap_function=True,  # to use structure-dependent gap penalty
                   alignment_type='PAIRWISE', align_block=aln_block,
                   feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
                   gap_penalties_1d=(-450, 0),
                   gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
                   similarity_flag=True)

        # Specify the output filename
        aln_multiple_ali = target_seq_code + '_mult.ali'
        aln_multiple_pap = target_seq_code + '_mult.pap'

        # Write alignment output files
        aln.write(file=aln_multiple_ali, alignment_format='PIR')
        aln.write(file=aln_multiple_pap, alignment_format='PAP', alignment_features='INDICES HELIX BETA')

        # Process the PIR file and write the output to a new file
        if ligand_presence_multiple:
            print('Perform a homology modeling with ligand/HETATM information')

            aln_multiple_ligand_ali = target_seq_code + '_mult_ligand.ali'
            process_pir_file(aln_multiple_ali, aln_multiple_ligand_ali)

            # Assign alignment file with ligand info as an input parameter
            aln_multiple_input = aln_multiple_ligand_ali

        else:
            print('Perform a homology modeling without ligand/HETATM information')

            # Assign alignment file without ligand info as an input parameter
            aln_multiple_input = aln_multiple_ali

"""
# Perform sequence alignment -> Error from infinite loop!
if template_auto_alignment:
    print('Perform automatic target-template alignment')
    # subprocess.run(["python", "01-template_alignment.py"])
    os.system("python 01-template_alignment.py")
else:
    print('Skip the automatic target-template alignment process')
"""

# Create a check point before proceed the downstream process
print('############################################################')
print('Recheck the alignment files before creating a homology model')
print('############################################################')
if not confirm_continue():
    print("Script terminated by user.")
    exit()  # Terminate the script

# *** Step 4: Create homology models *** #

if single_template_modeling:

    if ligand_presence_single:
        print('Perform a single-template homology modeling with ligand/HETATM information')
        aln_single_ligand_ali = target_seq_code + '_' + template_single + '_ligand' + '.ali'
        aln_single_input = aln_single_ligand_ali
    else:
        print('Perform a single-template homology modeling without ligand/HETATM information')
        aln_single_input = aln_single_ali

    # Specify input parameters
    alignment_file = aln_single_input
    template_code = template_single
    target_seq_code = target_seq_code
    start_index = single_start_index
    end_index = single_end_index
    include_ligand = ligand_presence_single
    num_cpus = num_cpus

    # Execute MODELLER software
    if not loop_model_single:
        print('Perform homology modeling using AutoModel')
        single_auto_model(alignment_file, template_code, target_seq_code, start_index, end_index,
                          include_ligand, num_cpus)
    else:
        print('Perform homology modeling using LoopModel')
        start_loop_index = single_loop_start_index
        end_loop_index = single_loop_end_index

        single_loop_model(alignment_file, template_code, target_seq_code, start_index, end_index, start_loop_index,
                          end_loop_index, include_ligand, num_cpus)

if multi_template_modeling:

    if ligand_presence_multiple:
        print('Perform a multi-template homology modeling with ligand/HETATM information')
        aln_multiple_ligand_ali = target_seq_code + '_mult_ligand.ali'
        aln_multiple_input = aln_multiple_ligand_ali
    else:
        print('Perform a multi-template homology modeling without ligand/HETATM information')
        aln_multiple_input = aln_multiple_ali

    # Specify input parameters
    alignment_file = aln_multiple_input
    target_seq_code = target_seq_code
    start_index = single_start_index
    end_index = single_end_index
    include_ligand = ligand_presence_single
    num_cpus = num_cpus

    # Check for the template tuple
    if os.path.isfile('template_multiple_tuple.pickle'):
        with open('template_multiple_tuple.pickle', "rb") as f:
            template_tuple = pickle.load(f)
    else:
        template_tuple = template_manual

    # Execute MODELLER software
    if not loop_model_multiple:
        print('Perform homology modeling using AutoModel')
        mult_auto_model(alignment_file, template_tuple, target_seq_code, start_index, end_index,
                        include_ligand, num_cpus)
    else:
        print('Perform homology modeling using LoopModel')
        start_loop_index = single_loop_start_index
        end_loop_index = single_loop_end_index

        mult_loop_model(alignment_file, template_tuple, target_seq_code, start_index, end_index, start_loop_index,
                        end_loop_index, include_ligand, num_cpus)

# *** Step 5: Perform molecular docking with AutodockFR *** #
import subprocess

# Change the working directory to the autodockfr folder
os.chdir("autodockfr")

# Execute the AutodockFR scripts sequentially
subprocess.run(["python", "00-summary_dope.py"])
subprocess.run(["python", "01-foldx_repair.py"])
subprocess.run(["python", "02-prepare_ligand_parallel.py"])
subprocess.run(["python", "03-reorganize_directory.py"])
subprocess.run(["python", "04-generate_affinity_map.py"])
subprocess.run(["python", "05-autodockfr.py"])
subprocess.run(["python", "06-summary_docking_results.py"])
subprocess.run(["python", "07-select_top_conformations.py"])
subprocess.run(["python", "08-convert_output_mol_mol2.py"])

# Change the working directory back to the original directory
os.chdir("..")
