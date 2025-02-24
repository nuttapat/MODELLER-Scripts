# MODELLER - Generate Multiple Homology Models for Ligand-guided Homology Modeling

This script automates the process of generating homology models using MODELLER. It supports both single-template and multi-template modeling, with options for loop refinement and ligand inclusion.

## Features

- Single-template and multi-template homology modeling.
- Automatic loop refinement using AutoLoop.
- Option to include ligand information in the modeling process.
- Parallel processing for faster model generation.
- Ranking of models based on DOPE score.
- Generation of summary output files.

## Requirements

- MODELLER
- AutodockFR
- Python 3.x

## Usage

1. **Prepare input files:**
   - **Target sequence:** Provide your target sequence in PIR format (e.g., `target.ali`).
   - **Template structures:** Provide the PDB files for your template structures.
   - **Ligand information:** If you want to include ligand information, provide the ligand PDB file.

2. **Configure script parameters:**
   - Edit the `execute_modeller.py` script to set the desired parameters, such as:
     - `single_template_modeling` and `multi_template_modeling`: Enable or disable single-template and multi-template modeling.
     - `ligand_presence_single` and `ligand_presence_multiple`: Enable or disable ligand inclusion for single-template and multi-template modeling.
     - `loop_model_single` and `loop_model_multiple`: Enable or disable loop refinement for single-template and multi-template modeling.
     - `template_auto_detect`: Enable or disable automatic detection of templates in the current directory.
     - `template_manual`: Manually specify template structures if `template_auto_detect` is disabled.
     - `num_cpus`: Set the number of CPUs to use for parallel processing.

3. **Run the script:**
   ```bash
   python execute_modeller.py

## Output
The script will generate the following output files:
   - Homology models: PDB files for the generated models.
   - Summary files: CSV files containing model evaluation scores.
   - Pickle files: Pickle files containing model objects and other data.

## Additional Notes
   - The script assumes that MODELLER and AutodockFR are installed and configured correctly.
   - The script uses the core_func.py file in the script_main subfolder for core functions.
   - The script can be further customized by modifying the core_func.py file.
