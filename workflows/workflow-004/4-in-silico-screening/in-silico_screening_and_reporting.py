#!/usr/bin/env python
# coding: utf-8

# # Step 1 Installation of dependencies
# Gnina will run within a linux environment provided by google colab virtual machine.
# 
# 1. `useful_rdkit_utils` is a Python package written and maintained by Pat Walters that contains useful RDKit functions. We will use it for the functions `mcs_rmsd` (explained later).
# 2. `py3Dmol` is used for molecular visualization.
# 3. The RDKit is a popular cheminiformatics package we will use for processing molecules.
# 

# In[ ]:


import subprocess
import os
import requests

# Install Python packages
#subprocess.run(["pip", "install", "useful_rdkit_utils", "py3Dmol", "rdkit"])



#import json

#with open("params.json") as f:
    #data = json.load(f)

#ligand_id = data["ligand_id"]
#print("Ligand ID parsed:", ligand_id)


# In[3]:


### Check for previous gnina installation ###
gnina_file = "gnina"

# Download URL
gnina_url = "https://github.com/gnina/gnina/releases/download/v1.3/gnina"

# Check if file exists
if os.path.exists(gnina_file):
    print(f"✅ {gnina_file} already exists, skipping download.")
else:
    print(f"⬇️ {gnina_file} not found, downloading...")
    subprocess.run([
        "wget",
        gnina_url,
        "-O", gnina_file
    ], check=True)
    
    # Make gnina executable
    subprocess.run(["chmod", "+x", "gnina"])
    print("✅ Download complete.")

#pdb_id = input("Enter PDB code: ") # The Protein ID we're looking at
#ligand_id = input("Enter ligand code: ") # The ID of the co-crystallized ligand

### Working directories ###
protein_directory = "molecular_docking/protein_files"
protein_filename = f"{pdb_id}.pdb"
ligand_directory = "molecular_docking/ligand_structures"
docking_results_directory = "molecular_docking/docking_results"


# Commands for running Gnina
# 
# ```
# ./gnina \
#   # Specify the receptor structure file (-r).
#   # This file (e.g. 7LME.pdbqt) should be prepared for docking (e.g., with hydrogens added).
#   -r docking_files/7LME_all_atom.pdbqt \
#   # Specify the ligand structure file (-l) to be docked.
#   # This file (Y6J_ideal.pdbqt) contains the 3D coordinates of the ligand.
#   -l docking_files/Y6J_ideal.pdbqt \
#   # Define the docking search box automatically (--autobox_ligand).
#   # The box will be centered around the coordinates of the ligand in the specified file
#   # (Y6J_corrected_pose.sdf), which is the known experimental pose in this redocking example.
#   # An optional padding (default 4Å) is added.
#   --autobox_ligand docking_files/Y6J_corrected_pose.sdf \
#   # Specify the output file path (-o) where the resulting docked poses will be saved.
#   # The output format will be SDF, containing multiple poses ranked by score.
#   -o docking_results/Y6J_docked_e12.sdf \
#   # Set the random number generator seed (--seed) to 0.
#   # Using a fixed seed makes the docking calculation reproducible.
#   --seed 0 \
#   # Set the exhaustiveness level (--exhaustiveness) to 12.
#   # This controls the number of Monte Carlo chains for the ligand.
#   # The default is 8
#   --exhaustiveness 16
#   # Run without Convolutional Neural Network (CNN) score
#   --cnn_scoring none
#   ```
# ```
# Full command list:
# 
# Input:
#   -r [ --receptor ] arg              rigid part of the receptor
#   --flex arg                         flexible side chains, if any (PDBQT)
#   -l [ --ligand ] arg                ligand(s)
#   --flexres arg                      flexible side chains specified by comma
#                                      separated list of chain:resid
#   --flexdist_ligand arg              Ligand to use for flexdist
#   --flexdist arg                     set all side chains within specified
#                                      distance to flexdist_ligand to flexible
#   --flex_limit arg                   Hard limit for the number of flexible
#                                      residues
#   --flex_max arg                     Retain at at most the closest flex_max
#                                      flexible residues
# 
# Search space (required):
#   --center_x arg                     X coordinate of the center
#   --center_y arg                     Y coordinate of the center
#   --center_z arg                     Z coordinate of the center
#   --size_x arg                       size in the X dimension (Angstroms)
#   --size_y arg                       size in the Y dimension (Angstroms)
#   --size_z arg                       size in the Z dimension (Angstroms)
#   --autobox_ligand arg               Ligand to use for autobox. A multi-ligand
#                                      file still only defines a single box.
#   --autobox_add arg                  Amount of buffer space to add to
#                                      auto-generated box (default +4 on all six
#                                      sides)
#   --autobox_extend arg (=1)          Expand the autobox if needed to ensure the
#                                      input conformation of the ligand being
#                                      docked can freely rotate within the box.
#   --no_lig                           no ligand; for sampling/minimizing
#                                      flexible residues
# 
# Covalent docking:
#   --covalent_rec_atom arg            Receptor atom ligand is covalently bound
#                                      to.  Can be specified as
#                                      chain:resnum:atom_name or as x,y,z
#                                      Cartesian coordinates.
#   --covalent_lig_atom_pattern arg    SMARTS expression for ligand atom that
#                                      will covalently bind protein.
#   --covalent_lig_atom_position arg   Optional.  Initial placement of covalently
#                                      bonding ligand atom in x,y,z Cartesian
#                                      coordinates.  If not specified,
#                                      OpenBabel's GetNewBondVector function will
#                                      be used to position ligand.
#   --covalent_fix_lig_atom_position   If covalent_lig_atom_position is
#                                      specified, fix the ligand atom to this
#                                      position as opposed to using this position
#                                      to define the initial structure.
#   --covalent_bond_order arg (=1)     Bond order of covalent bond. Default 1.
#   --covalent_optimize_lig            Optimize the covalent complex of ligand
#                                      and residue using UFF. This will change
#                                      bond angles and lengths of the ligand.
# 
# Scoring and minimization options:
#   --scoring arg                      specify alternative built-in scoring
#                                      function: ad4_scoring default dkoes_fast
#                                      dkoes_scoring dkoes_scoring_old vina
#                                      vinardo
#   --custom_scoring arg               custom scoring function file
#   --custom_atoms arg                 custom atom type parameters file
#   --score_only                       score provided ligand pose
#   --local_only                       local search only using autobox (you
#                                      probably want to use --minimize)
#   --minimize                         energy minimization
#   --randomize_only                   generate random poses, attempting to avoid
#                                      clashes
#   --num_mc_steps arg                 fixed number of monte carlo steps to take
#                                      in each chain
#   --max_mc_steps arg                 cap on number of monte carlo steps to take
#                                      in each chain
#   --num_mc_saved arg                 number of top poses saved in each monte
#                                      carlo chain
#   --temperature arg                  temperature for metropolis accept
#                                      criterion
#   --minimize_iters arg (=0)          number iterations of steepest descent;
#                                      default scales with rotors and usually
#                                      isn't sufficient for convergence
#   --accurate_line                    use accurate line search
#   --simple_ascent                    use simple gradient ascent
#   --minimize_early_term              Stop minimization before convergence
#                                      conditions are fully met.
#   --minimize_single_full             During docking perform a single full
#                                      minimization instead of a truncated
#                                      pre-evaluate followed by a full.
#   --approximation arg                approximation (linear, spline, or exact)
#                                      to use
#   --factor arg                       approximation factor: higher results in a
#                                      finer-grained approximation
#   --force_cap arg                    max allowed force; lower values more
#                                      gently minimize clashing structures
#   --user_grid arg                    Autodock map file for user grid data based
#                                      calculations
#   --user_grid_lambda arg (=-1)       Scales user_grid and functional scoring
#   --print_terms                      Print all available terms with default
#                                      parameterizations
#   --print_atom_types                 Print all available atom types
# 
# Convolutional neural net (CNN) scoring:
#   --cnn_scoring arg (=1)             Amount of CNN scoring: none, rescore
#                                      (default), refinement, metrorescore
#                                      (metropolis+rescore), metrorefine
#                                      (metropolis+refine), all
#   --cnn arg                          built-in model to use, specify
#                                      PREFIX_ensemble to evaluate an ensemble of
#                                      models starting with PREFIX:
#                                      all_default_to_default_1_3_1
#                                      all_default_to_default_1_3_2
#                                      all_default_to_default_1_3_3
#                                      crossdock_default2018
#                                      crossdock_default2018_1
#                                      crossdock_default2018_1_3
#                                      crossdock_default2018_1_3_1
#                                      crossdock_default2018_1_3_2
#                                      crossdock_default2018_1_3_3
#                                      crossdock_default2018_1_3_4
#                                      crossdock_default2018_2
#                                      crossdock_default2018_3
#                                      crossdock_default2018_4
#                                      crossdock_default2018_KD_1
#                                      crossdock_default2018_KD_2
#                                      crossdock_default2018_KD_3
#                                      crossdock_default2018_KD_4
#                                      crossdock_default2018_KD_5 default1.0
#                                      default2017 dense dense_1 dense_1_3
#                                      dense_1_3_1 dense_1_3_2 dense_1_3_3
#                                      dense_1_3_4 dense_1_3_PT_KD
#                                      dense_1_3_PT_KD_1 dense_1_3_PT_KD_2
#                                      dense_1_3_PT_KD_3 dense_1_3_PT_KD_4
#                                      dense_1_3_PT_KD_def2018
#                                      dense_1_3_PT_KD_def2018_1
#                                      dense_1_3_PT_KD_def2018_2
#                                      dense_1_3_PT_KD_def2018_3
#                                      dense_1_3_PT_KD_def2018_4 dense_2 dense_3
#                                      dense_4 fast general_default2018
#                                      general_default2018_1
#                                      general_default2018_2
#                                      general_default2018_3
#                                      general_default2018_4
#                                      general_default2018_KD_1
#                                      general_default2018_KD_2
#                                      general_default2018_KD_3
#                                      general_default2018_KD_4
#                                      general_default2018_KD_5
#                                      redock_default2018 redock_default2018_1
#                                      redock_default2018_1_3
#                                      redock_default2018_1_3_1
#                                      redock_default2018_1_3_2
#                                      redock_default2018_1_3_3
#                                      redock_default2018_1_3_4
#                                      redock_default2018_2 redock_default2018_3
#                                      redock_default2018_4 redock_default2018_KD
#                                      _1 redock_default2018_KD_2
#                                      redock_default2018_KD_3
#                                      redock_default2018_KD_4
#                                      redock_default2018_KD_5
#   --cnn_model arg                    torch cnn model file; if not specified a
#                                      default model ensemble will be used
#   --cnn_rotation arg (=0)            evaluate multiple rotations of pose (max
#                                      24)
#   --cnn_mix_emp_force                Merge CNN and empirical minus forces
#   --cnn_mix_emp_energy               Merge CNN and empirical energy
#   --cnn_empirical_weight arg (=1)    Weight for scaling and merging empirical
#                                      force and energy
#   --cnn_center_x arg                 X coordinate of the CNN center
#   --cnn_center_y arg                 Y coordinate of the CNN center
#   --cnn_center_z arg                 Z coordinate of the CNN center
#   --cnn_verbose                      Enable verbose output for CNN debugging
# 
# Output:
#   -o [ --out ] arg                   output file name, format taken from file
#                                      extension
#   --out_flex arg                     output file for flexible receptor residues
#   --log arg                          optionally, write log file
#   --atom_terms arg                   optionally write per-atom interaction term
#                                      values
#   --atom_term_data                   embedded per-atom interaction terms in
#                                      output sd data
#   --pose_sort_order arg (=0)         How to sort docking results: CNNscore
#                                      (default), CNNaffinity, Energy
#   --full_flex_output                 Output entire structure for out_flex, not
#                                      just flexible residues.
# 
# Misc (optional):
#   --cpu arg                          the number of CPUs to use (the default is
#                                      to try to detect the number of CPUs or,
#                                      failing that, use 1)
#   --seed arg                         explicit random seed
#   --exhaustiveness arg (=8)          exhaustiveness of the global search
#                                      (roughly proportional to time)
#   --num_modes arg (=9)               maximum number of binding modes to
#                                      generate
#   --min_rmsd_filter arg (=1)         rmsd value used to filter final poses to
#                                      remove redundancy
#   -q [ --quiet ]                     Suppress output messages
#   --addH arg                         automatically add hydrogens in ligands (on
#                                      by default)
#   --stripH arg                       remove polar hydrogens from molecule
#                                      _after_ performing atom typing for
#                                      efficiency (off by default - nonpolar are
#                                      always removed)
#   --device arg (=0)                  GPU device to use
#   --no_gpu                           Disable GPU acceleration, even if
#                                      available.
# 
# Configuration file (optional):
#   --config arg                       the above options can be put here
# 
# Information (optional):
#   --help                             display usage summary
#   --help_hidden                      display usage summary with hidden options
#   --version                          display program version
# 
# 
#   Execute the next cell to run gnina.```

# ## Functions (to be expanded)

# In[4]:


import useful_rdkit_utils as uru
from rdkit import Chem
from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.rdBase import BlockLogs
import pandas as pd

docking_results = ""
# ------------------------------------------------------------------------------
# 4.1 Docking modes
# ------------------------------------------------------------------------------

def single_ligand_docking():
    subprocess.run([
    "./gnina",
    "-r", f"{protein_directory}/{pdb_id}_A.pdbqt",
    "-l", f"{ligand_directory}/{ligand_id}.sdf",
    "--autobox_ligand", f"{ligand_directory}/{ligand_id}_corrected_pose.sdf",
    "-o", f"{docking_results_directory}/{ligand_id}_docked_{pdb_id}.sdf",
    "--seed", "0",
    "--exhaustiveness", f"{ex}",
    *gpu_flag,  # expands to nothing or ["--no_gpu"]
    *cnn_flag  # expands to nothing or ["--cnn_scoring=none"]
    ])

def batch_docking():
    subprocess.run([
    "./gnina",
    "-r", f"{protein_directory}/{pdb_id}_A.pdbqt",
    "-l", f"{ligand_directory}/ligands_to_dock.sdf",
    "--autobox_ligand", f"{ligand_directory}/{ligand_id}_corrected_pose.sdf",
    "-o", f"{docking_results_directory}/multiple_ligands_docked_{pdb_id}.sdf",
    "--seed", "0",
    "--exhaustiveness", f"{ex}",
    *gpu_flag,  # expands to nothing or ["--no_gpu"]
    *cnn_flag  # expands to nothing or ["--cnn_scoring=none"]
    ])

def flexible_docking():
    subprocess.run([
    "./gnina",
    "-r", f"{protein_directory}/{pdb_id}_A.pdbqt",
    "-l", f"{ligand_directory}/{ligand_id}.sdf",
    "--autobox_ligand", f"{ligand_directory}/{ligand_id}_corrected_pose.sdf",
    "-o", f"{docking_results_directory}/{ligand_id}_flex.sdf",
    "--flexdist_ligand", f"{ligand_directory}/{ligand_id}_corrected_pose.sdf",
    "--flexdist", "3.59",
    "--seed", "0",
    "--exhaustiveness 64"
    *gpu_flag,  # expands to nothing or ["--no_gpu"]
    *cnn_flag  # expands to nothing or ["--cnn_scoring=none"]
    ])

def unknown_site_docking():
    subprocess.run([
    "./gnina",
    "-r", f"{protein_directory}/{pdb_id}_A.pdbqt",
    "-l", f"{ligand_directory}/{ligand_id}.sdf",
    "--autobox_ligand", f"{protein_directory}/{pdb_id}_A.pdbqt",
    "-o", f"{docking_results_directory}/{ligand_id}_docked_whole_{pdb_id}.sdf",
    "--seed", "0",
    "--exhaustiveness", f"{ex}",
    *gpu_flag,  # expands to nothing or ["--no_gpu"]
    *cnn_flag  # expands to nothing or ["--cnn_scoring=none"]
    ])
    
# ------------------------------------------------------------------------------
# 4.2 Main functions
# ------------------------------------------------------------------------------

def docking_main():
    # Redocking with extracted ligand
    if selection == "a":
        print("\n === Single ligand docking ===")
        single_ligand_docking()
        rmsd_calculation()
        report(docking_results = f"{docking_results_directory}/{ligand_id}_docked_{pdb_id}.sdf")
    
    # Docking with multiple ligands
    elif selection == "b":
        print("\n === Batch docking ===")
        batch_docking()
        report(docking_results = f"{docking_results_directory}/multiple_ligands_docked_{pdb_id}.sdf")
    
    # Flexible docking
    elif selection == "c":
        print("Flexible docking")
        flexible_docking()
        report(docking_results = f"{docking_results_directory}/{ligand_id}_flex.sdf")
    
    # Docking on unknown site
    else:
        print("Docking on unknown site")
        unknown_site_docking()
        report(docking_results = f"{docking_results_directory}/{ligand_id}_docked_whole_{pdb_id}.sdf")
        
def report(docking_results):
    score_columns = [
        "minimizedAffinity",
        "CNNscore",
        "CNNaffinity",
        "CNN_VS",
        "CNNaffinity_variance",
    ]

    # Normalize input: accept either a string or a list of strings
    if isinstance(docking_results, str):
        sdf_paths = [docking_results]
    else:
        sdf_paths = docking_results
    
    df_list = []
    for filename in sdf_paths:
        with BlockLogs():
            df_list.append(PandasTools.LoadSDF(filename))
    
    combo_df = pd.concat(df_list)
    
    # Convert score columns to float if they exist and are not entirely empty
    for col in score_columns:
        if col in combo_df.columns:
            # If the column is completely empty, drop it
            if combo_df[col].isnull().all() or (combo_df[col] == "").all():
                combo_df.drop(columns=[col], inplace=True)
            else:
                combo_df[col] = combo_df[col].astype(float)
    
    # Extract SMILES from the RDKit molecule column
    combo_df["SMILES"] = combo_df["ROMol"].apply(lambda mol: Chem.MolToSmiles(mol) if mol is not None else None)
    
    # Save with SMILES included
    output_csv = f"{docking_results_directory}/docking_results_{pdb_id}.csv"
    combo_df.to_csv(output_csv, index=False)
    
    print(f"Docking results with SMILES saved to {output_csv}")
    
    combo_df.head()

# ------------------------------------------------------------------------------
# 4.3 Additional functions
# ------------------------------------------------------------------------------

def rmsd_calculation():
    print("\n === Running mcs (maximum common structure) rmsd calculation ===")
    cognate = Chem.MolFromMolFile(f"{ligand_directory}/{ligand_id}_corrected_pose.sdf")
    poses = Chem.SDMolSupplier(f"{docking_results_directory}/{ligand_id}_docked_{pdb_id}.sdf")
    for i, pose in enumerate(poses):
        RDLogger.DisableLog('rdApp.warning')
        n_match, rmsd = uru.mcs_rmsd(cognate, pose)
        print(f"{n_match}\t{rmsd:.2f}")


# In[ ]:


# Run gnina
import os

os.makedirs("molecular_docking/docking_results", exist_ok=True)

selection = input("\n === Welcome to gnina. Please select the docking mode: \
                    \n single docking (rmsd and cnn score will be calculated): a \
                    \n batch docking: b \
                    \n flexible docking (maximum exhaustiveness): c \
                    \n docking on unknown sites: d \
                ")

inp = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80]

ex = int(input("\n Define exhaustiveness (8, 16, 24, 32, 40, 48, 56, 64): "))
if ex in inp:
    print(f"You have selected  exhaustiveness level {ex}")
    #ex = int(input("\n Define exhaustiveness (8, 16, 24, 32, 40, 48, 56, 64): ")) 
    ask_gpu = input("\n Run with GPU? [y/n]: ").strip().lower()
    ask_cnn = input("\n Run with CNN? [y/n]: ").strip().lower()
    
    # Decide GPU flag
    gpu_flag = []
    if ask_gpu == 'y':
        gpu_flag = []  # run with GPU
    elif ask_gpu == 'n':
        gpu_flag = ["--no_gpu"]  # run without GPU
    else:
        print("\n Invalid input. Please try again with 'y' or 'n'.")
    
    # Decide CNN flag
    cnn_flag = []
    if ask_cnn == 'y':
        cnn_flag = []  # run with CNN
    elif ask_cnn == 'n':
        cnn_flag = ["--cnn_scoring=none"]  # run without CNN
    else:
        print("\n Invalid input. Please try again with 'y' or 'n'.")

    print(f"You have selected option {selection}.")

    if __name__ == "__main__":
        # pdb_id = "4OHU"
        pdb_id = os.getenv("PARAM_PDB_ID")
        # ligand_name = "2TK"
        ligand_name = os.getenv("PARAM_LIGAND_NAME")

        #with open("params.json") as f:
            #data = json.load(f)

        #ligand_id = data["ligand_id"]
        #print("Ligand ID parsed:", ligand_id)

        docking_main()
else:
    print("\n === Invalid input. Please try again with the given options only! ===")
    print("\n === Exiting gnina ===")
    ex = ''


# In[ ]:




