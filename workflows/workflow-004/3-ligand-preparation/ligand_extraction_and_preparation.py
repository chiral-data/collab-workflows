#!/usr/bin/env python
# coding: utf-8

# # Ligand preparation

# Please install the following packages:
# * `Biopython`
# * `MDAnalysis`
# * `RDKit`
# * `OpenMM` (and `OpenMMForceFields`)
# * `OpenBabel`
# * `Scrubber` (package: "`molscrub`")
# * `py3Dmol`
# * `rdkit utils`
# * `subprocess`

# In[ ]:


#pdb_id = input("Enter PDB code: ") # The Protein ID we're looking at
#pdb_id = os.getenv("PARAM_PDB_ID")
#pdb_id = "4OHU"


# In[15]:


import os
import requests

### Set protein directory ###
protein_directory = "molecular_docking/protein_files"
protein_filename = f"{pdb_id}_A.pdb"
protein_filepath = os.path.join(protein_directory, protein_filename)

### Set ligand directory ###
ligand_directory = "molecular_docking/ligand_structures"
ligand_filename = ""
ligand_filepath = os.path.join(ligand_directory, ligand_filename)

# Actually make the directory, the exist_ok flag lets the command execute even if the folder already exists. It does NOT overwrite existing data.
os.makedirs(ligand_directory, exist_ok=True)

print("Ligand directory: The following folder has been created: ", ligand_filepath)


# In[ ]:


import numpy as np
import MDAnalysis as mda
from molscrub import Scrub
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import subprocess

# ------------------------------------------------------------------------------
# 3.1 Select and extract ligand from PDB
# ------------------------------------------------------------------------------

def select_ligand_from_pdb():
    # Load the original PDB
    u = mda.Universe(f"{protein_directory}/{pdb_id}_A.pdb")
    
    ligands = u.select_atoms("not protein and not water")
    i = 0 # index
    ligand_residue_names = ligands.residues.resnames
    
    # Loop through all ligands present and prints out their code
    unique_ligands = list(dict.fromkeys(ligand_residue_names))
    print("Ligands found:")
    for i, lig in enumerate(unique_ligands):
        print(i, lig)
    
    ligand_id = unique_ligands[int(input('Enter index: '))]
    print(f"\n === You have selected {ligand_id} as your ligand ===")

    single_ligand = u.select_atoms(f"resname {ligand_id}")
    single_ligand.write(f"{ligand_directory}/{ligand_id}_fromPDB.pdb")
    print(f"Ligand {ligand_id} extracted from original PDB!")
    
    return ligand_id

# ------------------------------------------------------------------------------
# 3.2 Download ideal ligand from RCSB
# ------------------------------------------------------------------------------

def download_ideal_ligand():

    ideal_ligand_filename = f"{ligand_id}_ideal.sdf"
    print(f"Downloading ligand {ligand_id}...")
    ligand_url = f"https://files.rcsb.org/ligands/download/{ideal_ligand_filename}"
    ligand_request = requests.get(ligand_url)
    ligand_request.raise_for_status() # Check for errors

    ideal_ligand = f"{ligand_directory}/{ideal_ligand_filename}"

    with open(ideal_ligand, "w") as f:
        f.write(ligand_request.text)
    print(f"Saved ligand to {ideal_ligand}")
    return ideal_ligand

# ------------------------------------------------------------------------------
# 3.3 Fix extracted ligand and align it with ideal ligand
# ------------------------------------------------------------------------------
    
def fix_and_align(ideal_mol, pose_mol):
    
    print("\n === Fixing and correcting the pose of extracted ligand ===")
    
    ### Disconnect any organometal ###
    rdMolStandardize.DisconnectOrganometallicsInPlace(pose_mol)

    ### Remove disconnected fragments ###
    fragmenter = rdMolStandardize.FragmentRemover()
    pose_mol_f = fragmenter.remove(pose_mol)

    ### Choose largest fragment ###
    chooser = rdMolStandardize.LargestFragmentChooser()
    pose_mol_lf = chooser.choose(pose_mol_f)

    ### Assign bond orders from the template to the pose molecule ###
    corrected_pose = AllChem.AssignBondOrdersFromTemplate(ideal_mol, pose_mol_lf)

    ### Add hydrogens ###
    corrected_pose_with_H = Chem.AddHs(corrected_pose, addCoords=True)

    ### Save the corrected pose to an SDF file ###
    ligand_corrected_pose_file = f"{ligand_directory}/{ligand_id}_corrected_pose.sdf"
    writer = Chem.SDWriter(ligand_corrected_pose_file)
    writer.write(corrected_pose_with_H)
    writer.close()
    print("Extracted ligand fixed and aligned!")

    return corrected_pose_with_H

# ------------------------------------------------------------------------------
# 3.4 Scrub extracted ligand and list of ligands
# ------------------------------------------------------------------------------

def scrubbing_ligands():
    
    ### Extracted ligand ###
    
    print("\n === Initiating extracted ligand preparation ===")
    print("\n === Initiating scrubber ===")
    cmd = f"""scrub.py {ligand_directory}/{ligand_id}_corrected_pose.sdf -o {ligand_directory}/{ligand_id}.sdf"""
    subprocess.run(cmd, shell=True)
    print("\n === Writing scrubbed sdf file ===")
    print(f"{ligand_id}.sdf is ready for docking!")
    
    ### List of ligands ###
    
    print("\n === Initiating multiple ligand preparation ===")
    print("\n === Loading SMILES from csv file ===")
    smiles_supplier = Chem.SmilesMolSupplier(f"{ligand_directory}/ligands_to_dock.csv", delimiter=",")
    mols = []
    
    for mol in smiles_supplier:
        print(Chem.MolToSmiles(mol))
        mols.append(mol)
    
    ligands_to_dock_dirty = f"{ligand_directory}/ligands_to_dock_dirty.sdf"
    writer = Chem.SDWriter(ligands_to_dock_dirty)
    
    for m in mols:
        writer.write(m)
    writer.close()
    
    print("\n === SMILES from csv file were saved as ligands_to_dock.sdf ===")
    print("\n === Initiating scrubber ===")
    cmd = f"""scrub.py {ligand_directory}/ligands_to_dock_dirty.sdf -o {ligand_directory}/ligands_to_dock.sdf"""
    subprocess.run(cmd, shell=True)
    
    print("ligands_to_dock.sdf is ready for docking!")

    ### Summary ###
    print(f"Summary: {ligand_id}.sdf, {ligand_id}_corrected_pose.sdf, {ligand_id}_tautomers.sdf and ligands_to_dock.sdf have been generated and are ready for docking!")

if __name__ == "__main__":
    pdb_id = os.getenv("PARAM_PDB_ID") 
    # ligand_id = "2TK"
    ligand_id = select_ligand_from_pdb()
    download_ideal_ligand()
    corrected_pose_with_H = fix_and_align(ideal_mol = Chem.MolFromMolFile(f"{ligand_directory}/{ligand_id}_ideal.sdf", removeHs=True), pose_mol = Chem.MolFromPDBFile(f"{ligand_directory}/{ligand_id}_fromPDB.pdb", removeHs=True))
    scrubbing_ligands()



# In[20]:


import json
#ligand_id = ligand   # dynamically detected ligand ID
#pdb_id = pdb      # dynamically detected PDB ID

with open("params.json", "w") as f:
    json.dump({"ligand_id": ligand_id, 
              "pdb_id": pdb_id
              }, f)

print(f"✅ Saved {ligand_id} and {pdb_id} to params.json")


# In[ ]:


# Voila! The ligands are ready for docking!

