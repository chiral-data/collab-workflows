#!/usr/bin/env python
# coding: utf-8

# ## Protein preparation

# Please install the following packages:
# 
# - PDBFixer
# - Biopython
# - OpenMM (and OpenMMForceFields)
# - OpenBabel

# In[1]:


#pdb_id = input("Enter PDB code: ") # The Protein ID we're looking at
#pdb_id = os.getenv("PARAM_PDB_ID")
#pdb_id = "4OHU"

#protein_directory = os.getenv("PARAM_PROTEIN_DIRECTORY")


# In[ ]:


import os
import requests
from Bio.PDB import PDBParser, Select, PDBIO
from pdbfixer import PDBFixer
from openmm.app import PDBFile, ForceField, Simulation
from openmm import VerletIntegrator
import openmm.unit as unit

def select_chain(pdb_file):
    # From now on we will work with only one domain/chain of the target protein
    print("\n=== Selecting chain A if protein contains multiple chains ===")
    class ChainSelector(Select):
        def __init__(self, target_chain):
            self.target_chain = target_chain
        def accept_chain(self, chain):
            return chain.id == self.target_chain
    
    # Load structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    # Save each chain (monomer) as a separate PDB
    io = PDBIO()
    for model in structure:
        for chain in model:
            chain_id = chain.id
            io.set_structure(structure)
            io.save(f"{protein_directory}/{pdb_id}_{chain_id}.pdb", ChainSelector(chain_id))
            print(f"\n=== {pdb_id}_{chain_id}.pdb has been extracted.===")
    
    # Save chain A as a separate PDB file
    #io = PDBIO()
    #io.set_structure(structure)
    #io.save(f"{protein_directory}/{pdb_id}_A.pdb", ChainSelector("A"))
    
    print(f"\n=== By default chain A of {pdb_id} was selected for further processing!===")

def fixing_protein(pdb_file_chain_A):
    
    fixer = PDBFixer(filename = pdb_file_chain_A)
    
    print("\n === Starting PDBFixer...===")
    # Fixing the structure at pH 7.4
    fixer.findMissingResidues()
    fixer.missingResidues
    fixer.findNonstandardResidues()
    print(fixer.nonstandardResidues)
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    print(fixer.missingAtoms)
    print(fixer.missingTerminals)
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    print("\n === Loading force field (amber14-all.xml', 'amber14/tip3p.xml)... ===")
    # Load a force field (e.g., Amber)
    #forcefield = ForceField('amber19-all.xml', 'amber19/tip3pfb.xml')
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

    # Create OpenMM system for minimization
    #system = forcefield.createSystem(fixer.topology, ignoreExternalBonds=True)
    system = forcefield.createSystem(fixer.topology, ignoreExternalBonds=False)
    system.getForces()
    print("\n === Force field loaded ===")

    print("\n === Creating simulation for minimization ===")
    # Use a generic VerletIntegrator
    integrator = VerletIntegrator(0.001 * unit.picoseconds)
    # Optional, if you have access to a CUDA GPU, comment out the next line and uncomment the one after it
    platform = None
    # platform = Platform.getPlatformByName('CUDA')
    simulation = Simulation(fixer.topology, system, integrator, platform)
    simulation.context.setPositions(fixer.positions)
    print("\n Minimizing energy...")
    simulation.minimizeEnergy()
    minimized_positions = simulation.context.getState(getPositions=True).getPositions()
    
    # Write minimized structure to a PDB file
    with open(f"{protein_directory}/{pdb_id}_A_fixed.pdb", 'w') as output:
        PDBFile.writeFile(fixer.topology, minimized_positions, output)
    
    print(f"\n Minimization complete. Minimized structure saved to {protein_directory}/{pdb_id}_A_fixed_.pdb")

    #with open(f"{protein_directory}/{pdb_id}_A_fix_heavy.pdb", 'w') as f:
        #PDBFile.writeFile(fixer.topology, fixer.positions, f, True)

    print("\n === Generating pdbqt file ===")
    # Invoke OpenBabel's CLI from Python. Can also use subprocess as its safer, but os.system works fine here.
    receptor_pdbqt_path = f"{protein_directory}/{pdb_id}_A.pdbqt"
    receptor_fixed_path = f"{protein_directory}/{pdb_id}_A_fixed.pdb"
    
    # Generate the PDBQT file.
    # We could have "--partialcharge <method>" as a flag if we wanted to compute the partial charges, but this will just assume they are all "0"
    # The "-xh" flag preserves the hydrogens we worked so hard to get.
    os.system(f"obabel -ipdb {receptor_fixed_path} -opdbqt -O {receptor_pdbqt_path}")
    print(f" \n === {pdb_id}_A.pdbqt has been generated and saved ===")
    # If you get a status code "2" here, rerun it. You want status code 0
    
    print("\n === Fixing protein complete ===")

if __name__ == "__main__":
    pdb_id = os.getenv("PARAM_PDB_ID") 
    protein_directory = "molecular_docking/protein_files"
    protein_filename = f"{pdb_id}.pdb"
    select_chain(f"{protein_directory}/{pdb_id}.pdb")
    fixing_protein(f"{protein_directory}/{pdb_id}_A.pdb")

