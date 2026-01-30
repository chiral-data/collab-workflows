#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Protein Preparation Script
Converted from Protein_Preparation/Protein_Preparation.ipynb
"""

import os
import subprocess
import urllib.request

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select
from openmm.app import PDBFile
from pdbfixer import PDBFixer

# Read parameters from environment variables
pdb_id = os.environ.get("PARAM_PDB_ID", "5Y7J")
ligand_name = os.environ.get("PARAM_LIGAND_NAME", "8OL")


def main():
    """Main execution function"""
    print("Starting protein preparation...")
    print(f"PDB ID: {pdb_id}")
    print(f"Ligand name: {ligand_name}")

    # Download PDB file
    pdb_file = download_pdb_file()

    # Extract AB chains
    output_pdb_file = extract_ab_chains(pdb_file)

    # Calculate ligand center coordinates
    center_coords = calculate_ligand_center(output_pdb_file)

    # Generate docking config file
    generate_docking_config(center_coords)

    # Fix PDB structure with PDBFixer
    fixed_pdb_file = fix_pdb_structure(output_pdb_file)

    # Add AMBER charges with PDB2PQR
    pqr_file = add_amber_charges(fixed_pdb_file)

    # Copy files to results directory
    copy_to_visualization(fixed_pdb_file)

    print("Protein preparation completed.")


def download_pdb_file():
    """Download PDB file"""
    print("\n=== Downloading PDB file ===")

    pdb_file = f"./{pdb_id}.pdb"
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    print(f"Downloading PDB ID {pdb_id}...")
    try:
        urllib.request.urlretrieve(url, pdb_file)
        print(f"Download complete: {pdb_file}")
        return pdb_file
    except Exception as e:
        print(f"Error during download: {e}")
        exit()


def extract_ab_chains(pdb_file):
    """Check and extract AB chains"""
    print("\n=== Extracting AB chains ===")

    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    chains = [chain.id for model in structure for chain in model]
    unique_chains = sorted(list(set(chains)))

    print(f"Chains found in PDB file: {', '.join(unique_chains)}")

    # Check if both A and B chains exist
    if "A" in unique_chains and "B" in unique_chains:
        print("Chains A and B detected. Extracting these chains.")
        output_pdb_file = f"{os.path.splitext(pdb_file)[0]}_AB_chains.pdb"

        # Class to select AB chains and ligand
        class ABChainAndLigandSelect(Select):
            def accept_chain(self, chain):
                return chain.id in ["A", "B"]

            def accept_residue(self, residue):
                return (
                    residue.get_parent().id in ["A", "B"]
                    or residue.get_resname() == ligand_name
                )

        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb_file, ABChainAndLigandSelect())
        print(f"Extracted AB chains and ligand: {output_pdb_file}")

    else:
        print("Chains A and B not detected. Using original file.")
        output_pdb_file = pdb_file

    return output_pdb_file


def calculate_ligand_center(output_pdb_file):
    """Calculate ligand binding coordinates (center)"""
    print("\n=== Calculating ligand center coordinates ===")

    extracted_ligand_coords = []

    # Load extracted or original PDB file
    parser = PDBParser()
    target_structure = parser.get_structure("target_protein", output_pdb_file)

    for model in target_structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == ligand_name:
                    coords = [atom.get_coord() for atom in residue]
                    if coords:
                        extracted_ligand_coords.append(coords)

    # Use only the first ligand if multiple found
    if extracted_ligand_coords:
        coords_array = np.array(extracted_ligand_coords[0])
        center = np.mean(coords_array, axis=0)

        print("--- Ligand center coordinates ---")
        print(f"Ligand name: {ligand_name}")
        print(f"Center (x, y, z): {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}")

        return center, coords_array
    else:
        print(
            f"Error: Ligand '{ligand_name}' not found in file '{output_pdb_file}'."
        )
        exit()


def generate_docking_config(center_data):
    """Auto-calculate grid size and generate config file"""
    print("\n=== Generating docking config file ===")

    center, coords_array = center_data

    # Determine grid size from ligand extent
    lig_min = coords_array.min(axis=0)
    lig_max = coords_array.max(axis=0)
    extent = lig_max - lig_min  # molecular size in x,y,z
    padding = 8.0  # margin around ligand [A]
    min_size = 20.0  # practical minimum grid [A]
    size_vec = np.maximum(extent + padding, min_size)

    print(
        f"Recommended grid size (A): size_x={size_vec[0]:.1f}, size_y={size_vec[1]:.1f}, size_z={size_vec[2]:.1f}"
    )

    # Write config file
    config_path = "config.txt"

    config_lines = [
        f"center_x = {center[0]:.3f}",
        f"center_y = {center[1]:.3f}",
        f"center_z = {center[2]:.3f}",
        "size_x   = 15",
        "size_y   = 15",
        "size_z   = 15",
        "exhaustiveness = 8",
        "num_modes = 5",
        "energy_range = 4",
    ]

    with open(config_path, "w", encoding="utf-8") as f:
        f.write("\n".join(config_lines) + "\n")

    print(f"Docking config saved to '{config_path}'.")


def fix_pdb_structure(output_pdb_file):
    """Fix and cleanup structure with PDBFixer"""
    print(f"\n=== Processing {output_pdb_file} with PDBFixer...")

    fixed_pdb_file = f"{os.path.splitext(output_pdb_file)[0]}_fixed.pdb"
    fixer = PDBFixer(filename=output_pdb_file)

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.addMissingHydrogens()

    with open(fixed_pdb_file, "w") as fout:
        PDBFile.writeFile(fixer.topology, fixer.positions, fout)
    print(f"Fixed structure saved: {fixed_pdb_file}")

    return fixed_pdb_file


def add_amber_charges(fixed_pdb_file):
    """Add AMBER charges and generate PQR file with PDB2PQR"""
    print(
        f"\n=== Adding AMBER charges to {fixed_pdb_file} with PDB2PQR..."
    )

    pqr_file = f"{os.path.splitext(fixed_pdb_file)[0]}.pqr"

    try:
        subprocess.run(
            ["pdb2pqr", "--ff=AMBER", fixed_pdb_file, pqr_file],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        print(f"PQR file with AMBER charges created: {pqr_file}")
        return pqr_file
    except subprocess.CalledProcessError as e:
        print(f"Error running PDB2PQR: {e.stderr.decode()}")
        exit()
    except FileNotFoundError:
        print("Error: pdb2pqr not found in system path.")
        exit()


def copy_to_visualization(fixed_pdb_file):
    """Copy files to results directory"""
    print("\n=== Copying to results directory ===")

    import shutil

    results_dir = "./results"
    os.makedirs(results_dir, exist_ok=True)

    shutil.copy(fixed_pdb_file, f"{results_dir}/{os.path.basename(fixed_pdb_file)}")
    print(f"Copied {fixed_pdb_file} to {results_dir}.")


if __name__ == "__main__":
    main()
