#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Protein Preparation Script
Converted from Protein_Preparation/Protein_Preparation.ipynb
"""

import os
import subprocess

import numpy as np
from Bio.PDB import PDBParser
from openmm.app import PDBFile
from pdbfixer import PDBFixer


def main():
    """Main execution function"""
    pdb_id = os.getenv("PARAM_PDB_ID")
    pdb_file = f"./{pdb_id}.pdb"
    output_pdb_file = f"{os.path.splitext(pdb_file)[0]}_A_NAD.pdb"

    # Calculate ligand center (2TK coordinates)
    center_coords = calculate_ligand_center(pdb_file)  # use original for 2TK coords

    # Generate docking config
    generate_docking_config(center_coords)

    # Fix structure with PDBFixer (keep NAD)
    fixed_pdb_file = fix_pdb_structure(output_pdb_file)

    # Reattach NAD if PDBFixer removed it
    fixed_pdb_file = reattach_nad(output_pdb_file, fixed_pdb_file)

    # Add AMBER charges
    pqr_file = add_amber_charges(fixed_pdb_file)

    # Copy final receptor to results/
    copy_to_visualization(fixed_pdb_file)

    print("\n✅ Protein preparation completed successfully!")


# ------------------------------------------------------------------------------
# 3. Compute ligand box (from 2TK in original file)
# ------------------------------------------------------------------------------


def calculate_ligand_center(pdb_file):
    """Compute ligand (2TK) center coordinates"""
    print("\n=== Calculating ligand (2TK) center coordinates ===")

    ligand_name = "2TK"
    extracted_ligand_coords = []

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target_protein", pdb_file)

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == ligand_name:
                    coords = [atom.get_coord() for atom in residue]
                    if coords:
                        extracted_ligand_coords.append(coords)

    if not extracted_ligand_coords:
        print(f"❌ Ligand '{ligand_name}' not found in {pdb_file}")
        exit()

    coords_array = np.array(extracted_ligand_coords[0])
    center = np.mean(coords_array, axis=0)

    print(
        f"Center of {ligand_name}: ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})"
    )
    return center, coords_array


# ------------------------------------------------------------------------------
# 4. Generate config
# ------------------------------------------------------------------------------


def generate_docking_config(center_data):
    """Generate docking config file"""
    print("\n=== Generating docking configuration file ===")

    center, coords_array = center_data
    lig_min = coords_array.min(axis=0)
    lig_max = coords_array.max(axis=0)
    extent = lig_max - lig_min
    padding = 8.0
    min_size = 20.0
    size_vec = np.maximum(extent + padding, min_size)

    exhaustiveness = os.getenv("PARAM_EXHAUSTIVENESS")
    num_modes = os.getenv("PARAM_NUM_MODES")
    energy_range = os.getenv("PARAM_ENERGY_RANGE")

    config_path = "config.txt"
    config_lines = [
        f"center_x = {center[0]:.3f}",
        f"center_y = {center[1]:.3f}",
        f"center_z = {center[2]:.3f}",
        f"size_x   = {size_vec[0]:.1f}",
        f"size_y   = {size_vec[1]:.1f}",
        f"size_z   = {size_vec[2]:.1f}",
        f"exhaustiveness = {exhaustiveness}",
        f"num_modes = {num_modes}",
        f"energy_range = {energy_range}",
    ]

    with open(config_path, "w") as f:
        f.write("\n".join(config_lines) + "\n")

    print(f"✅ Wrote docking config: {config_path}")


# ------------------------------------------------------------------------------
# 5. Fix structure (keep NAD)
# ------------------------------------------------------------------------------


def fix_pdb_structure(output_pdb_file):
    """Fix structure using PDBFixer, keeping NAD"""
    print(f"\n=== Fixing {output_pdb_file} with PDBFixer ===")

    fixed_pdb_file = f"{os.path.splitext(output_pdb_file)[0]}_fixed.pdb"
    fixer = PDBFixer(filename=output_pdb_file)

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    # Keep cofactors (remove only water)
    fixer.removeHeterogens(keepWater=False)
    fixer.addMissingHydrogens()

    with open(fixed_pdb_file, "w") as fout:
        PDBFile.writeFile(fixer.topology, fixer.positions, fout)

    print(f"✅ Fixed structure saved: {fixed_pdb_file}")
    return fixed_pdb_file


# ------------------------------------------------------------------------------
# 6. Reattach NAD if missing
# ------------------------------------------------------------------------------


def reattach_nad(original_pdb, fixed_pdb):
    """Reattach NAD residues from original PDB if missing after fixing"""
    output_pdb = f"{os.path.splitext(fixed_pdb)[0]}_with_NAD.pdb"

    with (
        open(original_pdb) as orig,
        open(fixed_pdb) as fixed,
        open(output_pdb, "w") as out,
    ):
        fixed_lines = fixed.readlines()
        nad_lines = [
            line
            for line in orig
            if "NAD" in line and line.startswith(("HETATM", "ATOM"))
        ]

        # Only append NAD lines if not already present
        if not any("NAD" in line for line in fixed_lines):
            fixed_lines += nad_lines
            print("🔁 NAD was missing — reattached to fixed structure.")
        else:
            print("✅ NAD already present, no reattachment needed.")

        out.writelines(fixed_lines)

    print(f"✅ Final receptor with NAD: {output_pdb}")
    return output_pdb


# ------------------------------------------------------------------------------
# 7. PQR generation
# ------------------------------------------------------------------------------


def add_amber_charges(fixed_pdb_file):
    """Add AMBER charges using PDB2PQR (optional)"""
    print(f"\n=== Adding AMBER charges to {fixed_pdb_file} ===")
    pqr_file = f"{os.path.splitext(fixed_pdb_file)[0]}.pqr"

    try:
        subprocess.run(
            ["pdb2pqr", "--ff=AMBER", fixed_pdb_file, pqr_file],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        print(f"✅ Generated PQR file: {pqr_file}")
        return pqr_file
    except subprocess.CalledProcessError as e:
        print("⚠️ PDB2PQR failed (nonstandard residue like NAD detected). Skipping.")
        print(f"Error details:\n{e.stderr.decode()}")
        return None
    except FileNotFoundError:
        print("⚠️ pdb2pqr not found in PATH — skipping charge calculation.")
        return None


# ------------------------------------------------------------------------------
# 8. Copy to results/
# ------------------------------------------------------------------------------


def copy_to_visualization(fixed_pdb_file):
    """Copy file to results directory"""
    print("\n=== Copying results ===")
    import shutil

    results_dir = "./results"
    os.makedirs(results_dir, exist_ok=True)
    shutil.copy(
        fixed_pdb_file, os.path.join(results_dir, os.path.basename(fixed_pdb_file))
    )
    print(f"✅ Copied to {results_dir}/")


# ------------------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
