#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dock all ligands from ligand_library/ onto 2B35 receptor (with NAD cofactor)
"""

import os
import subprocess


def main():
    print("=== Starting docking on 2B35 receptor (with NAD) ===")

    # Check environment and receptor
    print(f"Current working directory: {os.getcwd()}")
    check_smina()
    pdb_id = os.getenv("PARAM_PDB_ID")
    receptor_pdbqt = prepare_receptor(pdb_id)

    # Prepare folders
    ligand_dir = "ligand_library"
    output_dir = "docking_results"
    os.makedirs(output_dir, exist_ok=True)

    if not os.path.isdir(ligand_dir):
        print(f"✗ Error: ligand folder '{ligand_dir}' not found.")
        exit()

    # Dock all ligands in ligand_library/
    ligands = [f for f in os.listdir(ligand_dir) if f.endswith(".sdf")]
    if not ligands:
        print(f"✗ No .sdf files found in {ligand_dir}/")
        exit()

    print(f"\nFound {len(ligands)} ligand(s) to dock:\n" + "\n".join(ligands))

    for ligand_file in ligands:
        ligand_path = os.path.join(ligand_dir, ligand_file)
        ligand_name = os.path.splitext(ligand_file)[0]
        run_docking(ligand_name, receptor_pdbqt, ligand_path, output_dir)

    print("\n=== Docking complete ===")
    print(f"Results saved in: {output_dir}/")


def check_smina():
    """Ensure smina command is available."""
    print("\n--- Checking smina installation ---")
    try:
        result = subprocess.run(
            ["smina", "--version"], capture_output=True, text=True, timeout=10
        )
        print("✔ smina found:")
        print(result.stdout.strip() or result.stderr.strip())
    except FileNotFoundError:
        print("✗ Error: smina not found in PATH. Please install Smina and try again.")
        exit()
    except Exception as e:
        print(f"✗ Unexpected error checking smina: {e}")
        exit()


def prepare_receptor(pdb_id):
    """Prepare receptor PDBQT (preserve NAD)."""
    print("\n--- Preparing receptor ---")

    receptor_pdb = f"./{pdb_id}_A_NAD_fixed_with_NAD.pdb"
    receptor_pdbqt = f"./{pdb_id}_A_NAD_fixed_with_NAD.pdbqt"

    if os.path.exists(receptor_pdbqt):
        print(f"✔ Using existing receptor file: {receptor_pdbqt}")
        return receptor_pdbqt

    cmd = [
        "prepare_receptor4.py",
        "-r",
        receptor_pdb,
        "-o",
        receptor_pdbqt,
        "-A",
        "hydrogens",
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f"Receptor prepared successfully: {receptor_pdbqt}")
    except FileNotFoundError:
        print("✗ prepare_receptor4.py not found. Install MGLTools/AutoDockTools.")
        exit()
    except subprocess.CalledProcessError as e:
        print(f"✗ Error preparing receptor:\n{e}")
        exit()

    return receptor_pdbqt


def run_docking(ligand_name, receptor_pdbqt, ligand_path, output_dir):
    """Dock a single ligand using Smina."""
    print(f"\n--- Docking {ligand_name} ---")

    config_file = "./config.txt"
    out_sdf = os.path.join(output_dir, f"{ligand_name}_docked.sdf")
    out_log = os.path.join(output_dir, f"{ligand_name}_docking.log")

    smina_cmd = [
        "smina",
        "-r",
        receptor_pdbqt,
        "-l",
        ligand_path,
        "--config",
        config_file,
        "-o",
        out_sdf,
        "--log",
        out_log,
        "--scoring",
        "vina",
        "--num_modes",
        "1",
    ]

    print("Command:")
    print(" ".join(smina_cmd))

    try:
        subprocess.run(
            smina_cmd, check=True, capture_output=True, text=True, timeout=600
        )
        print(f"✔ Docking complete: {ligand_name}")
    except subprocess.TimeoutExpired:
        print(f"⚠ Docking timed out for {ligand_name}")
    except subprocess.CalledProcessError as e:
        print(f"✗ Docking failed for {ligand_name}: {e}")
    except Exception as e:
        print(f"✗ Unexpected error during docking: {e}")

    affinity = extract_affinity(out_log)
    if affinity is not None:
        print(f"Predicted binding energy: {affinity:.2f} kcal/mol")
    else:
        print("⚠ Could not extract binding energy.")


def extract_affinity(log_file):
    """Extract top binding energy from a Smina log file."""
    try:
        with open(log_file, "r") as f:
            for line in f:
                if line.strip().startswith("1 "):
                    parts = line.split()
                    if len(parts) > 1:
                        return float(parts[1])
    except Exception as e:
        print(f"Could not read {log_file}: {e}")
    return None


if __name__ == "__main__":
    main()
