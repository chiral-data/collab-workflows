#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 11: SMINA - In-silico Screening - In-silico screening
Input: config.txt, 5Y7J_amber.pdb (or 5Y7J_amber.pqr), Ligands_select.sdf
Output: Docking score (docking result files)
"""

import glob
import os
import subprocess

DEFAULT_PDB_ID = "5Y7J"
CONFIG_FILE = "config.txt"
LIGANDS_SDF = "Ligands_select.sdf"
OUTPUT_DIR = "docking_results"

def main():
    """Main execution function"""
    print("=== Node 11: SMINA - In-silico screening ===")
    
    # Get PDB ID
    pdb_id = os.environ.get("PDB_ID", DEFAULT_PDB_ID)
    
    # Check input files (original code used 5Y7J_AB_chains_fixed.pdb)
    # Following workflow image, prioritize Node 10 output (amber.pqr)
    receptor_file = f"./{pdb_id}_amber.pqr"  # Prioritize PQR file
    if not os.path.exists(receptor_file):
        receptor_file = f"./{pdb_id}_amber.pdb"  # Fallback to PDB file
    if not os.path.exists(receptor_file):
        receptor_file = f"./{pdb_id}_AB_chains_fixed.pdb"  # Original code filename
    if not os.path.exists(receptor_file):
        receptor_file = f"./{pdb_id}_clean.pdb"  # Fallback to cleaned PDB
    
    if not os.path.exists(receptor_file):
        print(f"❌ Error: Receptor file not found.")
        print("   Please run Node 9 or Node 10 first.")
        exit(1)
    
    if not os.path.exists(CONFIG_FILE):
        print(f"❌ Error: {CONFIG_FILE} not found.")
        print("   Please run Node 8 (node_08_ligand_center.py) first.")
        exit(1)
    
    if not os.path.exists(LIGANDS_SDF):
        print(f"❌ Error: {LIGANDS_SDF} not found.")
        print("   Please run Node 5 (node_05_ligand_selection.py) first.")
        exit(1)
    
    print(f"Receptor file: {receptor_file}")
    print(f"Config file: {CONFIG_FILE}")
    print(f"Ligand file: {LIGANDS_SDF}")
    
    # Check smina command
    check_smina()
    
    # Create docking results directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}")
    
    # Process ligand files
    # Original code processed individual SDF files directly, but following workflow image,
    # use Ligands_select.sdf
    # If Ligands_select.sdf doesn't exist, fallback to original behavior (individual files)
    smina_path = "smina"
    
    if os.path.exists(LIGANDS_SDF):
        # Extract individual molecules from Ligands_select.sdf and dock
        from rdkit import Chem
        
        supplier = Chem.SDMolSupplier(LIGANDS_SDF)
        ligands = [mol for mol in supplier if mol is not None]
        
        print(f"Number of ligands to dock: {len(ligands)}")
        
        for i, mol in enumerate(ligands, 1):
            if mol is None:
                continue
            
            # Save molecule to temporary SDF file
            temp_ligand_file = f"{OUTPUT_DIR}/temp_ligand_{i}.sdf"
            writer = Chem.SDWriter(temp_ligand_file)
            writer.write(mol)
            writer.close()
            
            # Output file names
            fname = f"ligand_{i}"
            out_sdf = f"{OUTPUT_DIR}/{fname}_docked.sdf"
            out_log = f"{OUTPUT_DIR}/{fname}_log.txt"
            
            print(f"\n[{i}/{len(ligands)}] Docking {fname}...")
            
            try:
                result = subprocess.run(
                    [
                        smina_path,
                        "-r",
                        receptor_file,
                        "-l",
                        temp_ligand_file,
                        "--config",
                        CONFIG_FILE,
                        "-o",
                        out_sdf,
                        "--log",
                        out_log,
                        "--scoring",
                        "vina",
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                    timeout=300,
                )
                
                print(f"✓ Docking complete for {fname}.")
                
                affinity = extract_affinity_from_log(out_log)
                if affinity is not None:
                    print(f"  Binding energy: {affinity:.2f} kcal/mol")
            
            except subprocess.TimeoutExpired:
                print(f"✗ Docking timeout for {fname}.")
            except subprocess.CalledProcessError as e:
                print(f"✗ Error occurred during docking for {fname}: {e}")
                if e.stderr:
                    print(f"  Error details: {e.stderr[:200]}...")
            except Exception as e:
                print(f"✗ Unexpected error occurred during docking for {fname}: {e}")
            finally:
                if os.path.exists(temp_ligand_file):
                    os.remove(temp_ligand_file)
    else:
        # Fallback: Original code behavior (process individual SDF files directly)
        import glob
        print(f"Warning: {LIGANDS_SDF} not found. Searching for individual SDF files.")
        
        # For testing: first 11 ligands
        ligands = glob.glob("./constructed_library/clean_drug108*.sdf")
        # For production: all ligands
        # ligands = glob.glob("./constructed_library/clean_drug*.sdf")
        
        if not ligands:
            print("❌ Error: Ligand files not found.")
            exit(1)
        
        print(f"Number of ligand files found: {len(ligands)}")
        
        for i, lig in enumerate(ligands, 1):
            fname = os.path.splitext(os.path.basename(lig))[0]
            out_sdf = f"{OUTPUT_DIR}/{fname}_docked.sdf"
            out_log = f"{OUTPUT_DIR}/{fname}_log.txt"
            
            print(f"\n[{i}/{len(ligands)}] Docking {fname}...")
            
            try:
                result = subprocess.run(
                    [
                        smina_path,
                        "-r",
                        receptor_file,
                        "-l",
                        lig,
                        "--config",
                        CONFIG_FILE,
                        "-o",
                        out_sdf,
                        "--log",
                        out_log,
                        "--scoring",
                        "vina",
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                    timeout=300,
                )
                
                print(f"✓ Docking complete for {fname}.")
                
                affinity = extract_affinity_from_log(out_log)
                if affinity is not None:
                    print(f"  Binding energy: {affinity:.2f} kcal/mol")
            
            except subprocess.TimeoutExpired:
                print(f"✗ Docking timeout for {fname}.")
            except subprocess.CalledProcessError as e:
                print(f"✗ Error occurred during docking for {fname}: {e}")
                if e.stderr:
                    print(f"  Error details: {e.stderr[:200]}...")
            except Exception as e:
                print(f"✗ Unexpected error occurred during docking for {fname}: {e}")
    
    print(f"\n✅ Docking screening complete.")
    print(f"Results saved in {OUTPUT_DIR}/ directory.")


def check_smina():
    """Check smina command"""
    smina_path = "smina"
    
    try:
        result = subprocess.run(
            [smina_path, "--help"], capture_output=True, text=True, timeout=10
        )
        print("smina command is available.")
    except subprocess.TimeoutExpired:
        print("smina command check timed out.")
    except FileNotFoundError:
        print(f"❌ Error: {smina_path} not found.")
        print("Please verify that smina is correctly installed.")
        exit(1)
    except Exception as e:
        print(f"Error occurred while checking smina command: {e}")


def extract_affinity_from_log(log_file):
    """Extract binding energy from log file"""
    try:
        with open(log_file, "r") as f:
            for line in f:
                if line.strip().startswith("1 "):  # Get mode1 line
                    parts = line.split()
                    if len(parts) > 1:
                        try:
                            affinity = float(parts[1])  # Affinity value (kcal/mol)
                            return affinity
                        except ValueError:
                            pass
    except Exception as e:
        print(f"  Log file read error: {e}")
    return None


if __name__ == "__main__":
    main()

