#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SMINA - In-silico Screening - In-silico screening
Input: config.txt, cleaned PDB file, selected_compounds/ (individual SDF files)
Output: Docking score (docking result files)
"""

import glob
import json
import os
import subprocess

DEFAULT_PDB_ID = "5Y7J"
# Default chain ID (should match protein extraction)
DEFAULT_CHAIN_ID = "A"

def load_global_params():
    """Load parameters from global_params.json"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Try workflow-003/global_params.json first, then workflow-003/global_params.json
    global_params_file = os.path.join(script_dir, "..", "global_params.json")
    if not os.path.exists(global_params_file):
        global_params_file = os.path.join(script_dir, "..", "global_params.json")
    if os.path.exists(global_params_file):
        try:
            with open(global_params_file, "r") as f:
                params = json.load(f)
                return params.get("pdb_id", DEFAULT_PDB_ID)
        except Exception as e:
            print(f"⚠ Warning: Could not load global_params.json: {e}")
    return DEFAULT_PDB_ID
# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
# Output from this node should be in outputs/ directory
INPUT_DIR = os.path.join(SCRIPT_DIR, "input")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
# Output files go directly to outputs/ directory (not in a subdirectory)
DOCKING_RESULTS_DIR = OUTPUT_DIR

# Get files from three sources (all from input/ directory, copied by run.sh):
# 1. Protein from input/ (copied from protein_extraction/outputs/ by run.sh)
# 2. Config from input/ (copied from docking_setup/outputs/ by run.sh)
CONFIG_FILE = os.path.join(INPUT_DIR, "config.txt")
# Fallback to direct path from docking_setup/outputs/
if not os.path.exists(CONFIG_FILE):
    DOCKING_SETUP_DIR = os.path.join(SCRIPT_DIR, "..", "docking_setup", "outputs")
    CONFIG_FILE = os.path.join(DOCKING_SETUP_DIR, "config.txt")
# 3. Selected compounds from input/selected_compounds/ (copied from library_construction/outputs/selected_compounds/ by run.sh)
SELECTED_COMPOUNDS_DIR = os.path.join(INPUT_DIR, "selected_compounds")
# No fallback - must be in input/selected_compounds/ (copied by run.sh from library_construction/outputs/selected_compounds/)
# Optional: real_ligand.sdf from input/ (copied from docking_setup/outputs/ by run.sh)
REAL_LIGAND_INPUT = os.path.join(INPUT_DIR, "real_ligand.sdf")
# Protein extraction directory for fallback
PROTEIN_EXTRACTION_DIR = os.path.join(SCRIPT_DIR, "..", "protein_extraction", "outputs")

def main():
    """Main execution function"""
    print("=== SMINA - In-silico screening ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load parameters from global_params.json or environment variables
    global_pdb_id = load_global_params()
    pdb_id = os.environ.get("PDB_ID") or os.environ.get("PARAM_PDB_ID") or global_pdb_id
    chain_id = os.environ.get("CHAIN_ID", DEFAULT_CHAIN_ID)
    
    # Find receptor file (cleaned PDB) in input/ directory (copied from protein_extraction/outputs/ by run.sh)
    receptor_file = None
    pdb_patterns = []
    
    if chain_id:
        # Try chain-specific file first (e.g., 4OHU_chain_A_clean.pdb)
        pdb_patterns.append(f"{pdb_id}_chain_{chain_id}_clean.pdb")
        # Try general cleaned PDB
        pdb_patterns.append(f"{pdb_id}_clean.pdb")
    else:
        # Try general cleaned PDB
        pdb_patterns.append(f"{pdb_id}_clean.pdb")
    
    # Search for any matching PDB file in input/ directory first
    for pattern in pdb_patterns:
        candidate = os.path.join(INPUT_DIR, pattern)
        if os.path.exists(candidate):
            receptor_file = candidate
            break
    
    # If still not found, search for any PDB file containing pdb_id in input/ directory
    if not receptor_file:
        if os.path.exists(INPUT_DIR):
            all_pdb_files = glob.glob(os.path.join(INPUT_DIR, "*.pdb"))
            for pdb_file in all_pdb_files:
                basename = os.path.basename(pdb_file)
                if pdb_id in basename and ("clean" in basename or "fixed" in basename):
                    receptor_file = pdb_file
                    break
    
    # Fallback to protein_extraction/outputs/ directory
    if not receptor_file or not os.path.exists(receptor_file):
        for pattern in pdb_patterns:
            candidate = os.path.join(PROTEIN_EXTRACTION_DIR, pattern)
            if os.path.exists(candidate):
                receptor_file = candidate
                break
        
        # If still not found, search for any PDB file containing pdb_id in protein_extraction/outputs/
        if not receptor_file:
            if os.path.exists(PROTEIN_EXTRACTION_DIR):
                all_pdb_files = glob.glob(os.path.join(PROTEIN_EXTRACTION_DIR, "*.pdb"))
                for pdb_file in all_pdb_files:
                    basename = os.path.basename(pdb_file)
                    if pdb_id in basename and ("clean" in basename or "fixed" in basename):
                        receptor_file = pdb_file
                        break
    
    if not receptor_file or not os.path.exists(receptor_file):
        print(f"❌ Error: Receptor file not found.")
        print(f"   Searched in: {INPUT_DIR}")
        print(f"   Searched in: {PROTEIN_EXTRACTION_DIR}")
        if os.path.exists(INPUT_DIR):
            files = [f for f in os.listdir(INPUT_DIR) if f.endswith('.pdb')]
            if files:
                print(f"   Available PDB files in input: {files}")
        if os.path.exists(PROTEIN_EXTRACTION_DIR):
            files = [f for f in os.listdir(PROTEIN_EXTRACTION_DIR) if f.endswith('.pdb')]
            if files:
                print(f"   Available PDB files in protein_extraction/outputs: {files}")
        print(f"   Tried patterns: {pdb_patterns}")
        print("   Please ensure protein_extraction has completed and run.sh has copied files to input/.")
        exit(1)
    
    # Check config file from input/ (copied by run.sh)
    if not os.path.exists(CONFIG_FILE):
        print(f"❌ Error: config.txt not found.")
        print(f"   Searched in: {os.path.join(INPUT_DIR, 'config.txt')}")
        if os.path.exists(INPUT_DIR):
            files = [f for f in os.listdir(INPUT_DIR) if f.endswith('.txt')]
            if files:
                print(f"   Available text files in input: {files}")
        # Check fallback location
        DOCKING_SETUP_DIR = os.path.join(SCRIPT_DIR, "..", "docking_setup", "outputs")
        if os.path.exists(DOCKING_SETUP_DIR):
            files = [f for f in os.listdir(DOCKING_SETUP_DIR) if f.endswith('.txt')]
            if files:
                print(f"   Available text files in docking_setup/outputs: {files}")
        print("   Please ensure docking_setup has completed and run.sh has copied files to input/.")
        exit(1)
    
    # Check selected_compounds directory from input/selected_compounds/ (copied by run.sh from library_construction/outputs/selected_compounds/)
    if not os.path.exists(SELECTED_COMPOUNDS_DIR):
        print(f"❌ Error: selected_compounds directory not found.")
        print(f"   Searched in: {os.path.join(INPUT_DIR, 'selected_compounds')}")
        if os.path.exists(INPUT_DIR):
            dirs = [d for d in os.listdir(INPUT_DIR) if os.path.isdir(os.path.join(INPUT_DIR, d))]
            if dirs:
                print(f"   Available directories in input: {dirs}")
        # Check source location for debugging
        LIBRARY_CONSTRUCTION_DIR = os.path.join(SCRIPT_DIR, "..", "library_construction", "outputs", "selected_compounds")
        if os.path.exists(LIBRARY_CONSTRUCTION_DIR):
            print(f"   Source location exists: {LIBRARY_CONSTRUCTION_DIR}")
            print(f"   Please check if run.sh has copied files from library_construction/outputs/selected_compounds/ to input/selected_compounds/")
        else:
            print(f"   Source location does not exist: {LIBRARY_CONSTRUCTION_DIR}")
            print(f"   Please ensure library_construction has completed.")
        print("   Please ensure library_construction has completed and run.sh has copied files to input/.")
        exit(1)
    
    # Find all SDF files in selected_compounds directory (from input/selected_compounds/)
    sdf_pattern = os.path.join(SELECTED_COMPOUNDS_DIR, "*.sdf")
    ligand_files = sorted(glob.glob(sdf_pattern))
    
    # Check if real_ligand.sdf exists in selected_compounds directory (optional)
    # Note: real_ligand.sdf should already be in input/selected_compounds/ from library_construction
    real_ligand_in_selected = os.path.join(SELECTED_COMPOUNDS_DIR, "real_ligand.sdf")
    real_ligand_in_selected_exists = os.path.exists(real_ligand_in_selected)
    
    # Check if true_ligand.sdf exists in selected_compounds directory (optional)
    true_ligand_in_selected = os.path.join(SELECTED_COMPOUNDS_DIR, "true_ligand.sdf")
    true_ligand_in_selected_exists = os.path.exists(true_ligand_in_selected)
    
    # Note: We don't copy files here - all files should already be in input/selected_compounds/ from run.sh
    if not real_ligand_in_selected_exists and not true_ligand_in_selected_exists:
        print(f"⚠ Warning: real_ligand.sdf or true_ligand.sdf not found in {SELECTED_COMPOUNDS_DIR} (optional files)")
        print(f"   Searched in: {SELECTED_COMPOUNDS_DIR}")
    
    # Add real_ligand.sdf or true_ligand.sdf to ligand_files if it exists and is not already in the list
    if real_ligand_in_selected_exists:
        if real_ligand_in_selected not in ligand_files:
            ligand_files.append(real_ligand_in_selected)
            ligand_files = sorted(ligand_files)
            print(f"✓ Added real_ligand.sdf to docking list")
        else:
            print(f"✓ real_ligand.sdf is already in docking list")
    elif true_ligand_in_selected_exists:
        if true_ligand_in_selected not in ligand_files:
            ligand_files.append(true_ligand_in_selected)
            ligand_files = sorted(ligand_files)
            print(f"✓ Added true_ligand.sdf to docking list")
        else:
            print(f"✓ true_ligand.sdf is already in docking list")
    else:
        print(f"⚠ Warning: real_ligand.sdf or true_ligand.sdf will not be docked (file not found, optional)")
    
    if not ligand_files:
        print(f"❌ Error: No SDF files found in {SELECTED_COMPOUNDS_DIR}.")
        print("   Please run prepare_ligands and true_ligand_addition first.")
        exit(1)
    
    # Verify that real_ligand.sdf or true_ligand.sdf is in the list
    ligand_basenames = [os.path.basename(f) for f in ligand_files]
    if "real_ligand.sdf" in ligand_basenames or "true_ligand.sdf" in ligand_basenames:
        found_ligand = "real_ligand.sdf" if "real_ligand.sdf" in ligand_basenames else "true_ligand.sdf"
        print(f"✓ Verified: {found_ligand} is in the docking list (total: {len(ligand_files)} ligands)")
    else:
        print(f"⚠ Warning: real_ligand.sdf or true_ligand.sdf is not in the ligand files list after processing.")
        print(f"   Ligand files: {ligand_basenames[:10]}...")
    
    print(f"Receptor file: {receptor_file}")
    print(f"  Exists: {os.path.exists(receptor_file)}")
    if os.path.exists(receptor_file):
        print(f"  Size: {os.path.getsize(receptor_file) / 1024:.2f} KB")
    print(f"Config file: {CONFIG_FILE}")
    print(f"  Exists: {os.path.exists(CONFIG_FILE)}")
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, 'r') as f:
            print(f"  Content preview: {f.read()[:100]}...")
    print(f"Ligand directory: {SELECTED_COMPOUNDS_DIR}")
    print(f"  Exists: {os.path.exists(SELECTED_COMPOUNDS_DIR)}")
    print(f"Number of ligands to dock: {len(ligand_files)}")
    if ligand_files:
        print(f"  First ligand: {os.path.basename(ligand_files[0])}")
        print(f"    Exists: {os.path.exists(ligand_files[0])}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Docking results directory: {DOCKING_RESULTS_DIR}")
    
    # Check smina command
    check_smina()
    smina_path = get_smina_path()
    print(f"Using SMINA path: {smina_path}")
    
    # Process each ligand file
    successful_count = 0
    failed_count = 0
    
    for i, ligand_file in enumerate(ligand_files, 1):
        # Get base filename without extension for output naming
        base_name = os.path.splitext(os.path.basename(ligand_file))[0]
        out_sdf = os.path.join(DOCKING_RESULTS_DIR, f"{base_name}_docked.sdf")
        out_log = os.path.join(DOCKING_RESULTS_DIR, f"{base_name}_log.txt")
        
        print(f"\n[{i}/{len(ligand_files)}] Docking {base_name}...")
        print(f"  Ligand file: {ligand_file}")
        print(f"    Exists: {os.path.exists(ligand_file)}")
        print(f"  Output SDF: {out_sdf}")
        print(f"  Output log: {out_log}")
        
        # Build SMINA command
        # SMINA accepts PDB format for receptor, but we need to ensure paths are absolute
        abs_receptor = os.path.abspath(receptor_file)
        abs_ligand = os.path.abspath(ligand_file)
        abs_config = os.path.abspath(CONFIG_FILE)
        abs_out_sdf = os.path.abspath(out_sdf)
        abs_out_log = os.path.abspath(out_log)
        
        # Read config file to get center and size parameters
        config_params = {}
        try:
            with open(CONFIG_FILE, 'r') as f:
                for line in f:
                    line = line.strip()
                    if '=' in line:
                        key, value = line.split('=', 1)
                        key = key.strip()
                        value = value.strip()
                        config_params[key] = value
        except Exception as e:
            print(f"  ⚠ Warning: Could not read config file: {e}")
            print(f"  Will try to use --config option")
        
        # Get parameters from environment or config file or use defaults
        exhaustiveness = int(os.environ.get("PARAM_EXHAUSTIVENESS", config_params.get("exhaustiveness", "16")))
        num_modes = int(os.environ.get("PARAM_NUM_MODES", config_params.get("num_modes", "1")))
        energy_range = int(os.environ.get("PARAM_ENERGY_RANGE", config_params.get("energy_range", "3")))
        
        # Build SMINA command
        smina_cmd = [
            smina_path,
            "-r", abs_receptor,
            "-l", abs_ligand,
            "-o", abs_out_sdf,
            "--log", abs_out_log,
            "--scoring", "vina",
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--energy_range", str(energy_range),
        ]
        
        # Add center and size from config file if available
        if "center_x" in config_params and "center_y" in config_params and "center_z" in config_params:
            smina_cmd.extend([
                "--center_x", config_params["center_x"],
                "--center_y", config_params["center_y"],
                "--center_z", config_params["center_z"],
            ])
            print(f"  Using center: ({config_params['center_x']}, {config_params['center_y']}, {config_params['center_z']})")
        else:
            print(f"  ⚠ Warning: Center coordinates not found in config file")
        
        if "size_x" in config_params and "size_y" in config_params and "size_z" in config_params:
            smina_cmd.extend([
                "--size_x", config_params["size_x"],
                "--size_y", config_params["size_y"],
                "--size_z", config_params["size_z"],
            ])
            print(f"  Using size: ({config_params['size_x']}, {config_params['size_y']}, {config_params['size_z']})")
        else:
            print(f"  ⚠ Warning: Size parameters not found in config file, using defaults")
        
        print(f"  Command: {' '.join(smina_cmd)}")
        
        try:
            result = subprocess.run(
                smina_cmd,
                check=True,
                capture_output=True,
                text=True,
                timeout=300,
            )
            
            # Print SMINA output for debugging
            if result.stdout:
                stdout_preview = result.stdout[:500] if len(result.stdout) > 500 else result.stdout
                print(f"  SMINA stdout: {stdout_preview}")
            if result.stderr:
                stderr_preview = result.stderr[:500] if len(result.stderr) > 500 else result.stderr
                print(f"  SMINA stderr: {stderr_preview}")
            
            print(f"✓ Docking complete for {base_name}.")
            
            affinity = extract_affinity_from_log(out_log)
            if affinity is not None:
                print(f"  Binding energy: {affinity:.2f} kcal/mol")
            
            successful_count += 1
        
        except subprocess.TimeoutExpired:
            print(f"✗ Docking timeout for {base_name}.")
            failed_count += 1
        except subprocess.CalledProcessError as e:
            print(f"✗ Error occurred during docking for {base_name}: {e}")
            if e.stderr:
                print(f"  Error details: {e.stderr[:500]}...")
            if e.stdout:
                print(f"  Output: {e.stdout[:500]}...")
            failed_count += 1
        except Exception as e:
            print(f"✗ Unexpected error occurred during docking for {base_name}: {e}")
            failed_count += 1
    
    # Clean up any old docking_results subdirectory if it exists
    old_docking_results_dir = os.path.join(OUTPUT_DIR, "docking_results")
    if os.path.exists(old_docking_results_dir) and os.path.isdir(old_docking_results_dir):
        print(f"⚠ Removing old docking_results subdirectory: {old_docking_results_dir}")
        import shutil
        try:
            shutil.rmtree(old_docking_results_dir)
        except Exception as e:
            print(f"  Warning: Could not remove old directory: {e}")
    
    print(f"\n✅ Docking screening complete.")
    print(f"   Successful: {successful_count}/{len(ligand_files)}")
    if failed_count > 0:
        print(f"   Failed: {failed_count}/{len(ligand_files)}")
    print(f"Results saved in {OUTPUT_DIR}/ directory.")
    
    
def get_smina_path():
    """Get smina executable path"""
    # In Docker container (chiral.sakuracr.jp/smina:2025_11_06), smina is installed at /usr/local/bin/smina
    # Check if we're in a Docker container
    if os.path.exists("/.dockerenv") or os.path.exists("/usr/local/bin/smina"):
        # In Docker, use /usr/local/bin/smina
        if os.path.exists("/usr/local/bin/smina"):
            return "/usr/local/bin/smina"
        # Fallback to system smina
        return "smina"
    
    # For local testing on macOS, try local smina.osx.12 if it exists
    local_smina = os.path.join(SCRIPT_DIR, "smina.osx.12")
    if os.path.exists(local_smina):
        try:
            # Check if it's executable
            if os.access(local_smina, os.X_OK):
                return local_smina
        except:
            pass
    
    # Fallback to system smina command
    return "smina"


def check_smina():
    """Check smina command"""
    smina_path = get_smina_path()
    
    # Check if smina file exists
    if not os.path.exists(smina_path) and smina_path != "smina":
        print(f"⚠ Warning: {smina_path} does not exist, trying system smina...")
        smina_path = "smina"
    
    try:
        # Try to run smina with --version or just check if it's executable
        result = subprocess.run(
            [smina_path, "--version"], capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            print(f"✓ smina command is available at: {smina_path}")
            if result.stdout:
                print(f"  Version info: {result.stdout.strip()[:100]}")
        else:
            # Try --help as fallback
            result2 = subprocess.run(
                [smina_path, "--help"], capture_output=True, text=True, timeout=10
        )
            if result2.returncode == 0:
                print(f"✓ smina command is available at: {smina_path}")
            else:
                print(f"⚠ smina command returned non-zero exit code: {result.returncode}")
                if result.stderr:
                    print(f"  Error: {result.stderr[:200]}")
    except subprocess.TimeoutExpired:
        print("⚠ smina command check timed out.")
    except FileNotFoundError:
        print(f"❌ Error: {smina_path} not found.")
        print("Please verify that smina is correctly installed.")
        print("In Docker container, smina should be at /usr/local/bin/smina")
        exit(1)
    except Exception as e:
        print(f"⚠ Error occurred while checking smina command: {e}")
        print(f"  Will attempt to use smina anyway...")


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

