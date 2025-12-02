#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reporting - Generate report
Input: Docking score (docking result files from smina_screening)
Output: Score report (score report)
"""

import glob
import json
import os
import re
import shutil
from pathlib import Path

DEFAULT_PDB_ID = "5Y7J"
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
# Output from this node should be in outputs/results/ directory (as specified in job.toml)
INPUT_DIR = os.path.join(SCRIPT_DIR, "input")
# Base output directory
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
# docking_results from docking/outputs/ (files are directly in outputs/, not in a subdirectory)
DOCKING_OUTPUT_DIR = os.path.join(SCRIPT_DIR, "..", "docking", "outputs")
DOCKING_RESULTS_DIR = DOCKING_OUTPUT_DIR
# Fallback to input/ for backward compatibility
if not os.path.exists(DOCKING_RESULTS_DIR):
    DOCKING_RESULTS_DIR = INPUT_DIR
# Protein PDB from docking/input/ (the receptor used for docking)
DOCKING_INPUT_DIR = os.path.join(SCRIPT_DIR, "..", "docking", "input")

def main():
    """Main execution function"""
    global OUTPUT_DIR
    
    print("=== Generate report ===")
    
    # Reset directories to ensure correctness
    OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
    
    # Find docking results from docking/outputs/ (files are directly in outputs/, not in a subdirectory)
    if not os.path.exists(DOCKING_RESULTS_DIR):
        print(f"❌ Error: docking outputs directory not found.")
        print(f"   Searched in: {DOCKING_OUTPUT_DIR}")
        print(f"   Searched in: {INPUT_DIR}")
        if os.path.exists(os.path.join(SCRIPT_DIR, "..", "docking", "outputs")):
            files = [f for f in os.listdir(os.path.join(SCRIPT_DIR, "..", "docking", "outputs")) if f.endswith(('.sdf', '.txt'))]
            if files:
                print(f"   Available files in docking/outputs: {files[:10]}...")
        if os.path.exists(INPUT_DIR):
            dirs = [d for d in os.listdir(INPUT_DIR) if os.path.isdir(os.path.join(INPUT_DIR, d))]
            if dirs:
                print(f"   Available directories in input: {dirs}")
        print("   Please ensure docking has completed.")
        exit(1)
    
    print(f"✓ Found docking results in: {DOCKING_RESULTS_DIR}")
    
    # Find protein PDB file from docking/input/ (the receptor used for docking)
    # Fallback to input/ directory (copied by run.sh)
    receptor_input_dir = DOCKING_INPUT_DIR
    if not os.path.exists(receptor_input_dir):
        print(f"⚠ Warning: docking/input/ directory not found, trying input/ directory")
        receptor_input_dir = INPUT_DIR
    
    # Ensure OUTPUT_DIR is not nested (double-check and fix if needed)
    OUTPUT_DIR = os.path.abspath(OUTPUT_DIR)
    
    # Remove nested outputs directory if it exists
    nested_outputs = os.path.join(OUTPUT_DIR, "outputs")
    if os.path.exists(nested_outputs) and os.path.isdir(nested_outputs):
        print(f"⚠ Warning: Found nested outputs directory: {nested_outputs}")
        print(f"   Removing it to avoid confusion...")
        shutil.rmtree(nested_outputs)
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Verify that OUTPUT_DIR was created
    if not os.path.exists(OUTPUT_DIR):
        print(f"❌ Error: Failed to create output directory: {OUTPUT_DIR}")
        exit(1)
    
    print(f"Input directory (docking_results): {DOCKING_RESULTS_DIR}")
    print(f"Protein PDB directory (docking/input): {receptor_input_dir}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"✓ Output directory created: {OUTPUT_DIR}")
    
    # Generate docking result ranking
    try:
        results = generate_docking_ranking(DOCKING_RESULTS_DIR, OUTPUT_DIR)
    except Exception as e:
        print(f"❌ Error in generate_docking_ranking: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
    
    # Create protein-ligand complexes for top 3 and true_ligand
    try:
        copy_top_compounds(results, DOCKING_RESULTS_DIR, receptor_input_dir, OUTPUT_DIR)
    except Exception as e:
        print(f"❌ Error in copy_top_compounds: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
    
    # Verify output directory has complex PDB files and ranking file
    if os.path.exists(OUTPUT_DIR):
        pdb_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('.pdb')]
        txt_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('.txt')]
        print(f"\n✅ Report generation complete. Output directory contains:")
        print(f"   - {len(pdb_files)} complex PDB file(s)")
        print(f"   - {len(txt_files)} text file(s)")
        if pdb_files:
            print(f"\n   Complex PDB files:")
            for pdb_file in sorted(pdb_files):
                print(f"     - {pdb_file}")
        if txt_files:
            print(f"\n   Text files:")
            for txt_file in sorted(txt_files):
                print(f"     - {txt_file}")
    else:
        print(f"⚠ Warning: Output directory does not exist: {OUTPUT_DIR}")


def parse_smina_log(log_file):
    """
    Extract mode1 affinity value from Smina log file.
    
    Args:
        log_file (str): Path to log file
    
    Returns:
        float or None: Binding energy value (kcal/mol) or None
    """
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
        print(f"Warning: Error occurred while reading {log_file}: {e}")
    return None  # mode1 data not found


def generate_docking_ranking(input_dir, output_dir):
    """Generate docking result ranking and save to file"""
    print("\n--- Generate docking result ranking ---")
    
    # Get docking log files
    log_files = glob.glob(os.path.join(input_dir, "*_log.txt"))
    
    if not log_files:
        error_msg = f"No log files found in {input_dir}."
        print(f"❌ Error: {error_msg}")
        print(f"   Searched pattern: {os.path.join(input_dir, '*_log.txt')}")
        if os.path.exists(input_dir):
            available_files = os.listdir(input_dir)
            print(f"   Available files: {available_files[:20]}...")  # Show first 20 files
        raise FileNotFoundError(error_msg)
    
    # Extract affinity from each log
    results = []
    true_ligand_result = None
    for log_file in log_files:
        affinity = parse_smina_log(log_file)
        if affinity is not None:
            compound_name = os.path.basename(log_file).replace("_log.txt", "")
            results.append((compound_name, affinity))
            # Check if this is true_ligand
            if compound_name == "true_ligand":
                true_ligand_result = (compound_name, affinity)
    
    if not results:
        error_msg = "No valid docking results found."
        print(f"❌ Error: {error_msg}")
        raise ValueError(error_msg)
    
    # Sort by affinity (binding energy) in ascending order (lower = stronger binding)
    results.sort(key=lambda x: x[1])
    
    # Find true_ligand's rank in the sorted results
    true_ligand_rank = None
    if true_ligand_result:
        for rank, (compound, affinity) in enumerate(results, 1):
            if compound == "true_ligand":
                true_ligand_rank = rank
                break
    
    # Write ranking to file
    ranking_file = os.path.join(output_dir, "docking_ranking.txt")
    os.makedirs(output_dir, exist_ok=True)
    
    with open(ranking_file, "w", encoding="utf-8") as f:
        f.write("Docking Result Ranking (sorted by binding strength)\n")
        f.write("=" * 60 + "\n")
        f.write(f"Total compounds: {len(results)}\n")
        f.write(f"Binding energy unit: kcal/mol (lower = stronger binding)\n")
        f.write("=" * 60 + "\n\n")
        
        for rank, (compound, affinity) in enumerate(results, 1):
            marker = " [TRUE LIGAND]" if compound == "true_ligand" else ""
            f.write(f"Rank {rank:4d}: {compound:30s}  Binding energy: {affinity:8.2f} kcal/mol{marker}\n")
        
        # Add summary for true_ligand if it exists
        if true_ligand_rank:
            f.write("\n" + "=" * 60 + "\n")
            f.write(f"True Ligand Summary:\n")
            f.write(f"  Name: true_ligand\n")
            f.write(f"  Rank: {true_ligand_rank} / {len(results)}\n")
            f.write(f"  Binding energy: {true_ligand_result[1]:.2f} kcal/mol\n")
            f.write("=" * 60 + "\n")
    
    # Display results
    print(f"Ranked docking results for {len(results)} compounds")
    print(f"✓ Ranking saved to: {ranking_file}")
    if true_ligand_rank:
        print(f"✓ True ligand found: Rank {true_ligand_rank} / {len(results)} (Binding energy: {true_ligand_result[1]:.2f} kcal/mol)")
    
    # Display top 10 results
    print("\n--- Top 10 compounds ---")
    for rank, (compound, affinity) in enumerate(results[:10], 1):
        marker = " [TRUE LIGAND]" if compound == "true_ligand" else ""
        print(f"Rank {rank}: Compound {compound}, Binding energy: {affinity:.2f} kcal/mol{marker}")
    
    return results


def extract_best_pose_from_sdf(sdf_file, output_sdf_file):
    """
    Extract the best pose (first molecule) from an SDF file that may contain multiple poses.
    The first molecule in SMINA output SDF is the best pose (mode 1, highest binding affinity).
    
    Args:
        sdf_file: Path to input SDF file (may contain multiple poses)
        output_sdf_file: Path to output SDF file (will contain only the first pose)
    
    Returns:
        bool: True if successful, False otherwise
    """
    if not os.path.exists(sdf_file):
        print(f"  ⚠ Warning: SDF file not found: {sdf_file}")
        return False
    
    try:
        with open(sdf_file, "r") as f:
            lines = f.readlines()
        
        if len(lines) < 4:
            print(f"  ⚠ Warning: SDF file too short: {sdf_file}")
            return False
        
        # Find the first molecule (ends with "$$$$" or end of file)
        first_molecule_end = len(lines)
        for i, line in enumerate(lines):
            if line.strip() == "$$$$":
                first_molecule_end = i + 1
                break
        
        # Write the first molecule to output file
        with open(output_sdf_file, "w") as out_f:
            out_f.writelines(lines[:first_molecule_end])
        
        return True
    
    except Exception as e:
        print(f"  ⚠ Warning: Failed to extract best pose from {sdf_file}: {e}")
        return False


def create_protein_ligand_complex(receptor_file, docked_sdf, complex_pdb_path, ligand_resname="LIG", ligand_chain_id="L"):
    """
    Create a simple protein-ligand complex PDB by combining:
    - receptor_file: cleaned protein PDB from docking/input/
    - docked_sdf: docked ligand SDF (best pose, first molecule)
    
    Notes:
    - This extracts the best pose (first molecule) from the SDF file.
    - It parses the first molecule in the SDF (V2000 style) and writes HETATM records.
    """
    if not os.path.exists(receptor_file):
        print(f"  ⚠ Warning: Receptor file not found for complex generation: {receptor_file}")
        return False
    if not os.path.exists(docked_sdf):
        print(f"  ⚠ Warning: Docked SDF not found for complex generation: {docked_sdf}")
        return False
    
    try:
        # Read receptor PDB
        with open(receptor_file, "r") as f:
            receptor_lines = f.readlines()
        
        # Determine the last atom serial number in receptor (for nicer PDB)
        max_serial = 0
        for line in receptor_lines:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                try:
                    serial = int(line[6:11])
                    if serial > max_serial:
                        max_serial = serial
                except ValueError:
                    continue
        
        # Parse SDF (first molecule only - best pose)
        with open(docked_sdf, "r") as f:
            sdf_lines = f.readlines()
        
        if len(sdf_lines) < 4:
            print(f"  ⚠ Warning: SDF file too short for complex generation: {docked_sdf}")
            return False
        
        # Find the end of the first molecule (best pose)
        first_molecule_end = len(sdf_lines)
        for i, line in enumerate(sdf_lines):
            if line.strip() == "$$$$":
                first_molecule_end = i
                break
        
        # Extract first molecule lines
        first_molecule_lines = sdf_lines[:first_molecule_end]
        
        if len(first_molecule_lines) < 4:
            print(f"  ⚠ Warning: First molecule in SDF too short: {docked_sdf}")
            return False
        
        counts_line = first_molecule_lines[3]
        try:
            # V2000 counts line: columns (0-3)=natoms, (3-6)=nbonds
            natoms = int(counts_line[0:3])
            nbonds = int(counts_line[3:6])
        except Exception as e:
            print(f"  ⚠ Warning: Could not parse atom/bond count from SDF counts line: {counts_line.strip()}")
            return False
        
        atom_lines = first_molecule_lines[4:4 + natoms]
        if len(atom_lines) < natoms:
            print(f"  ⚠ Warning: SDF atom block shorter than expected in {docked_sdf}")
            return False
        
        # Parse bond lines (after atom lines)
        bond_lines = first_molecule_lines[4 + natoms:4 + natoms + nbonds]
        
        ligand_pdb_lines = []
        ligand_connect_lines = []
        serial = max_serial
        res_seq = 1
        
        # Map SDF atom index (1-based) to PDB serial number
        atom_index_to_serial = {}
        
        for atom_idx, atom_line in enumerate(atom_lines, 1):
            # SDF atom line format (V2000):
            # x(10.4) y(10.4) z(10.4) atom_symbol(>30-34) ...
            try:
                x = float(atom_line[0:10])
                y = float(atom_line[10:20])
                z = float(atom_line[20:30])
                symbol = atom_line[31:34].strip()
                if not symbol:
                    symbol = "C"
            except Exception:
                # Skip malformed atoms
                continue
            
            serial += 1
            atom_index_to_serial[atom_idx] = serial
            element = symbol[0].upper()
            # Build a simple HETATM line
            # Columns follow standard PDB formatting
            pdb_line = (
                f"HETATM{serial:5d} "
                f"{symbol:<4s}"
                f"{ligand_resname:>3s} "
                f"{ligand_chain_id:1s}"
                f"{res_seq:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{1.00:6.2f}{0.00:6.2f}          "
                f"{element:>2s}"
            )
            ligand_pdb_lines.append(pdb_line + "\n")
        
        if not ligand_pdb_lines:
            print(f"  ⚠ Warning: No atoms parsed from SDF for complex generation: {docked_sdf}")
            return False
        
        # Parse bond information and create CONECT records
        for bond_line in bond_lines:
            try:
                # SDF bond line format (V2000):
                # atom1(3) atom2(3) bond_type(3) ...
                # atom indices are 1-based
                atom1_idx = int(bond_line[0:3])
                atom2_idx = int(bond_line[3:6])
                bond_type = int(bond_line[6:9])
                
                # Convert to PDB serial numbers
                if atom1_idx in atom_index_to_serial and atom2_idx in atom_index_to_serial:
                    serial1 = atom_index_to_serial[atom1_idx]
                    serial2 = atom_index_to_serial[atom2_idx]
                    
                    # Create CONECT record
                    # CONECT format: CONECT serial1 serial2 [serial3] [serial4] [serial5]
                    conect_line = f"CONECT{serial1:5d}{serial2:5d}\n"
                    ligand_connect_lines.append(conect_line)
            except (ValueError, IndexError):
                # Skip malformed bonds
                continue
        
        # Write combined complex PDB
        with open(complex_pdb_path, "w") as out_f:
            # Write receptor PDB (excluding END)
            for line in receptor_lines:
                if not line.startswith("END"):
                    out_f.write(line)
            
            # Write TER to mark end of receptor
            out_f.write("TER\n")
            
            # Write ligand HETATM records
            for line in ligand_pdb_lines:
                out_f.write(line)
            
            # Write CONECT records for ligand bonds
            for line in ligand_connect_lines:
                out_f.write(line)
            
            # Write END
            out_f.write("END\n")
        
        return True
    
    except Exception as e:
        print(f"  ⚠ Warning: Failed to create complex PDB for {docked_sdf}: {e}")
        return False


def copy_top_compounds(results, docking_output_dir, docking_input_dir, output_dir):
    """
    Create protein-ligand complexes for top 1, 2, 3 and true_ligand.
    
    Steps:
    1. Get protein PDB from docking/input/
    2. Get top 1, 2, 3 and true_ligand SDF files from docking/outputs/
    3. Extract best pose from each SDF and create complexes
    4. Output complexes to outputs/
    """
    print("\n--- Create protein-ligand complexes for top 3 and true_ligand ---")
    
    if not results:
        print("Warning: No results available for copying.")
        return
    
    # Load parameters from global_params.json or environment variables
    global_pdb_id = load_global_params()
    pdb_id = os.environ.get("PDB_ID") or os.environ.get("PARAM_PDB_ID") or global_pdb_id
    chain_id = os.environ.get("CHAIN_ID", DEFAULT_CHAIN_ID)
    
    # Step 1: Find receptor file from docking/input/
    receptor_file = None
    
    # Priority 1: Chain-specific cleaned PDB (e.g., 4OHU_chain_A_clean.pdb)
    if chain_id:
        receptor_file = os.path.join(docking_input_dir, f"{pdb_id}_chain_{chain_id}_clean.pdb")
        if os.path.exists(receptor_file):
            print(f"✓ Found receptor file: {os.path.basename(receptor_file)}")
    
    # Priority 2: General cleaned PDB (e.g., 4OHU_clean.pdb)
    if not receptor_file or not os.path.exists(receptor_file):
        receptor_file = os.path.join(docking_input_dir, f"{pdb_id}_clean.pdb")
        if os.path.exists(receptor_file):
            print(f"✓ Found receptor file: {os.path.basename(receptor_file)}")
    
    # Priority 3: Search for any PDB file with _clean in the name
    if not receptor_file or not os.path.exists(receptor_file):
        clean_files = glob.glob(os.path.join(docking_input_dir, f"{pdb_id}*_clean.pdb"))
        if clean_files:
            receptor_file = clean_files[0]
            print(f"✓ Found receptor file: {os.path.basename(receptor_file)}")
    
    # Priority 4: Search for any PDB file starting with pdb_id
    if not receptor_file or not os.path.exists(receptor_file):
        pdb_files = glob.glob(os.path.join(docking_input_dir, f"{pdb_id}*.pdb"))
        if pdb_files:
            # Prefer files with "clean" in the name
            clean_files = [f for f in pdb_files if "clean" in os.path.basename(f).lower()]
            if clean_files:
                receptor_file = clean_files[0]
            else:
                receptor_file = pdb_files[0]
            print(f"✓ Found receptor file: {os.path.basename(receptor_file)}")
    
    if not receptor_file or not os.path.exists(receptor_file):
        print(f"❌ Error: Receptor file not found in {docking_input_dir}")
        print(f"   Searched for: {pdb_id}_chain_{chain_id}_clean.pdb, {pdb_id}_clean.pdb, {pdb_id}*_clean.pdb")
        if os.path.exists(docking_input_dir):
            available_files = [f for f in os.listdir(docking_input_dir) if f.endswith('.pdb')]
            print(f"   Available PDB files: {available_files[:10]}...")
        return
    
    # Destination directory
    dst_dir = Path(output_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)
    
    complex_count = 0
    processed_ligands = set()  # Track which ligands have been processed
    
    # Step 2 & 3: Process top 1, 2, 3
    top_n = min(3, len(results))
    top_compounds = results[:top_n]
    
    for rank, (ligand_name, affinity) in enumerate(top_compounds, 1):
        print(f"\nRank {rank}: {ligand_name} (Binding energy: {affinity:.2f} kcal/mol)")
        
        # Get SDF file from docking/outputs/
        src_sdf = Path(docking_output_dir) / f"{ligand_name}_docked.sdf"
        if not src_sdf.exists():
            print(f"  ⚠ Warning: Docked SDF file not found: {src_sdf}")
            continue
        
        # Step 3: Extract best pose and create complex
        complex_pdb = dst_dir / f"top{rank}_{ligand_name}_complex.pdb"
        if create_protein_ligand_complex(str(receptor_file), str(src_sdf), str(complex_pdb)):
            print(f"  ✓ Created complex PDB: {complex_pdb.name}")
            complex_count += 1
            processed_ligands.add(ligand_name)
    
    # Step 2 & 3: Process true_ligand (always output, even if it's in top 3)
    true_ligand_name = "true_ligand"
    true_ligand_src = Path(docking_output_dir) / f"{true_ligand_name}_docked.sdf"
    
    if true_ligand_src.exists():
        # Find true_ligand's rank in the full results
        true_ligand_rank = None
        for rank, (compound, affinity) in enumerate(results, 1):
            if compound == true_ligand_name:
                true_ligand_rank = rank
                break
        
        if true_ligand_rank:
            print(f"\n--- True Ligand (Rank {true_ligand_rank}) ---")
            true_ligand_affinity = next((aff for name, aff in results if name == true_ligand_name), None)
            if true_ligand_affinity:
                print(f"Rank {true_ligand_rank}: {true_ligand_name} (Binding energy: {true_ligand_affinity:.2f} kcal/mol)")
            
            # Step 3: Extract best pose and create complex
            # If true_ligand is already in top 3, the file already exists, so we skip creating it again
            if true_ligand_name in processed_ligands:
                print(f"  ℹ Note: true_ligand is already in top 3, complex PDB already created")
            else:
                # Create true_ligand complex if it's not in top 3
                complex_pdb = dst_dir / f"top{true_ligand_rank}_{true_ligand_name}_complex.pdb"
                if create_protein_ligand_complex(str(receptor_file), str(true_ligand_src), str(complex_pdb)):
                    print(f"  ✓ Created complex PDB: {complex_pdb.name}")
                    complex_count += 1
    else:
        print(f"\n⚠ Warning: true_ligand SDF file not found: {true_ligand_src}")
    
    # Summary
    expected_files = []
    for rank in range(1, min(4, len(results) + 1)):
        if rank <= len(top_compounds):
            ligand_name = top_compounds[rank - 1][0]
            expected_files.append(f"top{rank}_{ligand_name}_complex.pdb")
    
    # Add true_ligand if it's not in top 3
    if true_ligand_src.exists() and true_ligand_name not in processed_ligands:
        if true_ligand_rank:
            expected_files.append(f"top{true_ligand_rank}_{true_ligand_name}_complex.pdb")
    
    print(f"\n✅ Created {complex_count} protein-ligand complex PDB file(s) in {output_dir}/")
    print(f"   Expected output files:")
    for expected_file in expected_files:
        print(f"     - {expected_file}")




if __name__ == "__main__":
    main()

