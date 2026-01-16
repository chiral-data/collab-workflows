#!/usr/bin/env python3
"""
Node 5: Pocket Discovery
Predicts binding pockets using P2Rank and defines grid box
"""

import subprocess
import json
import csv
import os
import shutil
from pathlib import Path
from datetime import datetime
import sys

def run_p2rank():
    workdir = Path("./")
    outputs_dir = workdir / "outputs"
    outputs_dir.mkdir(exist_ok=True)
    
    print(f"[Pocket Discovery] Working directory: {workdir.resolve()}")
    
    print("[Pocket Discovery] ========== Flattening Inputs ==========")
    # Flatten inputs: Find all PDBs recursively
    for pdb in workdir.glob("**/*.pdb"):
        if pdb.parent != workdir:
            target = workdir / pdb.name
            if not target.exists():
                shutil.copy(pdb, target)
                print(f"[Pocket Discovery] Flattened input: {pdb} → {target}")

    # Flatten PDBQTs (for pass-through)
    for pdbqt in workdir.glob("**/*.pdbqt"):
        if pdbqt.parent != workdir:
            target = workdir / pdbqt.name
            if not target.exists():
                shutil.copy(pdbqt, target)
                print(f"[Pocket Discovery] Flattened input: {pdbqt} → {target}")

    print("[Pocket Discovery] ========== Finding Input Protein ==========")
    pdb_files = list(workdir.glob("*.pdb"))
    if not pdb_files:
        print("[Pocket Discovery] ERROR: No .pdb file found")
        sys.exit(1)

    protein = pdb_files[0]
    protein_abs = protein.resolve()
    print(f"[Pocket Discovery] Found protein file: {protein.name}")

    print("[Pocket Discovery] ========== Locating P2Rank ==========")
    prank_path = Path("/opt/p2rank/p2rank_2.5.1/prank")
    if not prank_path.exists():
        prank_path = Path("/opt/p2rank/prank")
        if not prank_path.exists():
            print(f"[Pocket Discovery] ERROR: prank executable not found")
            sys.exit(1)

    print("[Pocket Discovery] ========== Running P2Rank ==========")
    cmd = [
        str(prank_path), "predict",
        "-f", str(protein_abs),
        "-o", str(outputs_dir.resolve())
    ]

    print(f"[Pocket Discovery] Command: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[Pocket Discovery] ERROR: P2Rank failed: {e}")
        sys.exit(1)

    print("[Pocket Discovery] ========== Parsing Results ==========")
    # Identify P2Rank output files
    predictions_csv = outputs_dir / f"{protein.name}_predictions.csv"
    pockets_pdb = outputs_dir / f"{protein.name}_predictions.pdb"
    residues_csv = outputs_dir / f"{protein.name}_residues.csv"
    
    if not predictions_csv.exists():
        predictions_csv = outputs_dir / f"{protein.stem}_predictions.csv"
        pockets_pdb = outputs_dir / f"{protein.stem}_predictions.pdb"
        residues_csv = outputs_dir / f"{protein.stem}_residues.csv"

    pockets_data = []
    pocket_residues = {}
    
    if predictions_csv.exists():
        with open(predictions_csv, 'r') as f:
            reader = csv.DictReader(f)
            fieldnames = [field.strip() for field in (reader.fieldnames or [])]
            reader.fieldnames = fieldnames
            for row in reader:
                cleaned_row = {k.strip(): v.strip() for k, v in row.items()}
                pocket_data = {
                    "pocket_name": cleaned_row.get("name", ""),
                    "rank": int(cleaned_row.get("rank", 0)),
                    "score": float(cleaned_row.get("score", 0.0)),
                    "probability": float(cleaned_row.get("probability", 0.0)),
                    "residue_count": int(cleaned_row.get("residue_count", 0)),
                    "surface_atoms": int(cleaned_row.get("surf_atoms", 0)),
                    "center_x": float(cleaned_row.get("center_x", 0.0)) if "center_x" in cleaned_row else None,
                    "center_y": float(cleaned_row.get("center_y", 0.0)) if "center_y" in cleaned_row else None,
                    "center_z": float(cleaned_row.get("center_z", 0.0)) if "center_z" in cleaned_row else None
                }
                pockets_data.append(pocket_data)
        print(f"[Pocket Discovery] Successfully parsed {len(pockets_data)} pockets")
    else:
        print(f"[Pocket Discovery] WARNING: {predictions_csv} not found")

    if residues_csv.exists():
        with open(residues_csv, 'r') as f:
            reader = csv.DictReader(f)
            fieldnames = [field.strip() for field in (reader.fieldnames or [])]
            reader.fieldnames = fieldnames
            for row in reader:
                cleaned_row = {k.strip(): v.strip() for k, v in row.items()}
                pocket_name = cleaned_row.get("pocket_name", "")
                if pocket_name not in pocket_residues:
                    pocket_residues[pocket_name] = []
                pocket_residues[pocket_name].append({
                    "residue": cleaned_row.get("residue_name", ""),
                    "chain": cleaned_row.get("chain", ""),
                    "residue_number": cleaned_row.get("residue_number", "")
                })

    print("[Pocket Discovery] ========== Consolidating Outputs ==========")
    # Standardize names in outputs/
    if pockets_pdb.exists() and pockets_pdb != (outputs_dir / "pockets.pdb"):
        shutil.copy(pockets_pdb, outputs_dir / "pockets.pdb")
    
    shutil.copy(protein, outputs_dir / "protein.pdb")
    
    # Grid Box Definition (merged from old Node 6)
    print("[Pocket Discovery] ========== Defining Grid Box ==========")
    grid_params = {}
    if pockets_data:
        selected_pocket = next((p for p in pockets_data if p["rank"] == 1), pockets_data[0])
        print(f"[Pocket Discovery] ✓ Selected top pocket: {selected_pocket['pocket_name']}")
        print(f"[Pocket Discovery] ✓ Score: {selected_pocket['score']}")
        print(f"[Pocket Discovery] ✓ Center: ({selected_pocket['center_x']}, {selected_pocket['center_y']}, {selected_pocket['center_z']})")
        
        grid_params = {
            "center_x": selected_pocket["center_x"],
            "center_y": selected_pocket["center_y"],
            "center_z": selected_pocket["center_z"],
            "size_x": 20.0,
            "size_y": 20.0,
            "size_z": 20.0
        }
        
        grid_file = outputs_dir / "grid_config.json"
        with open(grid_file, "w") as f:
            json.dump(grid_params, f, indent=2)
        print(f"[Pocket Discovery] ✓ Generated grid_config.json")
    else:
        print("[Pocket Discovery] WARNING: No pockets found. Cannot define grid box.")
    
    # Create comprehensive metadata
    data = {
        "input_protein": protein.name,
        "pockets": pockets_data,
        "pocket_residues": pocket_residues,
        "total_pockets": len(pockets_data),
        "top_pocket": pockets_data[0] if pockets_data else None,
        "pockets_pdb": "pockets.pdb" if (outputs_dir / "pockets.pdb").exists() else None,
        "protein_pdb": "protein.pdb",
        "grid_selection": {
            "selected_pocket": pockets_data[0] if pockets_data else None,
            "grid_params": grid_params
        },
        "grid_config": grid_params,
        "timestamp": datetime.now().isoformat()
    }
    
    # Pass-through files
    print("[Pocket Discovery] ========== Pass-Through Files ==========")
    receptor_pdbqt = workdir / "receptor.pdbqt"
    if receptor_pdbqt.exists():
        shutil.copy(receptor_pdbqt, outputs_dir / "receptor.pdbqt")
        print("[Pocket Discovery] ✓ Passed receptor.pdbqt forward")
    
    for sdf in workdir.glob("*.sdf"):
        try:
            shutil.copy(sdf, outputs_dir / sdf.name)
            print(f"[Pocket Discovery] ✓ Passed {sdf.name} forward")
        except shutil.SameFileError:
            pass
            
    for pdbqt in workdir.glob("*.pdbqt"):
        try:
            if not (outputs_dir / pdbqt.name).exists():
                shutil.copy(pdbqt, outputs_dir / pdbqt.name)
                print(f"[Pocket Discovery] ✓ Passed {pdbqt.name} forward")
        except shutil.SameFileError:
            pass
    
    # Save metadata
    json_file = outputs_dir / "pocket_discovery_metadata.json"
    with open(json_file, "w") as f:
        json.dump(data, f, indent=2)
    
    print(f"[Pocket Discovery] ✓ Metadata saved to {json_file}")
    print(f"[Pocket Discovery] ✓ Consolidated all files to {outputs_dir}")
    print("[Pocket Discovery] ✅ Pocket prediction complete")

if __name__ == "__main__":
    run_p2rank()