#!/usr/bin/env python3
"""
Node 4: Ligand Preparation
Converts ligands from SDF to PDBQT format for docking
"""

import os
import sys
import json
import shutil
import subprocess
from pathlib import Path
from datetime import datetime

def convert_ligands():
    input_dir = Path("./")
    outputs_dir = input_dir / "outputs"
    outputs_dir.mkdir(exist_ok=True)
    
    print("[Ligand Preparation] ========== Flattening Inputs ==========")
    # Flatten inputs: Find all SDFs recursively and move to root
    for sdf in input_dir.glob("**/*.sdf"):
        if sdf.parent != input_dir:
            target = input_dir / sdf.name
            if not target.exists():
                shutil.copy(sdf, target)
                print(f"[Ligand Preparation] Flattened input: {sdf} → {target}")

    # Flatten PDBs (for pass-through)
    for pdb in input_dir.glob("**/*.pdb"):
        if pdb.parent != input_dir:
            target = input_dir / pdb.name
            if not target.exists():
                shutil.copy(pdb, target)
                print(f"[Ligand Preparation] Flattened input: {pdb} → {target}")

    # Flatten PDBQTs (for pass-through)
    for pdbqt in input_dir.glob("**/*.pdbqt"):
        if pdbqt.parent != input_dir:
            target = input_dir / pdbqt.name
            if not target.exists():
                shutil.copy(pdbqt, target)
                print(f"[Ligand Preparation] Flattened input: {pdbqt} → {target}")

    ligands = [f for f in input_dir.iterdir() if f.suffix.lower() == ".sdf"]
    print(f"[Ligand Preparation] Found {len(ligands)} ligands (.sdf)")

    if not ligands:
        print("[Ligand Preparation] ERROR: No SDF ligand files found.")
        sys.exit(1)

    print("[Ligand Preparation] ========== Converting Ligands ==========")
    ligands_data = []
    
    for lig in ligands:
        lig_name = lig.stem
        pdbqt_out = input_dir / f"{lig_name}.pdbqt"

        print(f" → Converting {lig.name} to PDBQT...", end=" ", flush=True)
        
        try:
            # Use OpenBabel to convert SDF to PDBQT, preserving all hydrogens
            cmd = [
                "obabel",
                str(lig),
                "-O", str(pdbqt_out),
                "-xh"  # Preserve input hydrogens
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print("Done")
            
            # Count atoms in PDBQT for validation
            with open(pdbqt_out) as f:
                pdbqt_lines = f.readlines()
            atom_count = sum(1 for line in pdbqt_lines if line.startswith("ATOM") or line.startswith("HETATM"))
            
            # Get rotatable bonds count from PDBQT
            torsion_count = 0
            for line in pdbqt_lines:
                if line.startswith("TORSDOF"):
                    torsion_count = int(line.split()[1])
                    break
            
            ligands_data.append({
                "ligand_name": lig_name,
                "input_sdf": lig.name,
                "output_pdbqt": pdbqt_out.name,
                "conversion_status": "success",
                "num_atoms": atom_count,
                "num_rotatable_bonds": torsion_count,
                "file_sizes": {
                    "sdf_bytes": lig.stat().st_size,
                    "pdbqt_bytes": pdbqt_out.stat().st_size
                },
                "note": "All explicit hydrogens preserved from SDF input"
            })
            
        except Exception as e:
            print(f"Failed: {e}")
            ligands_data.append({
                "ligand_name": lig_name,
                "input_sdf": lig.name,
                "output_pdbqt": pdbqt_out.name,
                "conversion_status": "failed",
                "error": str(e)
            })
    
    # Copy PDBQT files to outputs for visualization
    print("[Ligand Preparation] ========== Consolidating Outputs ==========")
    for ligand in ligands_data:
        if ligand["conversion_status"] == "success":
            pdbqt_src = input_dir / ligand["output_pdbqt"]
            pdbqt_dst = outputs_dir / ligand["output_pdbqt"]
            if pdbqt_src.exists():
                shutil.copy(pdbqt_src, pdbqt_dst)
                print(f"[Ligand Preparation] Copied {ligand['output_pdbqt']} to outputs/")
    
    # Pass-through: Ensure receptor files reach downstream nodes
    for pdb in input_dir.glob("*.pdb"):
        try:
            shutil.copy(pdb, outputs_dir / pdb.name)
            print(f"[Ligand Preparation] ✓ Passed {pdb.name} forward")
        except shutil.SameFileError:
            pass
        except Exception as e:
            print(f"[Ligand Preparation] Warning: {e}")
            
    for pdbqt in input_dir.glob("*.pdbqt"):
        if not (outputs_dir / pdbqt.name).exists():
            try:
                shutil.copy(pdbqt, outputs_dir / pdbqt.name)
                print(f"[Ligand Preparation] ✓ Passed {pdbqt.name} forward")
            except Exception as e:
                print(f"[Ligand Preparation] Warning: {e}")

    print(f"[Ligand Preparation] ✓ Consolidated files to {outputs_dir}")
    
    # Generate JSON metadata
    data = {
        "ligands": ligands_data,
        "total_count": len(ligands_data),
        "successful_conversions": sum(1 for lig in ligands_data if lig["conversion_status"] == "success"),
        "failed_conversions": sum(1 for lig in ligands_data if lig["conversion_status"] == "failed"),
        "conversion_method": "OpenBabel with -xh flag (preserves all explicit hydrogens)",
        "note": "PDBQT files preserve all atoms from SDF input for chemical accuracy",
        "timestamp": datetime.now().isoformat()
    }
    
    json_file = outputs_dir / "refined_ligands_metadata.json"
    with open(json_file, "w") as f:
        json.dump(data, f, indent=2)

    print(f"[Ligand Preparation] Metadata saved to {json_file}")
    print("[Ligand Preparation] ✅ Ligand preparation complete")


if __name__ == "__main__":
    convert_ligands()