#!/usr/bin/env python3
"""
Node 3: Receptor Preparation
Prepares receptor structure for docking (PDB → PDBQT)
"""

import sys
import json
import subprocess
import shutil
from pathlib import Path
from datetime import datetime
from pdbfixer import PDBFixer
from openmm.app import PDBFile


def main():
    workdir = Path("./")
    outputs_dir = workdir / "outputs"
    outputs_dir.mkdir(exist_ok=True)

    print("[Receptor Preparation] ========== Flattening Inputs ==========")
    # Flatten inputs: Find all PDBs recursively
    for pdb in workdir.glob("**/*.pdb"):
        if pdb.parent != workdir:
            target = workdir / pdb.name
            if not target.exists():
                shutil.copy(pdb, target)
                print(f"[Receptor Preparation] Flattened input: {pdb} → {target}")

    # Flatten inputs: Find all SDFs recursively
    for sdf in workdir.glob("**/*.sdf"):
        if sdf.parent != workdir:
            target = workdir / sdf.name
            if not target.exists():
                shutil.copy(sdf, target)
                print(f"[Receptor Preparation] Flattened input: {sdf} → {target}")

    pdb_files = list(workdir.glob("*.pdb")) + list(workdir.glob("*.cif"))
    if not pdb_files:
        print("[Receptor Preparation] ERROR: No PDB or CIF file found")
        sys.exit(1)

    input_structure = pdb_files[0]
    fixed_pdb = workdir / "protein_fixed.pdb"
    receptor_pdbqt = workdir / "receptor.pdbqt"

    print(f"[Receptor Preparation] Input structure: {input_structure}")

    # PDBFixer - Clean structure
    print("[Receptor Preparation] ========== Running PDBFixer ==========")
    try:
        fixer = PDBFixer(filename=str(input_structure))

        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

        # Keep crystallographic waters
        fixer.removeHeterogens(False)

        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        # Add hydrogens for docking
        fixer.addMissingHydrogens(7.0)  # pH 7.0

        with open(fixed_pdb, "w") as f:
            PDBFile.writeFile(
                fixer.topology,
                fixer.positions,
                f,
                keepIds=True
            )

        print(f"[Receptor Preparation] Fixed PDB written: {fixed_pdb}")
    except Exception as e:
        print(f"[Receptor Preparation] ERROR in PDBFixer: {e}")
        sys.exit(1)

    # OpenBabel - Convert to PDBQT
    print("[Receptor Preparation] ========== Converting to PDBQT ==========")
    cmd = [
        "obabel",
        str(fixed_pdb),
        "-O", str(receptor_pdbqt),
        "-xr",  # Rigid molecule (for receptor)
        "-xh"   # Preserve explicit hydrogens
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"[Receptor Preparation] Receptor PDBQT saved: {receptor_pdbqt}")
    except subprocess.CalledProcessError as e:
        print(f"[Receptor Preparation] ERROR: OpenBabel conversion failed")
        print(e.stderr)
        sys.exit(1)

    # Generate JSON metadata
    print("[Receptor Preparation] ========== Generating Metadata ==========")
    
    # Copy files to outputs
    shutil.copy(fixed_pdb, outputs_dir / fixed_pdb.name)
    shutil.copy(receptor_pdbqt, outputs_dir / receptor_pdbqt.name)

    data = {
        "input_pdb": input_structure.name,
        "fixed_pdb": fixed_pdb.name,
        "output_pdbqt": receptor_pdbqt.name,
        "preparation_steps": [
            "Found and replaced non-standard residues",
            "Added missing atoms",
            "Added hydrogens at pH 7.0",
            "Converted to PDBQT format"
        ],
        "file_sizes": {
            "input_pdb_bytes": input_structure.stat().st_size,
            "fixed_pdb_bytes": fixed_pdb.stat().st_size,
            "pdbqt_bytes": receptor_pdbqt.stat().st_size
        },
        "conversion_status": "success",
        "timestamp": datetime.now().isoformat()
    }

    json_file = outputs_dir / "refined_receptor_metadata.json"
    with open(json_file, "w") as f:
        json.dump(data, f, indent=2)

    # Pass-through: Ensure original receptor.pdb and ligands reach downstream nodes
    for pdb in workdir.glob("*.pdb"):
        if pdb.name != "protein_fixed.pdb":  # Don't pass the intermediate fixer file
            try:
                shutil.copy(pdb, outputs_dir / pdb.name)
                print(f"[Receptor Preparation] ✓ Passed {pdb.name} forward")
            except shutil.SameFileError:
                pass
            except Exception as e:
                print(f"[Receptor Preparation] Warning: Failed to pass {pdb.name}: {e}")
                
    for sdf in workdir.glob("*.sdf"):
        try:
            shutil.copy(sdf, outputs_dir / sdf.name)
            print(f"[Receptor Preparation] ✓ Passed {sdf.name} forward")
        except shutil.SameFileError:
            pass
        except Exception as e:
            print(f"[Receptor Preparation] Warning: Failed to pass {sdf.name}: {e}")

    print(f"[Receptor Preparation] ✓ Consolidated metadata and all files to {outputs_dir}")
    print("[Receptor Preparation] ✅ Protein preparation complete")


if __name__ == "__main__":
    main()