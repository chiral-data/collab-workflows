#!/usr/bin/env python3
"""
Node 04: Uni-Mol Docking V2.

Runs protein-ligand docking using Uni-Mol Docking V2, a diffusion-style
deep-learning docking engine that achieves 77.6% RMSD < 2 Å on the PoseBusters
benchmark — significantly outperforming traditional physics-based methods on
predicted structures.

Important constraints:
  - Outputs 3D poses (SDF) only. The internal prmsd_score used for pose
    ranking is not written to the output file (Uni-Mol Docking V2 design).
  - Model weights (464 MB) must be pre-downloaded; see Dockerfile for instructions.
  - RDKit sanitization in Uni-Mol can silently drop invalid ligands
    (Uni-Mol #281). This script validates the ligand before docking.

Invocation:
  Docking is run via interface/demo.py in the cloned Uni-Mol repository
  (/opt/unimol/unimol_docking_v2/interface/demo.py). The docking grid is
  passed as a JSON file via --input-docking-grid rather than as individual
  coordinate/radius flags.
"""

import argparse
import glob
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import List


PREDICT_SCRIPT = "/opt/unimol/unimol_docking_v2/interface/demo.py"


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

def validate_receptor(pdb_path: str) -> None:
    path = Path(pdb_path)
    if not path.exists():
        raise FileNotFoundError(f"Receptor PDB not found: {pdb_path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Receptor PDB is empty: {pdb_path}")
    # Basic check: file should contain ATOM records
    content = path.read_text()
    if "ATOM" not in content and "HETATM" not in content:
        raise ValueError(f"No ATOM/HETATM records in receptor PDB: {pdb_path}")
    print(f"  Receptor: {pdb_path} ({path.stat().st_size // 1024} KB)", flush=True)


def validate_ligand(sdf_path: str) -> str:
    """
    Validate ligand SDF with RDKit before passing to Uni-Mol.
    Uni-Mol #281: RDKit sanitization inside Uni-Mol can silently drop invalid
    ligands, producing empty output. Pre-validating here catches the problem
    with a clear error rather than silent failure.
    Returns the canonical SMILES of the parsed molecule.
    """
    from rdkit import Chem

    path = Path(sdf_path)
    if not path.exists():
        raise FileNotFoundError(f"Ligand SDF not found: {sdf_path}")

    suppl = Chem.SDMolSupplier(str(path), removeHs=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        raise ValueError(
            f"RDKit could not parse any molecule from {sdf_path}. "
            "Check that the SDF is valid — Uni-Mol will silently drop it otherwise (Uni-Mol #281)."
        )

    smiles = Chem.MolToSmiles(mols[0])
    print(f"  Ligand: {sdf_path}  SMILES={smiles}", flush=True)
    return smiles


def validate_weights(weights_path: str) -> str:
    """Validate weights directory exists and return the path to the weight file."""
    path = Path(weights_path)
    if not path.exists():
        print(
            f"\nERROR: Model weights not found at {weights_path}\n"
            "Uni-Mol Docking V2 weights (464 MB) must be downloaded separately.\n"
            "1. Visit: https://github.com/deepmodeling/Uni-Mol/tree/main/unimol_docking_v2\n"
            "2. Download the weights from the Dropbox link in the README.\n"
            "3. Place the files at the path given by PARAM_WEIGHTS_PATH.\n",
            file=sys.stderr,
            flush=True,
        )
        sys.exit(1)
    weight_files = list(path.glob("*.pt")) + list(path.glob("*.pkl"))
    if not weight_files:
        raise FileNotFoundError(f"No .pt or .pkl weight files found in {weights_path}")
    resolved = str(weight_files[0])
    print(f"  Weights: {resolved}", flush=True)
    return resolved


# ---------------------------------------------------------------------------
# Docking
# ---------------------------------------------------------------------------

def _build_docking_cmd(
    receptor_pdb: str,
    ligand_sdf: str,
    grid_json: str,
    num_poses: int,
    weights_path: str,
    output_dir: str,
) -> List[str]:
    """Build the Uni-Mol Docking V2 interface/demo.py command."""
    return [
        sys.executable, PREDICT_SCRIPT,
        "--mode",               "single",
        "--input-protein",      receptor_pdb,
        "--input-ligand",       ligand_sdf,
        "--input-docking-grid", grid_json,
        "--conf-size",          str(num_poses),
        "--model-dir",          weights_path,
        "--output-ligand-dir",  output_dir,
        "--steric-clash-fix",
        "--cluster",
    ]


def run_docking(
    receptor_pdb: str,
    ligand_sdf: str,
    grid_json: str,
    num_poses: int,
    weights_path: str,
    output_dir: str,
) -> None:
    os.makedirs(output_dir, exist_ok=True)

    cmd = _build_docking_cmd(
        receptor_pdb, ligand_sdf,
        grid_json,
        num_poses, weights_path, output_dir,
    )

    print(f"Running: {' '.join(cmd)}", flush=True)
    result = subprocess.run(cmd, capture_output=False, text=True)
    if result.returncode != 0:
        print(f"ERROR: Uni-Mol docking exited with code {result.returncode}", flush=True)
        sys.exit(1)


def collect_poses(output_dir: str, dest_sdf: str) -> int:
    """
    Gather Uni-Mol output SDF files into a single docked_poses.sdf.
    Returns the number of poses collected.
    """
    sdf_files = sorted(set(
        glob.glob(os.path.join(output_dir, "**", "*.sdf"), recursive=True)
    ))

    if not sdf_files:
        print(
            f"WARNING: No SDF output found in {output_dir}. "
            "Check that the docking ran successfully and the weights are correct.",
            flush=True,
        )
        # Write an empty placeholder so downstream nodes don't crash
        Path(dest_sdf).write_text("")
        return 0

    # Concatenate all SDF files into one
    with open(dest_sdf, "w") as out:
        for sdf in sdf_files:
            out.write(Path(sdf).read_text())

    # Count poses: prefer $$$$ delimiter, fall back to M  END blocks
    content = Path(dest_sdf).read_text()
    n_poses = content.count("$$$$")
    if n_poses == 0:
        n_poses = content.count("M  END")
    print(f"Collected {n_poses} pose(s) from {len(sdf_files)} SDF file(s) → {dest_sdf}", flush=True)
    return n_poses


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--num-poses",    type=int,   default=10)
    parser.add_argument("--weights-path", type=str,   default="/opt/unimol_weights")
    args = parser.parse_args()

    os.makedirs("./outputs", exist_ok=True)

    print("Validating inputs...", flush=True)
    validate_receptor("receptor.pdb")
    ligand_smiles = validate_ligand("ligand.sdf")
    weight_file = validate_weights(args.weights_path)

    with open("grid.json") as f:
        grid = json.load(f)
    if Path("pocket_qc.json").exists():
        with open("pocket_qc.json") as f:
            pocket_qc = json.load(f)
    else:
        pocket_qc = {}
    print(
        f"Grid: center=({grid['center_x']}, {grid['center_y']}, {grid['center_z']})  "
        f"box={grid['size_x']}Å",
        flush=True,
    )

    print(f"\nRunning Uni-Mol Docking V2 ({args.num_poses} poses)...", flush=True)
    run_docking(
        receptor_pdb="receptor.pdb",
        ligand_sdf="ligand.sdf",
        grid_json="grid.json",
        num_poses=args.num_poses,
        weights_path=weight_file,
        output_dir="./unimol_output",
    )

    n_poses = collect_poses("./unimol_output", "./outputs/docked_poses.sdf")

    # Write summary for the report node
    summary = {
        "receptor_pdb":       "receptor.pdb",
        "ligand_sdf":         "ligand.sdf",
        "ligand_smiles":      ligand_smiles,
        "grid":               grid,
        "num_poses_requested": args.num_poses,
        "num_poses_generated": n_poses,
        "output_sdf":         "docked_poses.sdf",
        "pocket_qc_passed":   pocket_qc.get("selected_pocket_passes_qc"),
        "selected_pocket_rank": pocket_qc.get("selected_pocket_rank"),
        "note": (
            "Uni-Mol Docking V2 outputs 3D poses only. "
            "No binding affinity or confidence score is written to the SDF — "
            "the internal prmsd_score used for pose ranking is not exposed in output. "
            "Use Boltz-2 affinity_*.json for a complementary affinity estimate."
        ),
    }
    with open("./outputs/docking_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("Wrote docking_summary.json", flush=True)


if __name__ == "__main__":
    main()
