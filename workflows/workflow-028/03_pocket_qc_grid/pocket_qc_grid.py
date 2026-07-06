#!/usr/bin/env python3
"""
Node 03: Pocket QC + Grid preparation.

Reads P2Rank pocket predictions alongside Boltz-2 pLDDT confidence to filter
pockets by per-residue structural quality, then prepares the three files that
Uni-Mol Docking V2 requires:

  receptor.pdb  — protein-only PDB  (gemmi mmCIF→PDB; ligand chains removed)
  ligand.sdf    — 3D conformer      (RDKit ETKDGv3 + MMFF from SMILES)
  grid.json     — docking grid      ({center_x/y/z, size_x/y/z} for Uni-Mol)

Also writes:
  pocket_qc.json — per-pocket pLDDT stats for the report node

pLDDT caveat (Eguida & Rognan 2023): threshold ≥70 is necessary but not
sufficient — high pLDDT does not guarantee correct side-chain geometry at the
pocket. Pockets below threshold are flagged in pocket_qc.json but the script
does not abort; the selected pocket's QC status is surfaced to the report.
"""

import argparse
import glob
import json
import math
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# pLDDT helpers
# ---------------------------------------------------------------------------

def load_plddt(model_id: str) -> np.ndarray:
    """
    Load the per-residue pLDDT array for the selected Boltz-2 model.
    Returns values on the 0-100 scale.
    """
    candidates = glob.glob(f"plddt_{model_id}.npz") + glob.glob("plddt_*.npz")
    if not candidates:
        raise FileNotFoundError("No plddt_*.npz file found in inputs/")

    npz = np.load(candidates[0])
    key = "plddt" if "plddt" in npz.files else npz.files[0]
    arr = npz[key].astype(float).flatten()

    # Boltz-2 writes pLDDT in 0-1 range; convert to 0-100 for readability
    if arr.max() <= 1.01:
        arr = arr * 100.0
    return arr


def residue_nums_from_ids(ids_str: str) -> list[int]:
    """
    Extract 1-indexed residue numbers from a P2Rank residue_ids string.
    Handles formats like 'A_123 A_124', 'A_123_ALA A_124_GLY', or '123 124'.
    """
    nums = []
    for token in str(ids_str).split():
        parts = token.split("_")
        for part in parts:
            if part.isdigit():
                nums.append(int(part))
                break
    return nums


def pocket_plddt_stats(residue_ids_str: str, plddt: np.ndarray) -> dict:
    """Compute mean/min/std pLDDT for the residues belonging to a pocket."""
    resnums = residue_nums_from_ids(residue_ids_str)
    # plddt is 0-indexed; residue numbers from structure are 1-indexed
    idxs = [r - 1 for r in resnums if 0 <= r - 1 < len(plddt)]
    if not idxs:
        return {"mean": None, "min": None, "std": None, "n_residues": 0}
    vals = plddt[idxs]
    return {
        "mean": round(float(vals.mean()), 2),
        "min":  round(float(vals.min()),  2),
        "std":  round(float(vals.std()),  2),
        "n_residues": len(idxs),
    }


# ---------------------------------------------------------------------------
# P2Rank predictions
# ---------------------------------------------------------------------------

def load_predictions() -> pd.DataFrame:
    """Load P2Rank predictions.csv; normalise column names."""
    df = pd.read_csv("predictions.csv", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    return df


# ---------------------------------------------------------------------------
# mmCIF → PDB conversion
# ---------------------------------------------------------------------------

def cif_to_pdb(cif_path: str, pdb_path: str, protein_chain_ids: list[str]) -> None:
    """
    Convert Boltz-2 mmCIF to PDB keeping only protein chains.
    Uses gemmi; avoids BioPython MMCIFIO which is incompatible with
    Boltz-2 CIF files (missing entity.id field).
    """
    import gemmi

    st = gemmi.read_structure(cif_path)
    model = st[0]

    # Remove any chain not in the protein chain list (e.g. ligand chain B)
    chains_to_drop = [c.name for c in model if c.name not in protein_chain_ids]
    for name in chains_to_drop:
        model.remove_chain(name)

    st.write_pdb(pdb_path)
    print(f"  Wrote receptor PDB: {pdb_path} ({len(list(model))} chains retained)", flush=True)


# ---------------------------------------------------------------------------
# SMILES → SDF conversion
# ---------------------------------------------------------------------------

def smiles_to_sdf(smiles: str, sdf_path: str) -> None:
    """
    Generate a 3D conformer from SMILES and write to SDF.
    Uses RDKit ETKDGv3 geometry embedding followed by MMFF optimisation.
    Raises RuntimeError if embedding fails (invalid SMILES or geometry failure).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result != 0:
        raise RuntimeError(f"3D embedding failed for SMILES: {smiles}")

    AllChem.MMFFOptimizeMolecule(mol)
    mol = Chem.RemoveHs(mol)

    with Chem.SDWriter(sdf_path) as w:
        w.write(mol)
    print(f"  Wrote ligand SDF: {sdf_path}", flush=True)


# ---------------------------------------------------------------------------
# Grid JSON
# ---------------------------------------------------------------------------

def write_grid(center_x: float, center_y: float, center_z: float,
               box_size: float, path: str) -> None:
    """Write Uni-Mol Docking V2 grid JSON."""
    grid = {
        "center_x": round(center_x, 3),
        "center_y": round(center_y, 3),
        "center_z": round(center_z, 3),
        "size_x":   round(box_size, 3),
        "size_y":   round(box_size, 3),
        "size_z":   round(box_size, 3),
    }
    with open(path, "w") as f:
        json.dump(grid, f, indent=2)
    print(f"  Wrote grid JSON: {path}  center=({center_x:.1f}, {center_y:.1f}, {center_z:.1f})  box={box_size}Å", flush=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plddt-threshold", type=float, default=70.0)
    parser.add_argument("--box-size",        type=float, default=22.5)
    parser.add_argument("--pocket-rank",     type=int,   default=1)
    args = parser.parse_args()

    os.makedirs("./outputs", exist_ok=True)

    # --- Model ID -----------------------------------------------------------
    model_id = Path("selected_model_id.txt").read_text().strip()
    print(f"Selected Boltz-2 model: {model_id}", flush=True)

    # --- pLDDT array --------------------------------------------------------
    plddt = load_plddt(model_id)
    print(f"Loaded pLDDT array: {len(plddt)} tokens, range [{plddt.min():.1f}, {plddt.max():.1f}]", flush=True)

    # --- Pocket predictions -------------------------------------------------
    pockets_df = load_predictions()
    print(f"P2Rank pockets: {len(pockets_df)}", flush=True)

    # --- Per-pocket pLDDT QC ------------------------------------------------
    qc_records = []
    for _, row in pockets_df.iterrows():
        residue_ids_str = row.get("residue_ids", row.get("residue_ids ", ""))
        stats = pocket_plddt_stats(str(residue_ids_str), plddt)
        mean_plddt = stats["mean"]
        passes = (mean_plddt is not None) and (mean_plddt >= args.plddt_threshold)
        record = {
            "rank":         int(row["rank"]),
            "name":         str(row["name"]),
            "p2rank_score": round(float(row["score"]), 4),
            "probability":  round(float(row["probability"]), 4),
            "center_x":     round(float(row["center_x"]), 3),
            "center_y":     round(float(row["center_y"]), 3),
            "center_z":     round(float(row["center_z"]), 3),
            "plddt_mean":   stats["mean"],
            "plddt_min":    stats["min"],
            "plddt_std":    stats["std"],
            "n_residues":   stats["n_residues"],
            "plddt_passes": passes,
            "plddt_threshold": args.plddt_threshold,
        }
        qc_records.append(record)
        status = "PASS" if passes else "FAIL"
        print(
            f"  Pocket {record['rank']:2d} ({record['name']}): "
            f"score={record['p2rank_score']:.3f}  "
            f"pLDDT mean={mean_plddt if mean_plddt is not None else 'N/A'}  [{status}]",
            flush=True,
        )

    # --- Select target pocket -----------------------------------------------
    # Use requested pocket_rank; warn if it fails QC but don't abort.
    if not qc_records:
        print(
            "ERROR: P2Rank found no pockets in the predicted structure. "
            "The protein may be too short or disordered for pocket detection.",
            flush=True,
        )
        sys.exit(1)

    target_rank = args.pocket_rank
    selected = next((r for r in qc_records if r["rank"] == target_rank), None)
    if selected is None:
        selected = qc_records[0]
        print(f"WARNING: pocket_rank={target_rank} not found; falling back to rank 1", flush=True)

    if not selected["plddt_passes"]:
        print(
            f"WARNING: Selected pocket (rank {selected['rank']}) has mean pLDDT "
            f"{selected['plddt_mean']} < threshold {args.plddt_threshold}. "
            f"Proceeding — see pocket_qc.json. Consider increasing pocket_rank "
            f"or reviewing Boltz-2 confidence before trusting docking results. "
            f"(Eguida & Rognan 2023: pLDDT ≥70 is necessary but not sufficient.)",
            flush=True,
        )

    print(
        f"Selected pocket: rank={selected['rank']}  "
        f"center=({selected['center_x']}, {selected['center_y']}, {selected['center_z']})",
        flush=True,
    )

    # Write pocket_qc.json
    qc_output = {
        "threshold": args.plddt_threshold,
        "selected_pocket_rank": selected["rank"],
        "selected_pocket_passes_qc": selected["plddt_passes"],
        "pockets": qc_records,
        "caveat": (
            "pLDDT ≥70 is necessary but not sufficient for reliable docking. "
            "High pLDDT filters disordered regions but does not guarantee correct "
            "side-chain geometry (Eguida & Rognan 2023, JCIM, PMC9852548)."
        ),
    }
    with open("./outputs/pocket_qc.json", "w") as f:
        json.dump(qc_output, f, indent=2)
    print("  Wrote pocket_qc.json", flush=True)

    # --- mmCIF → PDB --------------------------------------------------------
    with open("input_summary.json") as f:
        summary = json.load(f)
    protein_chains = [e["id"] for e in summary.get("entities", []) if e.get("type") == "protein"]
    if not protein_chains:
        protein_chains = ["A"]  # safe default
    print(f"Protein chains to retain: {protein_chains}", flush=True)

    cif_to_pdb("selected_structure.cif", "./outputs/receptor.pdb", protein_chains)

    # --- SMILES → SDF -------------------------------------------------------
    ligand_entities = [e for e in summary.get("entities", []) if e.get("type") == "ligand"]
    if not ligand_entities:
        print("ERROR: No ligand entity found in input_summary.json", flush=True)
        sys.exit(1)
    smiles = ligand_entities[0].get("smiles", "")
    if not smiles:
        print("ERROR: Ligand entity has no SMILES in input_summary.json", flush=True)
        sys.exit(1)
    print(f"Ligand SMILES: {smiles}", flush=True)
    smiles_to_sdf(smiles, "./outputs/ligand.sdf")

    # --- Grid JSON ----------------------------------------------------------
    write_grid(
        center_x=selected["center_x"],
        center_y=selected["center_y"],
        center_z=selected["center_z"],
        box_size=args.box_size,
        path="./outputs/grid.json",
    )

    print("Node 03 complete.", flush=True)


if __name__ == "__main__":
    main()
