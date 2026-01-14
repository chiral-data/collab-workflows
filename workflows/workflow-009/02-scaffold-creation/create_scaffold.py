#!/usr/bin/env python3
"""
Create scaffold molecule with attachment point from the ligand.
"""

import argparse
import logging
import pickle
import sys
from pathlib import Path

import fegrow
from rdkit import Chem

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def create_scaffold(ligand_path: str, attachment_id: int, output_path: str) -> bool:
    """
    Create scaffold from ligand with specified attachment point.

    Args:
        ligand_path: Path to input ligand SMILES file (with atom map numbers)
        attachment_id: Atom map number for attachment point (from Node 01 visualization)
        output_path: Path to save scaffold pickle file

    Returns:
        True if successful, False otherwise
    """
    logger.info(f"Creating scaffold from ligand: {ligand_path}")
    logger.info(f"Attachment point atom map number: {attachment_id}")

    ligand_file = Path(ligand_path)
    if not ligand_file.exists():
        logger.error(f"Ligand file not found: {ligand_path}")
        return False

    # Load SMILES to get atom map number -> index mapping
    with open(ligand_file, "r") as f:
        smiles = f.read().strip()
        logger.info(f"Loaded SMILES: {smiles}")

    smiles_mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if smiles_mol is None:
        logger.error("Failed to load molecule from SMILES file")
        return False
    Chem.SanitizeMol(
        smiles_mol,
        sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
        ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS,
    )

    # Find atom index by map number in SMILES molecule
    atom_idx = None
    for atom in smiles_mol.GetAtoms():
        if atom.GetAtomMapNum() == attachment_id:
            atom_idx = atom.GetIdx()
            break

    if atom_idx is None:
        logger.error(f"No atom found with map number {attachment_id}")
        return False

    logger.info(f"Found atom with map number {attachment_id} at index {atom_idx} (in SMILES)")

    # Load SDF for 3D coordinates (same directory, .sdf extension)
    sdf_file = ligand_file.with_suffix(".sdf")
    if not sdf_file.exists():
        logger.error(f"SDF file not found: {sdf_file}")
        return False

    suppl = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
    mol = next(iter(suppl), None)
    if mol is None:
        logger.error("Failed to load molecule from SDF file")
        return False

    logger.info(f"Loaded SDF with {mol.GetNumAtoms()} atoms and 3D coordinates")

    # Check if SDF has atom map numbers
    sdf_map_nums = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
    logger.info(f"SDF atom map numbers: {sdf_map_nums[:10]}... (first 10)")

    # Find the atom with the same map number in SDF
    sdf_atom_idx = None
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == attachment_id:
            sdf_atom_idx = atom.GetIdx()
            break

    if sdf_atom_idx is not None:
        logger.info(f"Found atom with map number {attachment_id} at index {sdf_atom_idx} (in SDF)")
        atom_idx = sdf_atom_idx  # Use SDF index
    else:
        logger.warning(f"Map number {attachment_id} not found in SDF, using SMILES index {atom_idx}")

    # Verify atom count matches
    if mol.GetNumAtoms() != smiles_mol.GetNumAtoms():
        logger.error(f"Atom count mismatch: SDF has {mol.GetNumAtoms()}, SMILES has {smiles_mol.GetNumAtoms()}")
        return False

    # Create FEgrow scaffold from SDF molecule (has 3D coords)
    scaffold = fegrow.RMol(mol)

    # Set attachment point (mark atom as dummy)
    scaffold.GetAtomWithIdx(atom_idx).SetAtomicNum(0)
    logger.info(
        f"Set attachment point at atom index {atom_idx} (map number {attachment_id})"
    )

    # Save scaffold as pickle
    output_file = Path(output_path)
    with open(output_file, "wb") as f:
        pickle.dump(scaffold, f)

    logger.info(f"Scaffold saved to: {output_path}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Create scaffold from ligand with attachment point"
    )
    parser.add_argument(
        "--ligand",
        required=True,
        help="Input ligand SMILES file (with atom map numbers)",
    )
    parser.add_argument(
        "--attachment-id",
        type=int,
        default=26,
        help="Atom index for attachment point",
    )
    parser.add_argument(
        "--output", default="scaffold.pkl", help="Output scaffold pickle file"
    )

    args = parser.parse_args()

    success = create_scaffold(args.ligand, args.attachment_id, args.output)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
