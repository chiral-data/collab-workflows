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

    # Load the ligand molecule from SMILES (preserves atom map numbers)
    with open(ligand_file, "r") as f:
        smiles = f.read().strip()
        logger.info(f"Loaded molecule with smiles: {smiles}")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.error("Failed to load molecule from SMILES file")
        return False

    # Add hydrogens
    mol = Chem.AddHs(mol)
    logger.info(f"Loaded molecule with {mol.GetNumAtoms()} atoms")

    # Create FEgrow scaffold
    scaffold = fegrow.RMol(mol)

    # Set attachment point (mark atom as dummy)
    scaffold.GetAtomWithIdx(attachment_id).SetAtomicNum(0)
    logger.info(f"Set attachment point at atom index {attachment_id}")

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
        default=27,
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
