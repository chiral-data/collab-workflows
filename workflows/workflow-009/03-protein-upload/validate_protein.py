#!/usr/bin/env python3
"""
Validate and prepare the protein PDB file.
"""

import argparse
import logging
import shutil
import sys
from pathlib import Path

import prody

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def validate_protein(input_path: str, output_path: str) -> bool:
    """
    Validate and copy the protein PDB file.

    Args:
        input_path: Path to input protein PDB file
        output_path: Path to output protein PDB file

    Returns:
        True if validation successful, False otherwise
    """
    logger.info(f"Validating protein file: {input_path}")

    input_file = Path(input_path)
    if not input_file.exists():
        logger.error(f"Protein file not found: {input_path}")
        return False

    # Load and validate the protein
    try:
        protein = prody.parsePDB(str(input_file))
    except Exception as e:
        logger.error(f"Failed to parse PDB file: {e}")
        return False

    if protein is None:
        logger.error("Failed to load protein from PDB file")
        return False

    num_atoms = protein.numAtoms()
    num_residues = protein.numResidues()
    logger.info(f"Loaded protein with {num_atoms} atoms and {num_residues} residues")

    # Check for protein atoms
    protein_atoms = protein.select("protein")
    if protein_atoms is None:
        logger.error("No protein atoms found in the PDB file")
        return False

    logger.info(f"Found {protein_atoms.numAtoms()} protein atoms")

    # Copy the validated protein
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(input_file, output_file)

    logger.info(f"Validated protein saved to: {output_path}")
    return True


def main():
    parser = argparse.ArgumentParser(description="Validate protein PDB file")
    parser.add_argument("--input", required=True, help="Input protein PDB file")
    parser.add_argument(
        "--output", default="protein.pdb", help="Output validated protein PDB file"
    )

    args = parser.parse_args()

    success = validate_protein(args.input, args.output)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
