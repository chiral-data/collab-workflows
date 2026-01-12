#!/usr/bin/env python3
"""
Prepare protein structure for docking.
Removes water, nucleic acids, and hetero atoms, then fixes the receptor.
"""

import argparse
import logging
import sys
from pathlib import Path

import fegrow
import prody

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def prepare_protein(input_path: str, output_path: str) -> bool:
    """
    Prepare protein structure for docking.

    Args:
        input_path: Path to input protein PDB file
        output_path: Path to output prepared receptor PDB file

    Returns:
        True if successful, False otherwise
    """
    logger.info(f"Preparing protein from: {input_path}")

    input_file = Path(input_path)
    if not input_file.exists():
        logger.error(f"Protein file not found: {input_path}")
        return False

    # Load protein structure
    try:
        structure = prody.parsePDB(str(input_file))
    except Exception as e:
        logger.error(f"Failed to parse PDB file: {e}")
        return False

    logger.info(f"Loaded structure with {structure.numAtoms()} atoms")

    # Select protein atoms only (remove nucleic, hetero, water)
    receptor = structure.select("not (nucleic or hetatm or water)")

    if receptor is None:
        logger.error("No protein atoms found after selection")
        return False

    logger.info(f"Selected {receptor.numAtoms()} receptor atoms")

    # Save intermediate receptor file
    intermediate_path = Path("rec_intermediate.pdb")
    prody.writePDB(str(intermediate_path), receptor)
    logger.info(f"Intermediate receptor saved to: {intermediate_path}")

    # Fix receptor using FEgrow
    output_file = Path(output_path)
    try:
        fegrow.fix_receptor(str(intermediate_path), str(output_file))
        logger.info(f"Fixed receptor saved to: {output_path}")
    except Exception as e:
        logger.error(f"Failed to fix receptor: {e}")
        return False

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Prepare protein structure for docking"
    )
    parser.add_argument("--input", required=True, help="Input protein PDB file")
    parser.add_argument(
        "--output", default="rec_final.pdb", help="Output prepared receptor PDB file"
    )

    args = parser.parse_args()

    success = prepare_protein(args.input, args.output)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
