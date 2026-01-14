#!/usr/bin/env python3
"""
Validate and prepare the ligand SDF file.
"""

import argparse
import logging
import sys
from pathlib import Path

from rdkit import Chem

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def validate_ligand(input_path: str, output_path: str) -> bool:
    """
    Validate and copy the ligand SDF file.

    Args:
        input_path: Path to input ligand SDF file
        output_path: Path to output ligand SDF file

    Returns:
        True if validation successful, False otherwise
    """
    logger.info(f"Validating ligand file: {input_path}")

    input_file = Path(input_path)
    if not input_file.exists():
        logger.error(f"Ligand file not found: {input_path}")
        return False

    # Load and validate the molecule
    suppl = Chem.SDMolSupplier(str(input_file), removeHs=False)
    molecules = [mol for mol in suppl if mol is not None]

    if len(molecules) == 0:
        logger.error("No valid molecules found in the SDF file")
        return False

    mol = molecules[0]
    logger.info(f"Loaded molecule with {mol.GetNumAtoms()} atoms")

    # Add hydrogens if not present
    mol = Chem.AddHs(mol)
    logger.info(f"After adding hydrogens: {mol.GetNumAtoms()} atoms")

    # Set atom map numbers to preserve atom indices
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())

    # Write the validated ligand as SDF
    output_file = Path(output_path)
    with Chem.SDWriter(str(output_file)) as writer:
        writer.write(mol)
    logger.info(f"Validated ligand saved to: {output_path}")

    # Also write SMILES with atom map numbers for Node 02
    smiles_path = output_file.with_suffix(".smi")
    smiles = Chem.MolToSmiles(mol)
    with open(smiles_path, "w") as f:
        f.write(smiles)
    logger.info(f"SMILES with atom indices saved to: {smiles_path}")

    return True


def main():
    parser = argparse.ArgumentParser(description="Validate ligand SDF file")
    parser.add_argument("--input", required=True, help="Input ligand SDF file")
    parser.add_argument(
        "--output", default="ligand.sdf", help="Output validated ligand SDF file"
    )

    args = parser.parse_args()

    success = validate_ligand(args.input, args.output)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
