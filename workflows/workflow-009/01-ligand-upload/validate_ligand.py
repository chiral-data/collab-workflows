#!/usr/bin/env python3
"""
Validate and prepare the ligand file, output SMILES with atom map numbers.
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
    Validate ligand and output SMILES with atom map numbers.

    Args:
        input_path: Path to input ligand SDF file
        output_path: Path to output SMILES file

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

    # Set atom map numbers (these are preserved in SMILES for attachment point selection)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())

    # Write SMILES with atom map numbers (for attachment point selection in visualization)
    output_file = Path(output_path)
    smiles = Chem.MolToSmiles(mol)
    with open(output_file, "w") as f:
        f.write(smiles)
    logger.info(f"SMILES with atom map numbers saved to: {output_path}")

    # Also save SDF with 3D coordinates (needed for FEgrow conformer generation)
    sdf_output = output_file.with_suffix(".sdf")
    with Chem.SDWriter(str(sdf_output)) as w:
        w.write(mol)
    logger.info(f"SDF with 3D coordinates saved to: {sdf_output}")

    return True


def main():
    parser = argparse.ArgumentParser(description="Validate ligand file")
    parser.add_argument("--input", required=True, help="Input ligand SDF file")
    parser.add_argument(
        "--output",
        default="ligand.smi",
        help="Output SMILES file with atom map numbers",
    )

    args = parser.parse_args()

    success = validate_ligand(args.input, args.output)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
