#!/usr/bin/env python3
"""
Generate report with top compounds and their binding energies.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def get_energy_status(binding_energy: float) -> str:
    """Get status label for binding energy."""
    if binding_energy < -8:
        return "EXCELLENT"
    elif binding_energy < -5:
        return "GOOD"
    elif binding_energy < 0:
        return "FAIR"
    else:
        return "POOR"


def generate_report(input_path: str, top_n: int, output_dir: str) -> bool:
    """
    Generate report with top compounds.

    Args:
        input_path: Path to evaluated chemical space SDF file
        top_n: Number of top compounds to extract
        output_dir: Output directory for report files

    Returns:
        True if successful, False otherwise
    """
    logger.info(f"Generating report for top {top_n} compounds...")

    input_file = Path(input_path)
    if not input_file.exists():
        logger.error(f"Input file not found: {input_path}")
        return False

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load molecules from SDF
    molecules = []
    with Chem.SDMolSupplier(str(input_file)) as SDF:
        for mol in SDF:
            if mol is None:
                continue
            props = mol.GetPropsAsDict()
            # Check for valid score (prioritize score over minimizedAffinity)
            score = props.get("score") or props.get("minimizedAffinity")
            if score and str(score) not in ["<NA>", "nan", "None", ""]:
                molecules.append(mol)

    logger.info(f"Loaded {len(molecules)} successful molecules")

    if len(molecules) == 0:
        logger.error("No successful molecules found")
        return False

    # Sort by score (binding energy) - lower is better
    def get_score(mol):
        props = mol.GetPropsAsDict()
        # Prioritize score (pK) and convert to kcal/mol, fallback to minimizedAffinity
        if props.get("score") and str(props.get("score")) not in ["<NA>", "nan", "None", ""]:
            try:
                pk_value = float(props.get("score"))
                # Convert pK to ΔG (kcal/mol): ΔG ≈ -1.36 × pK at 298K
                return -1.36 * pk_value
            except (ValueError, TypeError):
                pass

        # Fallback to minimizedAffinity if score is not available
        score = props.get("minimizedAffinity", "999")
        try:
            return float(score)
        except (ValueError, TypeError):
            return 999

    sorted_molecules = sorted(molecules, key=get_score, reverse=False)

    # Get top compounds
    top_compounds = sorted_molecules[:top_n]
    logger.info(f"Selected top {len(top_compounds)} compounds")

    # Save top compounds SDF
    top_sdf_path = output_path / "top_compounds.sdf"
    with Chem.SDWriter(str(top_sdf_path)) as SDF_OUT:
        for mol in top_compounds:
            SDF_OUT.write(mol)
    logger.info(f"Top compounds SDF saved to: {top_sdf_path}")

    # Create detailed report
    report_data = []
    for i, mol in enumerate(top_compounds, 1):
        props = mol.GetPropsAsDict()
        binding_energy = get_score(mol)
        report_data.append(
            {
                "Rank": i,
                "Index": props.get("index", "N/A"),
                "Binding_Energy": binding_energy,
                "Status": get_energy_status(binding_energy),
                "Success": props.get("Success", "N/A"),
                "SMILES": Chem.MolToSmiles(mol) if mol else "N/A",
                "Molecular_Weight": rdMolDescriptors.CalcExactMolWt(mol)
                if mol
                else "N/A",
                "LogP": rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
                if mol
                else "N/A",
            }
        )

    # Save CSV report
    report_df = pd.DataFrame(report_data)
    csv_path = output_path / "top_compounds_report.csv"
    report_df.to_csv(csv_path, index=False)
    logger.info(f"CSV report saved to: {csv_path}")

    # Create summary text
    summary_lines = [
        "=" * 80,
        "TOP COMPOUNDS SUMMARY (Lower binding energy is better)",
        "=" * 80,
        "",
    ]

    for data in report_data:
        binding_energy = data["Binding_Energy"]
        summary_lines.append(
            f"Rank {data['Rank']:2d}: Binding Energy = {binding_energy:8.3f} ({data['Status']}), "
            f"MW = {data['Molecular_Weight']:6.1f}, LogP = {data['LogP']:5.2f}"
        )

    summary_lines.extend(["", "=" * 80, f"Total molecules analyzed: {len(molecules)}"])

    # Calculate statistics
    all_scores = [get_score(m) for m in molecules if get_score(m) < 999]
    if all_scores:
        summary_lines.extend(
            [
                f"Best binding energy: {min(all_scores):.3f}",
                f"Worst binding energy: {max(all_scores):.3f}",
                f"Average binding energy: {sum(all_scores)/len(all_scores):.3f}",
            ]
        )

    summary_lines.append("=" * 80)

    summary_text = "\n".join(summary_lines)
    summary_path = output_path / "summary.txt"
    with open(summary_path, "w") as f:
        f.write(summary_text)
    logger.info(f"Summary saved to: {summary_path}")

    # Print summary to console
    print("\n" + summary_text + "\n")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Generate report with top compounds"
    )
    parser.add_argument(
        "--input", required=True, help="Input evaluated chemical space SDF file"
    )
    parser.add_argument(
        "--top-n", type=int, default=10, help="Number of top compounds to extract"
    )
    parser.add_argument(
        "--output-dir", default=".", help="Output directory for report files"
    )

    args = parser.parse_args()

    success = generate_report(args.input, args.top_n, args.output_dir)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
