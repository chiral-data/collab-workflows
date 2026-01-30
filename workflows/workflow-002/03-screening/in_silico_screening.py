#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
In Silico Screening Script
Converted from Screening/In_silico_Screening.ipynb
"""

import glob
import os
import subprocess

# Read parameters from environment variables
test_mode = os.environ.get("PARAM_TEST_MODE", "true").lower() == "true"
exhaustiveness = int(os.environ.get("PARAM_EXHAUSTIVENESS", "8"))


def main():
    """Main execution function"""
    print("Starting In Silico screening...")
    print(f"Test mode: {test_mode}")
    print(f"Exhaustiveness: {exhaustiveness}")

    # Check current directory
    print(f"Current directory: {os.getcwd()}")

    # Check smina
    check_smina()

    # Create docking results directory
    create_docking_results_directory()

    # Run docking
    run_docking_screening()

    print("In Silico screening completed.")


def check_smina():
    """Check smina command availability"""
    print("\n=== Checking smina command ===")

    smina_path = "smina"

    try:
        result = subprocess.run(
            [smina_path, "--help"], capture_output=True, text=True, timeout=10
        )
        print("smina command is available.")
        print("Help output (first 500 chars):")
        print(
            result.stdout[:500] + "..." if len(result.stdout) > 500 else result.stdout
        )
    except subprocess.TimeoutExpired:
        print("smina command check timed out.")
    except FileNotFoundError:
        print(f"Error: {smina_path} not found.")
        print("Please verify smina is correctly installed.")
    except Exception as e:
        print(f"Error checking smina command: {e}")


def create_docking_results_directory():
    """Create docking results directory"""
    print("\n=== Creating docking results directory ===")

    os.makedirs("docking_results", exist_ok=True)
    print("Created docking_results directory.")


def run_docking_screening():
    """Run docking screening"""
    print("\n=== Running docking screening ===")

    # Select ligands based on test_mode
    if test_mode:
        ligands = glob.glob("./constructed_library/clean_drug108*.sdf")
        print("Test mode: screening 11 compounds")
    else:
        ligands = glob.glob("./constructed_library/clean_drug*.sdf")
        print("Full mode: screening all compounds")

    print(f"Number of ligand files found: {len(ligands)}")

    smina_path = "smina"

    # Run docking for each ligand
    for i, lig in enumerate(ligands, 1):
        fname = os.path.splitext(os.path.basename(lig))[0]
        out_sdf = f"docking_results/{fname}_docked.sdf"
        out_log = f"docking_results/{fname}_log.txt"

        print(f"\n[{i}/{len(ligands)}] Docking {fname}...")

        try:
            result = subprocess.run(
                [
                    smina_path,
                    "-r",
                    "./5Y7J_AB_chains_fixed.pdb",
                    "-l",
                    lig,
                    "--config",
                    "config.txt",
                    "-o",
                    out_sdf,
                    "--log",
                    out_log,
                    "--scoring",
                    "vina",
                    "--exhaustiveness",
                    str(exhaustiveness),
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=300,
            )

            print(f"Docking {fname} completed.")

            # Extract and display binding energy from log
            affinity = extract_affinity_from_log(out_log)
            if affinity is not None:
                print(f"  Binding energy: {affinity:.2f} kcal/mol")

        except subprocess.TimeoutExpired:
            print(f"Docking {fname} timed out.")
        except subprocess.CalledProcessError as e:
            print(f"Error docking {fname}: {e}")
            if e.stderr:
                print(f"  Details: {e.stderr[:200]}...")
        except Exception as e:
            print(f"Unexpected error docking {fname}: {e}")

    print("\nDocking screening completed.")
    print("Results saved to docking_results/ directory.")


def extract_affinity_from_log(log_file):
    """Extract binding energy from log file"""
    try:
        with open(log_file, "r") as f:
            for line in f:
                if line.strip().startswith("1 "):
                    parts = line.split()
                    if len(parts) > 1:
                        try:
                            affinity = float(parts[1])
                            return affinity
                        except ValueError:
                            pass
    except Exception as e:
        print(f"  Error reading log file: {e}")
    return None


def get_docking_summary():
    """Get docking results summary"""
    print("\n=== Docking Results Summary ===")

    log_files = glob.glob("docking_results/*_log.txt")
    results = []

    for log_file in log_files:
        affinity = extract_affinity_from_log(log_file)
        if affinity is not None:
            compound_name = os.path.basename(log_file).replace("_log.txt", "")
            results.append((compound_name, affinity))

    results.sort(key=lambda x: x[1])

    print(f"Successful dockings: {len(results)}")
    print("\nTop 5 results:")
    for i, (compound, affinity) in enumerate(results[:5], 1):
        print(f"{i}: {compound} - {affinity:.2f} kcal/mol")

    return results


if __name__ == "__main__":
    main()
