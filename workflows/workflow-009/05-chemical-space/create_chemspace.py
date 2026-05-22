#!/usr/bin/env python3
"""
Create chemical space by combining scaffold with linkers and R-groups.
"""

import argparse
import logging
import cloudpickle
import sys
import os
from pathlib import Path

from dask.distributed import LocalCluster
from fegrow import ChemSpace, Linkers, RGroups

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def create_chemical_space(
    scaffold_path: str,
    protein_path: str,
    num_linkers: int,
    num_rgroups: int,
    output_path: str,
) -> bool:
    """
    Create chemical space from scaffold, protein, linkers and R-groups.

    Args:
        scaffold_path: Path to scaffold pickle file
        protein_path: Path to prepared protein PDB file
        num_linkers: Number of linkers to use
        num_rgroups: Number of R-groups to use
        output_path: Path to output chemical space pickle file

    Returns:
        True if successful, False otherwise
    """
    logger.info("Creating chemical space...")
    logger.info(f"Using {num_linkers} linkers and {num_rgroups} R-groups")

    # Load scaffold
    scaffold_file = Path(scaffold_path)
    if not scaffold_file.exists():
        logger.error(f"Scaffold file not found: {scaffold_path}")
        return False

    with open(scaffold_file, "rb") as f:
        scaffold = cloudpickle.load(f)
    logger.info("Scaffold loaded successfully")

    # Check protein file
    protein_file = Path(protein_path)
    if not protein_file.exists():
        logger.error(f"Protein file not found: {protein_path}")
        return False

    # Setup Dask cluster
    logger.info("Setting up Dask cluster...")
    n_workers = os.cpu_count() or 2  # Use all available CPU cores, fallback to 2
    logger.info(f"Configuring LocalCluster with {n_workers} workers")
    lc = LocalCluster(
        processes=True,
        n_workers=n_workers,
        threads_per_worker=1,
        memory_limit="2GB",
        dashboard_address=":0",  # Use random port to avoid conflicts
        silence_logs=False,
    )

    try:
        # Initialize chemical space
        cs = ChemSpace(dask_cluster=lc)
        cs.set_dask_caching()

        # Add scaffold and protein
        cs.add_scaffold(scaffold)
        cs.add_protein(str(protein_file))
        logger.info("Scaffold and protein added to chemical space")

        # Initialize RGroups and Linkers
        rgroups = RGroups()
        linkers = Linkers()

        # Build chemical space
        logger.info(
            f"Building chemical space with {num_linkers} linkers x {num_rgroups} R-groups..."
        )
        for i in range(num_linkers):
            if i % 10 == 0:
                logger.info(f"Processing linker {i}/{num_linkers}")
            for j in range(num_rgroups):
                cs.add_rgroups(linkers.Mol[i], rgroups.Mol[j])

        logger.info(f"Chemical space built with {len(cs)} molecules")

        # Save chemical space
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, "wb") as f:
            cloudpickle.dump(cs, f)
        logger.info(f"Chemical space saved to: {output_path}")

        return True

    except Exception as e:
        logger.error(f"Error creating chemical space: {e}")
        return False

    finally:
        # Cleanup
        try:
            lc.close(timeout=5)
        except Exception as cleanup_error:
            logger.warning(f"Error during cleanup: {cleanup_error}")
            try:
                lc.shutdown()
            except:
                pass


def main():
    parser = argparse.ArgumentParser(
        description="Create chemical space from scaffold and protein"
    )
    parser.add_argument("--scaffold", required=True, help="Input scaffold pickle file")
    parser.add_argument("--protein", required=True, help="Input protein PDB file")
    parser.add_argument(
        "--num-linkers", type=int, default=10, help="Number of linkers to use"
    )
    parser.add_argument(
        "--num-rgroups", type=int, default=10, help="Number of R-groups to use"
    )
    parser.add_argument(
        "--output", default="chemspace.pkl", help="Output chemical space pickle file"
    )

    args = parser.parse_args()

    success = create_chemical_space(
        args.scaffold, args.protein, args.num_linkers, args.num_rgroups, args.output
    )
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
