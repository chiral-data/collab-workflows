#!/usr/bin/env python3
"""
Run active learning cycles with molecular docking.
"""

import argparse
import logging
import os
import sys
from pathlib import Path

import cloudpickle
import fegrow
from dask.distributed import LocalCluster
from fegrow.al import Model, Query
from rdkit import Chem

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("active_learning.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)


def setup_gnina():
    """Setup gnina binary path."""
    logger.info("Setting up gnina...")

    # First, try system-installed gnina
    system_gnina = Path("/usr/local/bin/gnina")
    if system_gnina.exists():
        fegrow.RMol.set_gnina(str(system_gnina))
        logger.info(f"Using system-installed gnina: {system_gnina}")
        return

    # Fall back to local gnina for backward compatibility
    local_gnina = Path("gnina")
    if local_gnina.exists():
        fegrow.RMol.set_gnina(str(local_gnina))
        logger.info(f"Using local gnina: {local_gnina.absolute()}")
        return

    # If neither found, log warning
    logger.warning("gnina not found in system path or local directory")
    logger.warning("Searched: /usr/local/bin/gnina and ./gnina")
    for f in os.listdir("."):
        if "gnina" in f.lower():
            logger.info(f"Found gnina-like file: {f}")


def setup_active_learning(cs, model_type: str, query_type: str):
    """
    Setup active learning model and query strategy.

    Args:
        cs: Chemical space object
        model_type: Type of model to use
        query_type: Type of query strategy to use
    """
    logger.info(
        f"Setting up active learning with {model_type} model and {query_type} query..."
    )

    # Set model
    if model_type == "gaussian_process":
        cs.model = Model.gaussian_process()
    elif model_type == "random_forest":
        cs.model = Model.random_forest()
    elif model_type == "linear":
        cs.model = Model.linear()
    else:
        raise ValueError(f"Unknown model type: {model_type}")

    # Set query strategy
    if query_type == "UCB":
        cs.query = Query.UCB(beta=1)
    elif query_type == "Greedy":
        cs.query = Query.Greedy()
    elif query_type == "PI":
        cs.query = Query.PI(tradeoff=0)
    elif query_type == "EI":
        cs.query = Query.EI(tradeoff=0)
    else:
        raise ValueError(f"Unknown query type: {query_type}")

    logger.info("Active learning setup complete")


def run_active_learning(
    chemspace_path: str,
    initial_molecules: int,
    num_cycles: int,
    molecules_per_cycle: int,
    model_type: str,
    query_type: str,
    output_path: str,
) -> bool:
    """
    Run active learning docking cycles.

    Args:
        chemspace_path: Path to chemical space pickle file
        initial_molecules: Number of initial random molecules
        num_cycles: Number of active learning cycles
        molecules_per_cycle: Molecules per cycle
        model_type: Type of surrogate model
        query_type: Active learning query strategy
        output_path: Path to output SDF file

    Returns:
        True if successful, False otherwise
    """
    logger.info("Starting active learning docking...")

    # Load chemical space
    chemspace_file = Path(chemspace_path)
    if not chemspace_file.exists():
        logger.error(f"Chemical space file not found: {chemspace_path}")
        return False

    with open(chemspace_file, "rb") as f:
        cs = cloudpickle.load(f)
    logger.info(f"Loaded chemical space with {len(cs)} molecules")

    # Debug: Check if protein/receptor is attached
    if hasattr(cs, "receptor") and cs.receptor is not None:
        logger.info(f"Receptor attached: {type(cs.receptor)}")
    else:
        logger.warning("No receptor attached to ChemSpace - will reload from file")
        # Reload protein from rec_final.pdb (passed as input to this node)
        protein_file = Path("inputs/rec_final.pdb")
        if protein_file.exists():
            cs.add_protein(str(protein_file))
            logger.info(f"Reloaded protein from: {protein_file}")
        else:
            logger.error("rec_final.pdb not found - cannot dock without receptor!")
            return False
    if hasattr(cs, "scaffold") and cs.scaffold is not None:
        logger.info(f"Scaffold attached: {type(cs.scaffold)}")

    # Setup Dask cluster for evaluation
    n_workers = 4  # Use 4 cores
    logger.info(
        f"Setting up Dask LocalCluster with {n_workers} workers for evaluation..."
    )
    eval_cluster = LocalCluster(
        processes=True,
        n_workers=n_workers,
        threads_per_worker=1,
        memory_limit="4GB",
        dashboard_address=":0",  # Use random port to avoid conflicts
        silence_logs=False,
    )

    try:
        # Attach cluster to ChemSpace
        cs.dask_cluster = eval_cluster
        logger.info("Dask cluster attached to ChemSpace")

        # Setup gnina
        setup_gnina()

        # Initial random selection
        logger.info(
            f"Performing initial random selection of {initial_molecules} molecules..."
        )
        random_molecules = cs.active_learning(initial_molecules, first_random=True)

        # Evaluate initial selection
        logger.info("Evaluating initial selection...")
        cs.evaluate(
            random_molecules,
            num_conf=5,
            gnina_gpu=False,
            penalty=0.0,
            al_ignore_penalty=False,
        )

        computed = cs.df[~cs.df.score.isna()]
        logger.info(f"Computed cases in total: {len(computed)}")

        # Setup active learning
        setup_active_learning(cs, model_type, query_type)

        # Run active learning cycles
        logger.info(f"Starting {num_cycles} active learning cycles...")
        for cycle in range(num_cycles):
            logger.info(f"Starting cycle {cycle + 1}/{num_cycles}...")
            try:
                # Check available molecules
                available_molecules = len(cs.df[~cs.df.score.isna()])
                if available_molecules < molecules_per_cycle:
                    logger.warning(
                        f"Not enough available molecules ({available_molecules}). "
                        f"Reducing molecules_per_cycle to {available_molecules}"
                    )
                    molecules_per_cycle = max(1, available_molecules)

                picks = cs.active_learning(molecules_per_cycle)
                logger.info(f"Selected {len(picks)} molecules for evaluation")

                picks_results = cs.evaluate(
                    picks,
                    num_conf=10,
                    gnina_gpu=False,
                    penalty=0.0,
                    al_ignore_penalty=False,
                )

                # Save cycle results
                results_file = f"iteration_{cycle}_results.csv"
                picks_results.to_csv(results_file)
                logger.info(
                    f"Cycle {cycle + 1} completed. Results saved to {results_file}"
                )

            except Exception as e:
                logger.error(f"Error in cycle {cycle + 1}: {str(e)}")
                logger.info("Continuing with next cycle...")
                continue

        # Save chemical space to SDF
        logger.info(f"Saving evaluated chemical space to {output_path}...")
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with Chem.SDWriter(str(output_file)) as SD:
            columns = cs.df.columns.to_list()
            columns.remove("Mol")

            for i, row in cs.df.iterrows():
                if row.Success is False:
                    continue

                mol = row.Mol
                mol.SetIntProp("index", i)
                for column in columns:
                    value = getattr(row, column)
                    mol.SetProp(column, str(value))

                mol.ClearProp("attachement_point")
                SD.write(mol)

        logger.info(f"Evaluated chemical space saved to {output_path}")
        return True

    finally:
        # Cleanup Dask cluster
        logger.info("Closing Dask cluster...")
        try:
            eval_cluster.close(timeout=5)
            logger.info("Dask cluster closed successfully")
        except Exception as cleanup_error:
            logger.warning(f"Error during cluster cleanup: {cleanup_error}")
            try:
                eval_cluster.shutdown()
            except:
                pass


def main():
    parser = argparse.ArgumentParser(description="Run active learning docking cycles")
    parser.add_argument(
        "--chemspace", required=True, help="Input chemical space pickle file"
    )
    parser.add_argument(
        "--initial-molecules",
        type=int,
        default=10,
        help="Number of initial random molecules",
    )
    parser.add_argument(
        "--num-cycles", type=int, default=3, help="Number of active learning cycles"
    )
    parser.add_argument(
        "--molecules-per-cycle",
        type=int,
        default=50,
        help="Molecules to evaluate per cycle",
    )
    parser.add_argument(
        "--model-type",
        default="gaussian_process",
        choices=["gaussian_process", "random_forest", "linear"],
        help="Surrogate model type",
    )
    parser.add_argument(
        "--query-type",
        default="UCB",
        choices=["UCB", "Greedy", "PI", "EI"],
        help="Active learning query strategy",
    )
    parser.add_argument(
        "--output",
        default="chemspace_evaluated.sdf",
        help="Output evaluated chemical space SDF file",
    )

    args = parser.parse_args()

    success = run_active_learning(
        args.chemspace,
        args.initial_molecules,
        args.num_cycles,
        args.molecules_per_cycle,
        args.model_type,
        args.query_type,
        args.output,
    )
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
