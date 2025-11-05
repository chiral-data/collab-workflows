🧬 Workflow: Protein–Ligand Virtual Screening (AutoDock Vina + P2Rank)
Overview

This workflow automates protein–ligand virtual screening using AutoDock Vina, with optional binding pocket prediction via P2Rank.
It is divided into three modular jobs that can be run independently or sequentially through Silva, ensuring clarity, reproducibility, and flexibility.

The workflow includes:

1. 01_download – retrieves receptor and ligand structures.

2. 02_prepare – prepares receptor and ligand files for docking.

3. 03_virtual_screening – performs docking using AutoDock Vina and ranks ligands automatically by binding affinity.

Each job runs in a separate container with its own dependencies to keep the environment clean and efficient.

🧩 Job Structure
01_download

Purpose:
Downloads protein and ligand structures from public databases (RCSB, PubChem, etc.).

Main Files:

- job.toml – defines the container and commands.

- Dockerfile – installs Python and libraries (requests, biopython).

- download_job.py – script for downloading receptor and ligand files.

Inputs:

Protein ID (e.g., PDB code).

List of ligand IDs.

Outputs:

receptor.pdb

ligands/ (folder of ligand files)

02_prepare

Purpose:
Prepares the receptor and ligand structures for docking (protonation, cleaning, and conversion to .pdbqt).

Main Files:

Dockerfile – installs MGLTools and Open Babel.

prepare_protein.py – prepares protein.

prepare_ligands.py – prepares ligands.
Inputs:

receptor.pdb

ligands/

Outputs:

receptor.pdbqt

ligands_prepared/

03_virtual_screening

Purpose:
Performs docking using AutoDock Vina and optionally identifies binding pockets using P2Rank.
It then parses Vina log files to generate an Excel sheet ranking ligands by binding affinity.

Main Files:

Dockerfile – installs AutoDock Vina, Python (pandas, openpyxl), and optionally Java for P2Rank.

run_vina.py – performs docking and saves log files.

rank_vina.py – ranks ligands automatically and saves results as results.xlsx.

p2rank_helper.py (optional) – allows selecting different binding pockets for docking.

Inputs:

receptor.pdbqt

ligands_prepared/

Outputs:

vina_results/ – docking poses and logs.

results.xlsx – ranked binding affinities.

🧠 How to Run

Clone the repository:

git clone https://github.com/yourusername/vina_workflow.git
cd vina_workflow


Set your Silva home directory:

export SILVA_HOME_DIR=/path/to/silva_home


Run each job in sequence using Silva:

 01_download
 02_prepare
 03_virtual_screening


The ranked results will appear in 03_virtual_screening/results.xlsx.

⚙️ Dependencies

Each job uses a minimal environment, installing only what’s needed:

Python 3 – scripting, automation, and data handling

Requests / Biopython – downloading and handling biological data

Open Babel / MGLTools – molecular file preparation

AutoDock Vina – docking engine

Pandas / OpenPyXL – ranking and Excel export

Java + P2Rank (optional) – pocket prediction

📚 References

AutoDock Vina — https://github.com/ccsb-scripps/AutoDock-Vina

Open Babel — https://github.com/openbabel/openbabel

MGLTools — https://ccsb.scripps.edu/mgltools/

P2Rank — https://github.com/rdk/p2rank

Biopython — https://biopython.org/

Pandas — https://pandas.pydata.org/

OpenPyXL — https://openpyxl.readthedocs.io/

Silva Workflow Platform — https://github.com/chiral-data/silva