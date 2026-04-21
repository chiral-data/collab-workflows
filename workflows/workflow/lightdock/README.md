# Create a silva-runnable workflow for LightDock from the references and node structure

- Node structure '~/dev/collab-workflows/workflows/workflow-016'
- Draft workflow 'https://github.com/chiral-data/collab-workflows/issues/90' and Sample Directory Structure (below)
- Reference for overall job script: https://github.com/chiral-data/container-images-for-potter/blob/v0.3.0/apps/docking/lightdock/job_script.sh = Job script that covers the full setup → swarm generation → docking → clustering pipeline end-to-end
- Reference for Visualizing and Creating Outputs: https://github.com/chiral-data/container-images-for-potter/blob/v0.3.0/reports/docking_report/lightdock_docking_dashboard.py = Plotly dashboard with CAPRI metrics, cluster ranking, and interface analysis 

- Workflow references
    - '~/dev/collab-workflows/workflows/workflow-007'
    - '~/dev/collab-workflows/workflows/workflow-014/'
    - '~/dev/collab-workflows/workflows/workflow-015/'
- Silva source code '~/dev/silva' 
- Silva migration guide 'https://github.com/chiral-data/collab-workflows/blob/main/SILVA_MIGRATION_GUIDE.md#5-step-3--implementation' (5. Step 3 — Implementation)
- LightDock project: 'https://github.com/lightdock/lightdock' 

## Tasks 
- [] Investigate the LightDock issue and repo.
- [] Use the existing examples and references, and the sample node structure, to build one node at a time to create a workflow. Make the workflow silva runnable. 
- [] Run 'silva ~/dev/collab-workflows/workflows/workflow/lightdock, debug and fix.

## Sample Node Structure
workflows/workflow/lightdock
├── .chiral/workflow.toml
├── Dockerfile
├── input_files/
│   └──protein1_barstar.pdb
│   └──protein2_barnase.pdb
├── 01_validate_inputs/
│   ├── .chiral/job.toml
│   ├── .chiral/test_inputs/
│   │   └──receptor_structure.pdb
│   │   └──ligand_structure.pdb
│   ├── run.sh
│   └── validate.py
├── 02_prepare/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── generate_swarms.py
├── 03_compute/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── run_lightdock.py
├── 04_visualize/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── generate_report.py
└── README.md