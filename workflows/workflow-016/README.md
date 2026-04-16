# Create a silva-runnable workflow for LightDock from the references and node structure

- Node structure '~/dev/collab-workflows/workflows/workflow-016'
- Draft workflow 'https://github.com/chiral-data/collab-workflows/issues/90' and Sample Directory Structure (below)
- LightDock workflow reference: https://github.com/chiral-data/container-images-for-potter (should include scripts for running lightdock + a lightdock report)
- Workflow references
    - '~/dev/collab-workflows/workflows/workflow-007'
    - '~/dev/collab-workflows/workflows/workflow-014/'
    - '~/dev/collab-workflows/workflows/workflow-015/'
- Silva source code '~/dev/silva' 
- Silva migration guide 'https://github.com/chiral-data/collab-workflows/blob/main/SILVA_MIGRATION_GUIDE.md#5-step-3--implementation' (5. Step 3 — Implementation)
- LightDock project: 'https://github.com/lightdock/lightdock' 

## Tasks 
- [] Investigate the LightDock issue and repo; create a summary of what's happening and how it works.
- [] Use the existing examples and references, and the sample node structure, to build one node at a time to create a workflow. Make the workflow silva runnable. 
- [] Run 'silva ~/dev/collab-workflows/workflows/workflow-016, debug and fix.

## Sample Node Structure
workflows/workflow-016
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