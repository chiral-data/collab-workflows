# Create a silva-runnable workflow for ADMET-AI Prediction from the node structure

- Node structure '~/dev/collab-workflows/workflows/workflow-014/'
- Draft workflow 'https://github.com/chiral-data/collab-workflows/issues/81' and Sample Directory Structure (below)
- Workflow references
    - ‘~/dev/collab-workflows/workflows/workflow-007’
    - '~/dev/collab-workflows/workflows/workflow-005/'
    - '~/dev/collab-workflows/workflows/workflow-011/'
- Silva source code ‘~/dev/silva’ 
- Silva migration guide 'https://github.com/chiral-data/collab-workflows/blob/main/SILVA_MIGRATION_GUIDE.md#5-step-3--implementation' (5. Step 3 — Implementation)
- ADMET-AI project: ‘https://github.com/swansonk14/admet_ai' 

## Tasks 
- [] Investigate the ADMET Prediction Pipeline (Training Assignment) issue and ADMET-AI repo; create a summary.
- [] Using existing examples and references, and the sample node structure, create a workflow by building one node at a time. Make it silva runnable. 
- [] Run ‘silva ~/dev/workflow-sc-rna’, debug and fix.

**Sample Directory Structure**
```
workflows/workflow-014/
├── .chiral/workflow.toml
├── input_files/sample_molecules.csv
├── 01_validate_inputs/
│   ├── .chiral/job.toml
│   ├── .chiral/test_inputs/sample_molecules.csv
│   ├── run.sh
│   └── validate.py
├── 02_compute/
│   ├── .chiral/job.toml
│   ├── pre_run.sh          # pip install admet-ai chemprop torch
│   ├── run.sh
│   └── compute_admet.py
├── 03_analyze/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── analyze.py
├── 04_visualize/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── generate_report.py
└── README.md
```

## Outputs