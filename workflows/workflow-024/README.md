# Create a silva-runnable workflow for Docking Comparison of Autodock Vina vs gnina from the node structure

Key Files and References: 
- Node structure | '~/dev/collab-workflows/workflows/workflow-024/'
- GitHub Issue | 'https://github.com/chiral-data/collab-workflows/issues/138' 
- Workflow references
    - ‘~/dev/collab-workflows/workflows/workflow-014’
    - '~/dev/collab-workflows/workflows/workflow-018/'
- Silva source code | ‘~/dev/silva’ 
- Silva migration guide | 'https://github.com/chiral-data/collab-workflows/blob/main/SILVA_MIGRATION_GUIDE.md#5-step-3--implementation' (5. Step 3 — Implementation)

Key Tool and Paper References: 
- AutoDock Vina: https://github.com/ccsb-scripps/AutoDock-Vina | pip install vina | Apache 2.0
- GNINA: https://github.com/gnina/gnina | prebuilt binary | Apache 2.0
- Open Babel: https://github.com/openbabel/openbabel | conda install -c conda-forge openbabel | GPL 2.0
- xTB: https://github.com/grimme-lab/xtb | conda install -c conda-forge xtb | LGPL 3.0 
- Buccheri et al., 2025 | 'https://pmc.ncbi.nlm.nih.gov/articles/PMC12388557/'

Docker Files
- AutoDock Vina Dockerfile | 'https://github.com/chiral-data/collab-workflows/blob/main/workflows/workflow-004/Dockerfile'
- Gnina Dockerfile | 'https://github.com/chiral-data/collab-workflows/blob/main/apps/g/gnina_2025_12_04/Dockerfile'

Reference Config File Formats: 
- **workflow.toml**: `workflows/workflow-018/.chiral/workflow.toml`
- **job.toml**: `workflows/workflow-018/00-download/.chiral/job.toml`
- Dependencies are declared centrally in `workflow.toml`, not in individual `job.toml` files
- Parameters use `env = "PARAM_..."` to map to environment variables
- Only nodes with upstream dependencies appear in `[dependencies]` (the first node is omitted)

## Tasks 
- [] Investigate the GitHub issue and fully the provided read references of files, tools, and papers; create a summary of key ideas and goals. 
- [] Rewrite the Sample Directory Structure in the README to match the Proposed Node Structure described in the GitHub Issue, keeping a similar format. 
- [] Using existing examples and references, and the sample directory structure, create a workflow by building one node at a time. Make it silva runnable.

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
│   ├── pre_run.sh          
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