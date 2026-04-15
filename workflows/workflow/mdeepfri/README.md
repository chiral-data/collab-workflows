# Create a silva-runnable workflow for mDeepFri from the node structure

- Node structure '~/dev/collab-workflows/workflows/workflow/mdeepfri/'
- Draft workflow 'https://github.com/chiral-data/collab-workflows/issues/90' and Sample Directory Structure (below)
- Workflow references
    - '~/dev/collab-workflows/workflows/workflow-007'
    - '~/dev/collab-workflows/workflows/workflow-014/'
- Silva source code '~/dev/silva' 
- Silva migration guide 'https://github.com/chiral-data/collab-workflows/blob/main/SILVA_MIGRATION_GUIDE.md#5-step-3--implementation' (5. Step 3 — Implementation)
- mDeepFri project: 'https://github.com/bioinf-mcb/Metagenomic-DeepFRI.git' 

## Tasks 
- [x] Investigate the mDeepFRI issue and repo; create a summary of what's happening and how it works.  
- [x] Using existing examples and references, and the sample node structure, create a workflow by building one node at a time. Make it silva runnable. 
- [x] Run 'silva ~/dev/collab-workflows/workflows/workflow/mdeepfri, debug and fix.
- [x] Simplify to 3 nodes: validate inputs, predict (align + predict combined), visualize. Model weights baked into Dockerfile.

## Investigation Summary

**What mDeepFRI does**: Metagenomic-DeepFRI annotates protein sequences with Gene Ontology (GO) terms using a deep learning approach that combines structural and sequence information. Input is a FASTA file; output is a TSV of GO term predictions with confidence scores.

**How it works (pipeline)**:
1. **MMseqs2 alignment**: Queries each input sequence against one or more FoldComp structural databases (e.g., AlphaFold, PDB100, ESMFold) to find similar proteins with known structures.
2. **Contact map generation**: For structural hits, downloads compressed structure entries and derives residue contact maps.
3. **GCN prediction**: Graph Convolutional Network model uses the contact map to make structure-aware GO predictions (more accurate for structurally matched proteins).
4. **CNN prediction**: For proteins with no structural hit, a sequence-only Convolutional Neural Network model is used as a fallback.
5. **Output**: `results.tsv` (protein, network type, mode, GO ID, confidence score) and `alignment_summary.tsv`.

**Key CLI**:
```bash
# Download model weights (baked into Docker image at /opt/mdeepfri-weights)
mDeepFRI get-models -o /opt/mdeepfri-weights -v 1.1

# Run prediction
mDeepFRI predict-function \
  -i sequences.fasta \
  -w /opt/mdeepfri-weights \   # model weights directory
  -d ./foldcomp_db/ \          # optional FoldComp database(s)
  -o ./output/ \
  -p mf -p bp -p cc \          # prediction modes
  --skip-pdb                   # skip auto-download of PDB100 (large)
```

**Prediction modes**: `mf` (Molecular Function), `bp` (Biological Process), `cc` (Cellular Component), `ec` (Enzyme Commission; v1.0 only).

**Installation**: `pip install mdeepfri` (Python 3.11–3.12). Requires ~1GB for model weights. FoldComp databases are optional but improve accuracy (PDB100 is ~50 GB).

**Why 3 nodes**: The mDeepFRI `predict-function` command internally runs alignment + prediction in one call, so they are kept together in node 02. Model weights are downloaded once at image build time (Dockerfile) rather than at runtime, keeping node 02 to a single responsibility.

---

**Sample Directory Structure**
```
workflows/workflow/mdeepfri/
├── .chiral/workflow.toml
├── Dockerfile
├── input_files/sample_proteins.fasta
├── 01_validate_inputs/
│   ├── .chiral/job.toml
│   ├── .chiral/test_inputs/sample_proteins.fasta
│   ├── run.sh
│   └── validate.py
├── 02_predict/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── predict.py
├── 03_visualize/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── generate_report.py
└── README.md
```

## Outputs

- `report.html` — self-contained HTML dashboard with GO term predictions per protein, score distribution chart, and alignment metadata
- `results.tsv` — raw mDeepFRI output: protein ID, network (gcn/cnn), mode, GO ID, confidence score
- `alignment_summary.tsv` — per-query structural alignment statistics

## Build the Docker Image

```bash
cd ~/dev/collab-workflows/workflows/workflow/mdeepfri
docker build -t mdeepfri:latest .
```

The image installs `mdeepfri` via pip (Python 3.12) and downloads model weights (~1 GB) to `/opt/mdeepfri-weights` at build time.
