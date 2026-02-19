# ABB3 Antibody Structure Prediction Workflow

A four-node pipeline that converts paired heavy/light chain FASTA files
into predicted 3D antibody structures and an interactive HTML report —
mirroring the steps in `example.ipynb`.

---

## Directory Layout

```
workflow/
├── data/
│   ├── heavy.fasta          ← your heavy-chain sequences
│   └── light.fasta          ← your light-chain sequences
├── node1/
│   └── node1.py             ← read, validate, pair chains
├── node2/
│   └── node2.py             ← ProtT5 PLM embeddings (optional)
├── node3/
│   └── node3.py             ← ABB3 structure prediction
├── node4/
│   └── node4.py             ← HTML visualization report
├── run1.sh
├── run2.sh
├── run3.sh
├── run4.sh
├── run_all.sh               ← runs all nodes end-to-end
└── results/
    ├── node1/               ← paired .pt files
    ├── node2/               ← .pt files + PLM embeddings
    ├── node3/               ← predicted .pdb files
    └── node4/
        ├── report.html      ← interactive viewer (open in browser)
        └── *.pdb            ← copies of predicted structures
```

---

## Node Overview

| Node | Script | Input | Output | Purpose |
|------|--------|-------|--------|---------|
| 1 | `node1/node1.py` | FASTA files | `results/node1/*.pt` | Read, validate, pair chains |
| 2 | `node2/node2.py` | `results/node1/*.pt` | `results/node2/*.pt` | ProtT5 embeddings *(optional)* |
| 3 | `node3/node3.py` | `results/node1/` or `node2/` | `results/node3/*.pdb` | ABB3 inference |
| 4 | `node4/node4.py` | `results/node3/*.pdb` | `results/node4/report.html` | Interactive HTML report |

---

## Quick Start

### 1. Prepare your data

```
data/heavy.fasta   — one entry per antibody (FASTA format)
data/light.fasta   — same number of entries, in matching order
```

Example `data/heavy.fasta`:
```
>Ab1_H
QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSLAISWVRQAPGQGLEWMGG...
>Ab2_H
EVQLVESGGGLVQPGGSLRLSCAASGFTFS...
```

### 2. Run the full workflow (plain ABB3)

```bash
chmod +x run_all.sh run1.sh run2.sh run3.sh run4.sh

HEAVY_FASTA=data/heavy.fasta \
LIGHT_FASTA=data/light.fasta \
CHECKPOINT_PATH=/path/to/best_second_stage.ckpt \
./run_all.sh
```

### 3. Run with ABB3-LM (ProtT5 language model)

```bash
# Option A: compute embeddings on the fly (requires GPU, ~3 GB model)
HEAVY_FASTA=data/heavy.fasta \
LIGHT_FASTA=data/light.fasta \
CHECKPOINT_PATH=/path/to/lm_best_second_stage.ckpt \
USE_PLM=1 \
./run_all.sh

# Option B: use pre-computed embeddings (fast, no GPU needed for Node 2)
HEAVY_FASTA=data/heavy.fasta \
LIGHT_FASTA=data/light.fasta \
CHECKPOINT_PATH=/path/to/lm_best_second_stage.ckpt \
USE_PLM=1 \
PRECOMPUTED_DIR=data/embeddings \
./run_all.sh
```

### 4. View results

Open `results/node4/report.html` in any modern browser.
Each predicted structure gets an interactive py3Dmol panel.

---

## Running Nodes Individually

Each node can be run on its own via its `runN.sh` wrapper.
All options are passed as environment variables.

```bash
# Node 1 only
HEAVY_FASTA=data/heavy.fasta LIGHT_FASTA=data/light.fasta ./run1.sh

# Node 3 only (pointing at Node 1 output)
INPUTS_DIR=results/node1 CHECKPOINT_PATH=/path/to/model.ckpt ./run3.sh

# Node 4 only
INPUTS_DIR=results/node3 REPORT_TITLE="My run" ./run4.sh
```

---

## Environment Variables Reference

| Variable | Default | Description |
|----------|---------|-------------|
| `HEAVY_FASTA` | `data/heavy.fasta` | Path to heavy-chain FASTA |
| `LIGHT_FASTA` | `data/light.fasta` | Path to light-chain FASTA |
| `CHECKPOINT_PATH` | *(required)* | ABB3 model checkpoint `.ckpt` |
| `USE_PLM` | `0` | Set to `1` to run Node 2 (ABB3-LM) |
| `PRECOMPUTED_DIR` | *(unset)* | Directory with pre-computed PLM `.pt` files |
| `DEVICE` | `auto` | `cpu`, `cuda`, or `auto` |
| `REPORT_TITLE` | `ABB3 Structure Predictions` | HTML report title |
| `RESULTS_ROOT` | `results` | Root directory for all outputs |

---

## Requirements

- Python ≥ 3.9
- `abodybuilder3` (and its dependencies, including PyTorch)
- `py3Dmol` (only needed for the notebook; Node 4 uses the CDN version)
- For ABB3-LM: `transformers`, `sentencepiece` (ProtT5 deps)

Install:
```bash
pip install abodybuilder3 py3Dmol
```