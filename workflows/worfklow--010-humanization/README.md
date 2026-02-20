# workflow-010-antibody-structure-prediction

Predicts Fv antibody structure from paired heavy/light chain FASTA sequences using [ABodyBuilder3 (ABB3)](https://github.com/oxpig/ABodyBuilder3).

Supports both plain **ABB3** and the language-model variant **ABB3-LM** (ProtT5 embeddings).

## Node Overview

| Node | Directory | Description |
|------|-----------|-------------|
| 01 | `01-fasta-validation/` | Read, validate, and pair heavy/light FASTA sequences |
| 02 | `02-plm-embedding/` | *(Optional)* Generate ProtT5 PLM embeddings for ABB3-LM |
| 03 | `03-structure-prediction/` | Run ABB3 forward pass; output one PDB per pair |
| 04 | `04-visualization-report/` | Generate self-contained interactive HTML report |

## Data Flow

```
inputs/heavy.fasta  ─┐
inputs/light.fasta  ─┘→ [01] → outputs/*.pt
                                     │
                              (optional) [02] → outputs/*.pt (+ PLM embedding)
                                                      │
                              inputs/checkpoint.ckpt ─┤
                                                   [03] → outputs/*.pdb
                                                               │
                                                           [04] → outputs/report.html
```

## Running Locally (Docker)

```bash
# Build
docker build -t chiral/workflow-010-antibody-structure-prediction:latest .

# Node 01
docker run --rm \
  -v $(pwd)/data:/workflow/01-fasta-validation/inputs \
  -v $(pwd)/results/01:/workflow/01-fasta-validation/outputs \
  chiral/workflow-010-antibody-structure-prediction:latest \
  bash 01-fasta-validation/run.sh

# Node 03 (skip 02 for plain ABB3)
docker run --rm \
  -e CHECKPOINT_PATH=/workflow/03-structure-prediction/inputs/checkpoint.ckpt \
  -v $(pwd)/results/01:/workflow/03-structure-prediction/inputs \
  -v $(pwd)/results/03:/workflow/03-structure-prediction/outputs \
  chiral/workflow-010-antibody-structure-prediction:latest \
  bash 03-structure-prediction/run.sh
```

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `HEAVY_FASTA` | `inputs/heavy.fasta` | Path to heavy chain FASTA |
| `LIGHT_FASTA` | `inputs/light.fasta` | Path to light chain FASTA |
| `CHECKPOINT_PATH` | `inputs/checkpoint.ckpt` | Path to ABB3 model checkpoint |
| `USE_PLM` | `0` | Set `1` to run Node 02 (ABB3-LM) |
| `DEVICE` | `cpu` | `cpu` or `cuda` |
| `REPORT_TITLE` | `ABB3 Structure Predictions` | Title for HTML report |