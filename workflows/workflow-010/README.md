# workflow-010-antibody-structure-prediction

Predicts Fv antibody structure from paired heavy/light chain FASTA sequences using [ABodyBuilder3 (ABB3)](https://github.com/Exscientia/abodybuilder3).

Supports both plain **ABB3** and the language-model variant **ABB3-LM** (ProtT5 embeddings).

## Node Overview

| Node | Directory | Description |
|------|-----------|-------------|
| 01 | `01-input-preparation/` | Read, validate, and pair heavy/light FASTA sequences |
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

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `PARAM_DEVICE` | `cpu` | Compute device: `cpu` or `cuda` |
| `PARAM_USE_PLM` | `0` | Set `1` to run Node 02 (ABB3-LM) |
| `PARAM_REPORT_TITLE` | `ABB3 Structure Predictions` | Title for HTML report |

## Running Locally (Docker)

```bash
# Build
docker build -t abodybuilder3:latest .

# Node 01 – Input Preparation
docker run --rm \
  -v $(pwd)/data:/workflow/01-input-preparation/inputs \
  -v $(pwd)/results/01:/workflow/01-input-preparation/outputs \
  -w /workflow/01-input-preparation \
  abodybuilder3:latest \
  bash run.sh

# Node 03 – Structure Prediction (skip 02 for plain ABB3)
docker run --rm \
  -v $(pwd)/results/01:/workflow/03-structure-prediction/inputs \
  -v $(pwd)/results/03:/workflow/03-structure-prediction/outputs \
  -w /workflow/03-structure-prediction \
  abodybuilder3:latest \
  bash run.sh

# Node 04 – Visualization Report
docker run --rm \
  -v $(pwd)/results/03:/workflow/04-visualization-report/inputs \
  -v $(pwd)/results/04:/workflow/04-visualization-report/outputs \
  -w /workflow/04-visualization-report \
  abodybuilder3:latest \
  bash run.sh
```
