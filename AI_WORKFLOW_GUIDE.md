# AI-Accelerated Workflow Building Guide

How to use Claude Code to build Silva workflows faster — a companion to the [Silva Migration Guide](./SILVA_MIGRATION_GUIDE.md).

**Authors:** Dev-App Team, Chiral Inc.
**Last Updated:** June 2026

---

## How This Guide Relates to the Migration Guide

The [Silva Migration Guide](./SILVA_MIGRATION_GUIDE.md) defines **what** to build — the structure, conventions, and manual process. This guide defines **how** to use AI (Claude Code) to accelerate each step.

Read the migration guide first. Then use this guide as an overlay.

**Core principle:** AI output quality is proportional to context quality. The most important thing you can do is feed Claude the right context before asking it to generate anything.

### The Development Loop

```
Prepare Context → Scaffold → Node-by-Node → Docker Debug → Review → Report → E2E Test
```

Each section below maps to a step in this loop.

---

## 1. What Context Does AI Need

Before asking Claude to build anything, prepare this checklist of context. Copy-paste or point Claude to these files at the start of your conversation.

### 1.1 Similar Examples from This Repo

Pick 1–2 existing workflows with similar characteristics as references. Feed Claude their `workflow.toml`, one representative `job.toml`, and the README.

| Your tool looks like... | Reference workflow | Key files to share |
|---|---|---|
| Simple pip-installable Python package | workflow-014 (ADMET-AI) | `workflows/workflow-014/.chiral/workflow.toml`, `workflows/workflow-014/Dockerfile` |
| Complex conda environment | workflow-004 (AutoDock Vina) | `workflows/workflow-004/Dockerfile` |
| GPU + large model weights | workflow-012 (Boltz-2) | `workflows/workflow-012/.chiral/workflow.toml`, `workflows/workflow-012/README.md` |
| Package with upstream bugs to patch | workflow-015 (mDeepFRI) | `workflows/workflow-015/Dockerfile` |
| Fan-out comparison (multiple models) | workflow-018 (Boltz-2 vs Chai-1) | `workflows/workflow-018/.chiral/workflow.toml`, `workflows/workflow-018/README.md` |

### 1.2 Workflow Proposal (Node Structure)

Write your proposed node structure as a markdown table before starting implementation. If you have a Google Doc proposal, ask Claude to convert it:

> "Convert this proposal document to a node structure table with columns: Node number, Name, Description, Inputs, Outputs, GPU required. Follow the format used in `workflows/workflow-018/README.md`."

Example node structure:

| Node | Name | Description | Inputs | Outputs | GPU |
|------|------|-------------|--------|---------|-----|
| 00 | Download | Fetch sequence from UniProt | (user params) | `sequence.fasta` | No |
| 01 | Validate | Check input format | `sequence.fasta` | `validated.fasta` | No |
| 02 | Predict | Run structure prediction | `validated.fasta` | `structure.cif`, `confidence.json` | Yes |
| 03 | Report | Generate HTML dashboard | `structure.cif`, `confidence.json` | `report.html` | No |

### 1.3 GitHub Issue Number

Create a GitHub issue **before** you start building (see [Migration Guide Section 8.2](./SILVA_MIGRATION_GUIDE.md#82-step-1--open-a-github-issue)). Give Claude the issue number so it can include `Closes #NNN` in commits and PR descriptions.

### 1.4 The Tool's Original Repo + References

Feed Claude as much of the following as you can find:

- **GitHub URL** and README of the tool
- **Installation instructions** (pip install, conda, etc.)
- **CLI help output** (`tool --help`) or Python API usage example
- **Paper** (if applicable) — at minimum the abstract and methods section
- **Example input/output** — what the tool actually takes and produces

This is what Claude needs to write correct prediction scripts. Without it, Claude will guess the API and get it wrong.

### 1.5 The .chiral Config Format

Share one example of each config file so Claude knows the exact format:

- **workflow.toml**: `workflows/workflow-018/.chiral/workflow.toml`
- **job.toml**: `workflows/workflow-018/00-download/.chiral/job.toml`

Key things Claude needs to know about the format:
- Dependencies are declared centrally in `workflow.toml`, not in individual `job.toml` files
- Parameters use `env = "PARAM_..."` to map to environment variables
- Only nodes with upstream dependencies appear in `[dependencies]` (the first node is omitted)

---

## 2. How to Prepare a Dockerfile

There are 4 common Dockerfile patterns in this repo. Tell Claude which pattern fits your tool.

### Pattern 1: Minimal pip (simplest)

For tools available as a Python package. Example: `workflows/workflow-014/Dockerfile`

```dockerfile
FROM continuumio/miniconda3
ENV USER=silva HOME=/tmp MPLCONFIGDIR=/tmp/matplotlib
RUN pip install --no-cache-dir admet-ai
WORKDIR /workspace
```

**Prompt:**
> "Create a Dockerfile for [tool]. It's installable via `pip install [package]`. Use the same pattern as `workflows/workflow-014/Dockerfile`. Pin the package version."

### Pattern 2: Complex conda environment

For tools with heavy native dependencies. Example: `workflows/workflow-016/Dockerfile`

**Prompt:**
> "Create a Dockerfile for [tool] based on `nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu22.04`. It needs a conda environment with [packages]. Install Miniconda, create the environment from an environment.yml, and set the conda env as default. Follow the pattern in `workflows/workflow-016/Dockerfile`."

### Pattern 3: GPU + model weight caching

For ML models that download large weights at runtime. Pre-download them at build time so prediction runs offline.

**Prompt:**
> "Create a Dockerfile for [tool]. It downloads model weights from HuggingFace on first run — we need to pre-cache them at build time. Use `huggingface-cli download [model]` or the tool's download command. Set `HF_HOME=/models`. Add a smoke test: `CMD python -c 'import [tool]; print(\"OK\")'`."

### Pattern 4: Package patching

For tools with known bugs in installed packages. Example: `workflows/workflow-015/Dockerfile`

**Prompt:**
> "Create a Dockerfile for [tool]. After pip install, we need to patch [file] in the installed package to fix [bug]. Use a Python heredoc or `sed` to apply the fix. See `workflows/workflow-015/Dockerfile` for the pattern."

### Dockerfile Checklist

- [ ] Pin all package versions (`pip install package==1.2.3`)
- [ ] Use `--no-cache-dir` with pip
- [ ] Clean apt caches (`rm -rf /var/lib/apt/lists/*`)
- [ ] Set `WORKDIR /workspace`
- [ ] Pre-download any model weights or data at build time

---

## 3. How to Test a Dockerfile Works

### Build

```bash
docker build -t my-tool:test -f workflows/workflow-XXX/Dockerfile .
```

### Smoke test

```bash
docker run --rm my-tool:test python -c "import my_tool; print('OK')"
```

### GPU test (if applicable)

```bash
docker run --rm --gpus all my-tool:test nvidia-smi
docker run --rm --gpus all my-tool:test python -c "import torch; print(torch.cuda.is_available())"
```

### Debugging build failures with AI

Copy the full error output and paste it to Claude:

> "This Dockerfile build failed with the following error: [paste]. Here is the Dockerfile: [paste]. Fix the issue."

---

## 4. How to Build Node Structure with AI

### Step 1: Scaffold

Ask Claude to generate the entire directory structure in one shot:

> "Scaffold a Silva workflow called 'workflow-XXX' with the following node structure: [paste your table from Section 1.2]. Use `.chiral/workflow.toml` for the DAG and `.chiral/job.toml` per node. Follow the conventions in `workflows/workflow-018/`. Include a README and example input file."

This creates all directories, config files, `run.sh` scripts, and the README in one pass.

### Step 2: Implement node by node

After scaffolding, implement each node as a separate step. **Do not ask Claude to write all nodes at once** — build and test one at a time.

> "Implement node 02-predict for workflow-XXX. It reads `validated.fasta` from `./inputs/` and writes `structure.cif` to `./outputs/`. The Docker image is `my-tool:test`. Here is the tool's Python API: [paste]. Read parameters from `PARAM_*` env vars as defined in the job.toml."

### Step 3: Commit each node separately

Follow the commit pattern used in this repo:

```
feat(workflow-XXX): scaffold [name] workflow
feat(workflow-XXX): add 00-download node
feat(workflow-XXX): add 01-validate-inputs node
feat(workflow-XXX): add 02-predict node
feat(workflow-XXX): add 03-report node
```

---

## 5. How to Debug Each Node with Docker + Sample Input

### Running a single node locally

```bash
# Create test directories
mkdir -p /tmp/test-node/inputs /tmp/test-node/outputs

# Copy test input files
cp workflows/workflow-XXX/input_files/test.yaml /tmp/test-node/inputs/

# Run the node
docker run --rm \
  -v /tmp/test-node/inputs:/workspace/inputs \
  -v /tmp/test-node/outputs:/workspace/outputs \
  -v $(pwd)/workflows/workflow-XXX/02-predict:/workspace \
  -w /workspace \
  -e PARAM_NUM_LOOPS=3 \
  -e PARAM_METHOD=default \
  my-tool:test bash run.sh
```

Key points:
- Mount `inputs/` and `outputs/` as volumes
- Mount the node directory to `/workspace` so `run.sh` and scripts are available
- Set `PARAM_*` environment variables manually (Silva does this automatically in production)
- For GPU nodes, add `--gpus all`

### The debug loop

1. Run the node with Docker
2. If it fails, copy the error output
3. Paste to Claude: "This node failed with: [error]. Here is run.sh: [paste]. Here is the Python script: [paste]. Fix the issue."
4. Apply the fix, re-run
5. Repeat until the node produces correct output files in `/tmp/test-node/outputs/`

### Chaining nodes

To test node B that depends on node A's output:

```bash
# Run node A
docker run --rm \
  -v /tmp/node-a/inputs:/workspace/inputs \
  -v /tmp/node-a/outputs:/workspace/outputs \
  ...

# Feed A's output as B's input
cp /tmp/node-a/outputs/* /tmp/node-b/inputs/

# Run node B
docker run --rm \
  -v /tmp/node-b/inputs:/workspace/inputs \
  -v /tmp/node-b/outputs:/workspace/outputs \
  ...
```

---

## 6. Code Review — Keep It Simple and Working

### Self-review with AI

Before pushing, ask Claude to check your workflow:

> "Review this workflow for correctness. Check:
> 1. All `outputs` declared in each `job.toml` are actually written by the scripts
> 2. All `inputs` declared match what the upstream node produces
> 3. The `[dependencies]` in `workflow.toml` match the actual data flow
> 4. All shell scripts have `set -e`
> 5. All Python scripts create `./outputs/` before writing
> 6. No hardcoded absolute paths
> 7. Parameters are read from `PARAM_*` env vars, not hardcoded"

### Handling PR review feedback

When you receive review comments on your PR:

1. Copy the review comments
2. Paste to Claude: "Here are the PR review comments: [paste]. Apply these fixes to the relevant files."
3. Commit with: `fix(workflow-XXX): address PR review comments`

---

## 7. How to Prepare a Report Template

Most workflows end with a visualization node that generates a self-contained HTML report. The report is a single `.html` file with all libraries loaded from CDN — no local assets needed.

### Recommended libraries

| Library | CDN | Use when |
|---|---|---|
| **Bootstrap 5** | `cdn.jsdelivr.net/npm/bootstrap@5.3.3` | Always — responsive layout |
| **Plotly.js** | `cdn.plot.ly/plotly-2.35.2.min.js` | Charts, heatmaps, radar plots |
| **3Dmol.js** | `3Dmol.org/build/3Dmol-min.js` | Simple protein/molecule 3D viewer (lightweight, easy API) |
| **Mol\*** (Molstar) | `cdn.jsdelivr.net/npm/molstar@latest` | Advanced 3D viewer with toggle, superposition, multi-model support |
| **SMILES-drawer** | `cdn.jsdelivr.net/npm/smiles-drawer` | 2D molecule structure rendering from SMILES strings |

Pick the libraries that fit your workflow's output. Not every report needs a 3D viewer.

### Recommended conventions for new workflows

> **Note:** Existing workflows in this repo vary in how they handle report I/O. The conventions below are the recommended standard for all new workflows.

1. **Read inputs from `./inputs/`** using hardcoded paths (not argparse). This aligns with how Silva populates the inputs directory from upstream node outputs.
2. **Write the report to `./outputs/report.html`**. This follows the same `./outputs/` convention as every other node and matches the `outputs` declaration in `job.toml`.
3. **Read parameters from `PARAM_*` environment variables**, not hardcoded values or argparse flags.
4. **Embed all data inline** as JavaScript variables in the HTML. The report must be fully self-contained.

### Prompt for report generation

> "Generate a Python script `generate_report.py` that reads [output files] from `./inputs/` and writes a self-contained HTML report to `./outputs/report.html`. Use Bootstrap 5 for layout and Plotly.js for charts, both loaded from CDN. Embed all data inline as JavaScript variables. Read inputs from `./inputs/` directly (no argparse). Create `./outputs/` with `os.makedirs`."

### Reference report scripts

- **ADMET dashboard** (radar charts + data tables): `workflows/workflow-014/04_visualize/generate_report.py`
- **Structure comparison** (Mol* 3D viewers + metric plots): `workflows/workflow-018/05-visualize/generate_report.py`
- **Docking results** (3Dmol.js 3D viewer): `workflows/workflow-017/04_visualize/generate_report_3dmol.py`
- **ML prediction** (SMILES-drawer molecules): `workflows/workflow-005/04-prediction/generate_prediction_report.py`

---

## 8. End-to-End Testing

### Phase 1: Docker chain test

Run all nodes sequentially using Docker (as described in [Section 5](#5-how-to-debug-each-node-with-docker--sample-input)). This validates the full file flow without needing Silva installed.

Verify:
- [ ] Each node produces the expected output files
- [ ] Files pass correctly between nodes (output of node N = input of node N+1)
- [ ] The final report.html opens in a browser and looks correct

### Phase 2: Silva test

Once the Docker chain works, test with Silva:

```bash
export SILVA_WORKFLOW_HOME=/path/to/directory/containing/your/workflow
./silva
```

Select your workflow and run. All nodes should complete with green checkmarks.

### Troubleshooting with AI

> "The Silva workflow runs nodes 00 through 02 successfully, but node 03 fails with: [error]. Here is the workflow.toml: [paste]. Here is node 03's job.toml: [paste]. Here is run.sh: [paste]. The inputs directory contains: [ls output]. What is wrong?"

---

## 9. Prompt Templates Quick Reference

| Phase | Prompt |
|-------|--------|
| **Scaffold** | "Scaffold a Silva workflow for [tool] with nodes: [table]. Use `.chiral/` config format. Follow `workflows/workflow-018/` conventions." |
| **Dockerfile** | "Create a Dockerfile for [tool] based on [base image]. Pin versions, use `--no-cache-dir`, set `WORKDIR /workspace`." |
| **Node implementation** | "Implement node [NN]-[name]. It reads [inputs] from `./inputs/` and writes [outputs] to `./outputs/`. Here is the tool's API: [paste]." |
| **Debug** | "This node failed with: [error]. Here is run.sh and the Python script. Fix the issue." |
| **Report** | "Generate `generate_report.py` that reads [files] and writes a self-contained HTML report with Bootstrap 5 + Plotly.js." |
| **Review** | "Review this workflow: check outputs match job.toml, dependencies match data flow, no hardcoded paths." |
| **PR description** | "Write a PR description for workflow-XXX. Closes #[issue]. Include summary, workflow structure, and verification checklist." |
| **Commit message** | "Write a commit message for this change. Use the format `feat(workflow-XXX): ...` or `fix(workflow-XXX): ...`." |

---

## 10. What Else AI Can Help With

- **Generating test input data**: "Generate a small valid [YAML/FASTA/SDF/CSV] test file for [tool]. Include 2–3 entries with realistic data."
- **Writing the GitHub issue**: "Draft a GitHub issue using the template from the migration guide for [tool]. Here is the tool's README: [paste]."
- **Converting existing scripts**: "Refactor this script to use `./inputs/` and `./outputs/` conventions. Read parameters from `PARAM_*` env vars instead of hardcoded values."
- **Troubleshooting dependency conflicts**: "This pip install failed with: [error]. Here is the Dockerfile. Fix the dependency conflict."
- **Updating applications.json**: "Add an entry to `apps/applications.json` for [tool] following the existing format."

---

## Appendix: Reference Workflow Map

| Characteristic | Workflow | Path |
|---|---|---|
| Simple pip tool, linear pipeline | workflow-014 (ADMET-AI) | `workflows/workflow-014/` |
| Pip + package patching | workflow-015 (mDeepFRI) | `workflows/workflow-015/` |
| Complex conda + GPU | workflow-016 (DiffDock-PP) | `workflows/workflow-016/` |
| Fan-out comparison (parallel models) | workflow-018 (Boltz-2 vs Chai-1) | `workflows/workflow-018/` |
| Multi-step ML pipeline | workflow-005 (QSAR) | `workflows/workflow-005/` |
| Protein docking | workflow-017 (LightDock) | `workflows/workflow-017/` |
