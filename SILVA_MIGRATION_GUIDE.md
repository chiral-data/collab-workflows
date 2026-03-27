# Silva Workflow Migration Guide

A step-by-step guide for transforming open-source computational biology tools into Silva-runnable workflows for the Potter platform.

**Authors:** Dev-App Team, Chiral Inc.
**Last Updated:** March 2026

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Key Concepts](#2-key-concepts)
3. [Step 1 — Tool Discovery & Selection](#3-step-1--tool-discovery--selection)
4. [Step 2 — Tool Analysis & Node Design](#4-step-2--tool-analysis--node-design)
5. [Step 3 — Implementation](#5-step-3--implementation)
6. [Step 4 — Testing & Verification](#6-step-4--testing--verification)
7. [Step 5 — Documentation & Submission](#7-step-5--documentation--submission)
8. [GitHub Workflow — Issues, Branches & Pull Requests](#8-github-workflow--issues-branches--pull-requests)
9. [Appendix — Common Pitfalls & Tips](#9-appendix--common-pitfalls--tips)

---

## 1. Introduction

### What is Silva?

Silva is Chiral's workflow execution system. It runs multi-step computational pipelines inside isolated Docker containers. Each pipeline step is defined as a **numbered job folder** containing a configuration file (`@job.toml`) and shell scripts. Silva manages execution order, dependency resolution, and file passing between steps.

- **GitHub:** https://github.com/chiral-data/silva
- **Reference workflows:** https://github.com/chiral-data/collab-workflows

> **Tip:** Before building your first workflow, browse the existing workflows in the `collab-workflows` repo. They serve as concrete, runnable examples of everything described in this guide.

### What is Potter?

Potter is Chiral's no-code drag-and-drop platform for molecular simulation and drug discovery workflows. Once a workflow runs successfully in Silva, it can be migrated to Potter using our automated migration tool, making it accessible to non-technical users through a visual interface.

### Your Mission

Take a popular open-source computational biology tool (from GitHub, a paper, etc.) and transform it into a properly structured, end-to-end runnable Silva workflow. The pipeline looks like this:

```
Discovery → Analysis → Silva Workflow → Verification → Potter Migration
```

---

## 2. Key Concepts

### Workflow Structure

A Silva workflow is a directory containing numbered job folders. Each job folder has a `@job.toml` config and associated scripts.

```
my-workflow/
├── input_files/            ← User drops input files here
├── global_params.json      ← Shared parameters across all jobs
├── 01_first_step/
│   ├── @job.toml
│   ├── run.sh
│   └── outputs/
├── 02_second_step/
│   ├── @job.toml
│   ├── params.json
│   ├── pre_run.sh
│   ├── run.sh
│   ├── my_script.py
│   └── outputs/
├── 03_third_step/
│   └── ...
└── README.md
```

### @job.toml — The Job Configuration File

Every job folder must contain a `@job.toml` file. This declares what the job needs and produces.

```toml
# Dependencies: which previous jobs must complete first
depends_on = ["01_first_step"]

# What files this job reads from previous jobs' outputs
inputs = ["data.csv"]

# What files this job produces
outputs = ["result.csv", "plot.png"]

# Docker image to run inside
[container]
docker_image = "python:3.11-slim"

# Scripts to execute (in order: pre → run → post)
[scripts]
pre = "pre_run.sh"    # Optional: install dependencies
run = "run.sh"         # Required: main execution
# post = "post_run.sh" # Optional: cleanup
```

**Key fields explained:**

| Field | Required | Description |
|-------|----------|-------------|
| `depends_on` | No | List of job folder names that must finish before this one starts |
| `inputs` | No | List of filenames this job needs. Silva copies them from dependency outputs into `./inputs/` |
| `outputs` | Yes | List of filenames this job produces. Must be written to `./outputs/` |
| `[container].docker_image` | Yes | Docker image name and tag |
| `[scripts].pre` | No | Runs before `run`. Typically installs pip packages |
| `[scripts].run` | Yes | Main execution script |
| `[scripts].post` | No | Runs after `run`. Cleanup or post-processing |

### File Passing Between Jobs

Silva uses a simple convention for passing data between jobs:

- Each job writes its output files to `./outputs/`
- Downstream jobs declare `depends_on` and `inputs` in their `@job.toml`
- Silva automatically copies the matching output files into `./inputs/` of the dependent job

```
01_prepare/outputs/data.csv  →  02_process/inputs/data.csv
```

### Parameters: global_params.json vs params.json

- **`global_params.json`** — Placed in the workflow root directory. Contains parameters shared across all jobs (e.g., input filenames, atom selection strings). Every job can read this from `./inputs/global_params.json` or the workflow root.
- **`params.json`** — Placed inside a specific job folder. Contains parameters specific to that job (e.g., reference frame index, figure dimensions).

Both are read by your Python scripts at runtime using `json.load()`.

---

## 3. Step 1 — Tool Discovery & Selection

### Where to Find Tools

- **GitHub:** Search for repositories with keywords relevant to your domain (e.g., "ADMET prediction", "protein structure", "molecular docking")
- **Papers with Code:** https://paperswithcode.com — Find tools associated with published papers
- **BioTools Registry:** https://bio.tools — Curated bioinformatics tool registry
- **Awesome lists:** Search for "awesome-computational-biology", "awesome-cheminformatics", etc.
- **Twitter/X, Reddit r/bioinformatics, r/cheminformatics** — Community buzz about new tools

### Selection Criteria Checklist

Before committing to a tool, evaluate it against these criteria:

| Criterion | ✅ Good Sign | ❌ Red Flag |
|-----------|-------------|------------|
| **Installable** | `pip install`, conda, or simple `git clone` | Complex build system, many system-level deps |
| **Docker-friendly** | Runs on Linux, no GUI required | Requires graphical interface, Windows-only |
| **Clear I/O** | Takes files in, produces files out | Interactive-only, no CLI or API |
| **Documented** | README with usage examples, clear API | No docs, unclear what it does |
| **Maintained** | Recent commits, responsive issues | Abandoned (no activity in 2+ years) |
| **Scientifically useful** | Cited in papers, solves a real problem | Toy project, no scientific application |
| **License** | MIT, Apache, BSD, GPL | No license, or restrictive commercial terms |
| **Test data available** | Bundled test files, or small example dataset | Requires proprietary/large data to run |

**Pro tip:** Tools that are already a Python package with a CLI or a clear Python API are the easiest to convert.

### Discovery Documentation

When you find a candidate tool, document the following before proceeding:

1. **Tool name & URL**
2. **What it does** (1-2 sentences)
3. **Input files** (formats, types)
4. **Output files** (formats, types)
5. **Key parameters** the user can configure
6. **Dependencies** (Python packages, system libraries)
7. **Why it's useful** for Potter users

---

## 4. Step 2 — Tool Analysis & Node Design

### How to Decompose a Tool into Nodes

The goal is to break the tool's pipeline into **logically separable steps**, where each step:

- Does one clear thing
- Produces a reusable intermediate output
- Could be skipped or re-run independently

### Node Splitting Principles

1. **Input validation** is always Node 1. Check that user-provided files exist and are valid before any computation.
2. **Separate expensive computation from lightweight post-processing.** If a tool runs a long calculation and then plots results, split them — users shouldn't re-run expensive steps just to change a plot color.
3. **Separate data preparation from analysis.** If a tool needs preprocessing (alignment, filtering, format conversion), make it its own node.
4. **Keep nodes that can be reused.** An alignment step might be useful for multiple downstream analyses.
5. **Don't over-split.** If two operations are tightly coupled and always run together with no intermediate output worth saving, keep them in one node.

### Typical Node Patterns

Most comp bio workflows follow one of these patterns:

**Pattern A — Linear Pipeline (most common):**
```
Validate → Prepare → Compute → Visualize
```

**Pattern B — Fan-out Comparison:**
```
Validate → Prepare → Model A ─┐
                    → Model B ─┼→ Compare → Visualize
                    → Model C ─┘
```

**Pattern C — Iterative Refinement:**
```
Validate → Initial Run → Evaluate → Refine → Final Output
```

### Designing Your Node Structure

For each node, document:

1. **Node name** (e.g., `01_validate_inputs`)
2. **Purpose** (one sentence)
3. **Inputs** (files from upstream nodes or user)
4. **Outputs** (files produced)
5. **Dependencies** (which previous nodes)
6. **Docker image** needed
7. **Python/system packages** needed
8. **Key parameters** (what goes in params.json)

**Write this out as a document and get it reviewed before you start coding.** This is the most important step — getting the node design right saves a lot of rework later.

---

## 5. Step 3 — Implementation

### Step-by-Step Implementation Process

#### 5.1 Create the Folder Structure

```bash
mkdir my-workflow
cd my-workflow
mkdir input_files
mkdir 01_validate_inputs 02_prepare 03_compute 04_visualize
```

#### 5.2 Write global_params.json

Place shared parameters in the workflow root:

```json
{
    "input_file": "input.sdf",
    "output_format": "csv"
}
```

#### 5.3 Implement Node 1 — Input Validation

This node typically uses a lightweight image (like `ubuntu:22.04`) and just checks that files exist.

**`01_validate_inputs/@job.toml`:**
```toml
outputs = ["input.sdf"]

[container]
docker_image = "ubuntu:22.04"

[scripts]
run = "run.sh"
```

**`01_validate_inputs/run.sh`:**
```bash
#!/bin/bash
set -e

echo "Validating input files..."

INPUT="./inputs/input.sdf"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Input file not found at $INPUT"
    exit 1
fi

if [ ! -s "$INPUT" ]; then
    echo "ERROR: Input file is empty"
    exit 1
fi

echo "Input file OK ($(du -h $INPUT | cut -f1))"
cp "$INPUT" ./outputs/input.sdf
echo "Input validation complete."
```

#### 5.4 Implement Computation Nodes

For nodes that run Python code:

**`02_compute/@job.toml`:**
```toml
depends_on = ["01_validate_inputs"]
inputs = ["input.sdf"]
outputs = ["results.csv"]

[container]
docker_image = "python:3.11-slim"

[scripts]
pre = "pre_run.sh"
run = "run.sh"
```

**`02_compute/pre_run.sh`:**
```bash
#!/bin/bash
set -e
pip install some-package another-package
```

**`02_compute/run.sh`:**
```bash
#!/bin/bash
set -e
echo "Running computation..."
python3 compute.py
echo "Done. Results written to outputs/"
```

**`02_compute/params.json`:**
```json
{
    "threshold": 0.5,
    "method": "default"
}
```

#### 5.5 Writing Standalone Python Scripts

Your Python scripts should be self-contained files that:

1. Read parameters from `global_params.json` and/or `params.json`
2. Read input files from `./inputs/`
3. Do the computation
4. Write output files to `./outputs/`

**Template:**

```python
#!/usr/bin/env python3
"""Brief description of what this script does."""

import json
import os

def load_params():
    """Load global and local parameters."""
    global_params = {}
    if os.path.exists("./inputs/global_params.json"):
        with open("./inputs/global_params.json") as f:
            global_params = json.load(f)

    local_params = {}
    if os.path.exists("params.json"):
        with open("params.json") as f:
            local_params = json.load(f)

    return global_params, local_params

def main():
    global_params, local_params = load_params()

    # Read inputs
    input_file = f"./inputs/{global_params.get('input_file', 'input.sdf')}"

    # Do computation
    # ... your code here ...

    # Write outputs
    os.makedirs("./outputs", exist_ok=True)
    # ... write results to ./outputs/ ...

if __name__ == "__main__":
    main()
```

### Docker Image Selection Guide

| Use Case | Recommended Image |
|----------|-------------------|
| File validation only (bash) | `ubuntu:22.04` |
| Python-based computation | `python:3.11-slim` |
| Needs conda packages | `continuumio/miniconda3` |
| Heavy scientific stack (numpy, scipy, etc.) | `python:3.11-slim` + pip install |
| GPU-accelerated computation | `nvidia/cuda:12.x-runtime-ubuntu22.04` |
| Custom/complex environment | Write a `Dockerfile` in the job folder |

**Prefer `python:3.11-slim`** unless you have a specific reason for another image. It's small, fast to pull, and works with pip for most Python packages.

### pre_run.sh Best Practices

- Always start with `set -e` to fail fast on errors
- Install packages with `pip install` — avoid `pip install --upgrade pip` unless needed
- For GitHub-hosted packages: `pip install git+https://github.com/user/repo.git`
- Pin versions when possible: `pip install numpy==1.26.4 pandas==2.2.0`
- Keep it minimal — only install what this specific node needs

---

## 6. Step 4 — Testing & Verification

### Setting Up Silva Locally

1. **Install Silva:**
   ```bash
   git clone https://github.com/chiral-data/silva.git
   cd silva
   cargo build --release
   ```
   (Requires Rust toolchain and Docker)

2. **Set the workflow home directory:**
   ```bash
   export SILVA_WORKFLOW_HOME=/path/to/directory/containing/your/workflow
   ```

3. **Launch Silva:**
   ```bash
   ./target/release/silva
   ```

4. Navigate to the **Workflows** tab, select your workflow, and press **Enter** to run.

### Providing Test Data

Every workflow must include working test input files in the `input_files/` directory so it can run out of the box.

Options for test data:
- **Bundled test files** from the tool's own test suite (e.g., `MDAnalysis.tests.datafiles`)
- **Small example datasets** from the tool's documentation
- **Generated synthetic data** if no test files are available
- **Public datasets** (PDB structures, ChEMBL molecules, etc.)

The test data should be small enough to run in under a few minutes.

### Verification Checklist

Before submitting your workflow, confirm all of the following:

- [ ] All jobs complete with a **green checkmark** in Silva
- [ ] All declared output files exist after the run
- [ ] Output files contain meaningful, correct results (not empty or garbage)
- [ ] The workflow runs from scratch with no modifications (clean `input_files/` → green checkmarks)
- [ ] No hardcoded absolute paths in any scripts
- [ ] All file paths use `./inputs/` and `./outputs/` conventions
- [ ] `global_params.json` and `params.json` values are actually used (not ignored)

### Deliverable Proof

Take a **screenshot** of:
1. The Silva terminal UI showing all jobs completed with green checkmarks
2. A directory listing (`ls -la` or `tree`) showing the output files

---

## 7. Step 5 — Documentation & Submission

### README Template

Every workflow must include a `README.md` at the workflow root. Use this template:

```markdown
# [Workflow Name]

## Overview
Brief description of what this workflow computes and why it's useful.

**Original Tool:** [Link to GitHub repo or paper]

## Pipeline

| Node | Name | Description | Outputs |
|------|------|-------------|---------|
| 01 | Validate Inputs | Check input files | (staged inputs) |
| 02 | ... | ... | ... |
| 03 | ... | ... | ... |
| 04 | ... | ... | ... |

## Input Files
- `file.ext` — Description and format requirements

## Parameters

### global_params.json
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ... | ... | ... | ... |

### Node XX — params.json
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ... | ... | ... | ... |

## Output Files
- `result.csv` — Description
- `plot.png` — Description

## Running
1. Place input files in `input_files/`
2. Set `SILVA_WORKFLOW_HOME` to the parent directory
3. Launch Silva and select this workflow
4. Press Enter to run

## Requirements
- Docker
- Silva (https://github.com/chiral-data/silva)
```

### Submission Process

See the next section for the full GitHub workflow — issues, branches, PRs, and code review.

---

## 8. GitHub Workflow — Issues, Branches & Pull Requests

We work in the **[chiral-data/collab-workflows](https://github.com/chiral-data/collab-workflows)** repository. This section explains the full development workflow from start to finish. If you're new to Git/GitHub collaboration, read this carefully.

### 8.1 Overview of the Flow

```
Open Issue → Create Branch → Do Work → Push → Open PR → Code Review → Merge
```

Every workflow you build follows this cycle. Here's what each step means and how to do it.

### 8.2 Step 1 — Open a GitHub Issue

Before you start coding, create an issue to track the work. This lets the team know what you're working on and provides a place for discussion.

**How to create an issue:**
1. Go to https://github.com/chiral-data/collab-workflows/issues
2. Click **"New issue"**
3. Use this title format: `[Workflow] <Tool Name> — <Brief Description>`
4. In the body, include:

```markdown
## Tool
- **Name:** [tool name]
- **GitHub:** [link to tool repo]
- **Paper:** [link if applicable]

## What it does
[1-2 sentence description]

## Proposed Node Structure
- 01_validate_inputs — Check input files
- 02_... — ...
- 03_... — ...
- 04_... — ...

## Input/Output
- **Inputs:** [file types]
- **Outputs:** [file types]

## Notes
[Any questions, concerns, or things to discuss before starting]
```

5. Assign yourself to the issue
6. **Wait for approval from Rongjun before starting implementation** — he may have feedback on the node design

### 8.3 Step 2 — Clone the Repo & Create a Branch

If you haven't cloned the repo yet:

```bash
git clone https://github.com/chiral-data/collab-workflows.git
cd collab-workflows
```

Create a branch for your work. **Never work directly on `main`.**

```bash
# Make sure you're on the latest main
git checkout main
git pull origin main

# Create a new branch
git checkout -b workflow/<tool-name>
```

**Branch naming convention:** `workflow/<tool-name>` (lowercase, hyphens for spaces)

Examples:
- `workflow/admet-prediction`
- `workflow/boltz-chai-comparison`
- `workflow/ermsf`

### 8.4 Step 3 — Do Your Work

Build your workflow following the implementation steps in Sections 3-7 of this guide. Your work goes under `workflows/` in the repo.

```bash
# Your workflow lives here
workflows/workflow-XXX/
├── input_files/
├── global_params.json
├── 01_.../
├── 02_.../
└── README.md
```

**Commit frequently** as you make progress. Each commit should be a logical chunk of work.

```bash
# Stage your changes
git add .

# Commit with a descriptive message
git commit -m "Add node 01: input validation"

# More work...
git commit -m "Add node 02: alignment step"
git commit -m "Add node 03: main computation"
git commit -m "Add node 04: visualization"
git commit -m "Add README and test data"
```

**Good commit messages:**
- `Add node 01: input validation for ADMET workflow`
- `Fix pre_run.sh: pin numpy version to avoid compatibility issue`
- `Add verification screenshots`

**Bad commit messages:**
- `update`
- `fix`
- `wip`
- `asdfgh`

### 8.5 Step 4 — Push Your Branch

Push your branch to GitHub:

```bash
git push origin workflow/<tool-name>
```

If this is the first push for this branch, Git may ask you to set the upstream:

```bash
git push --set-upstream origin workflow/<tool-name>
```

### 8.6 Step 5 — Open a Pull Request (PR)

Once your workflow is complete and verified (all green checkmarks in Silva), open a PR to merge your work into `main`.

**How to open a PR:**
1. Go to https://github.com/chiral-data/collab-workflows/pulls
2. Click **"New pull request"**
3. Set **base:** `main` ← **compare:** `workflow/<tool-name>`
4. Use this title format: `Add workflow: <Tool Name>`
5. In the PR description, include:

```markdown
## Summary
[Brief description of what this workflow does]

Closes #[issue number]

## Workflow Structure
- 01_... — [description]
- 02_... — [description]
- 03_... — [description]
- 04_... — [description]

## Verification
- [ ] All jobs complete with green checkmarks in Silva
- [ ] All output files present and correct
- [ ] README complete

## Screenshots
[Paste your Silva terminal screenshot showing green checkmarks]
[Paste your output file listing]
```

**Important:** The `Closes #XX` line automatically closes the linked issue when the PR is merged.

6. Assign **Rongjun** as reviewer
7. Wait for review feedback

### 8.7 Step 6 — Code Review & Iteration

Rongjun will review your PR and may leave comments or request changes. This is normal and expected — it's how we maintain quality.

**When you receive review feedback:**

1. Read all comments carefully
2. Make the requested changes in your local branch
3. Commit and push again:

```bash
# Make fixes
git add .
git commit -m "Address review: fix params.json handling in node 03"
git push origin workflow/<tool-name>
```

4. The PR automatically updates with your new commits
5. Reply to each comment to explain what you changed (or discuss if you disagree)
6. Once Rongjun approves, the PR gets merged

### 8.8 Step 7 — After Merge

After your PR is merged, clean up:

```bash
# Switch back to main and pull the latest
git checkout main
git pull origin main

# Delete your local branch (optional, keeps things tidy)
git branch -d workflow/<tool-name>
```

Then start the cycle again for your next workflow.

### 8.9 Common Git Situations

**"I made changes but I'm on the wrong branch"**
```bash
# Stash your changes temporarily
git stash

# Switch to the correct branch (or create it)
git checkout -b workflow/<tool-name>

# Apply your stashed changes
git stash pop
```

**"I need to update my branch with the latest `main`"**
```bash
git checkout main
git pull origin main
git checkout workflow/<tool-name>
git merge main
```

**"I accidentally committed to `main`"**
```bash
# Undo the last commit (keeps your changes)
git reset HEAD~1

# Now create a branch and commit there
git checkout -b workflow/<tool-name>
git add .
git commit -m "Your message"
```

**"I have a merge conflict"**
1. Git will tell you which files have conflicts
2. Open those files — look for `<<<<<<<`, `=======`, `>>>>>>>` markers
3. Edit the file to resolve the conflict (keep the code you want)
4. Remove the conflict markers
5. `git add <file>` and `git commit`
6. Ask for help if you're unsure

### 8.10 Quick Reference

| Action | Command |
|--------|---------|
| Clone repo | `git clone https://github.com/chiral-data/collab-workflows.git` |
| Create branch | `git checkout -b workflow/<name>` |
| Check current branch | `git branch` |
| Stage changes | `git add .` |
| Commit | `git commit -m "message"` |
| Push | `git push origin workflow/<name>` |
| Switch branch | `git checkout <branch-name>` |
| Update from main | `git pull origin main` |
| See status | `git status` |
| See commit history | `git log --oneline` |

---

## 9. Appendix — Common Pitfalls & Tips

### Common Mistakes

1. **Forgetting `set -e` in shell scripts.** Without this, a failed command won't stop execution and you'll get confusing downstream errors.

2. **Hardcoding file paths.** Always use `./inputs/` and `./outputs/` relative paths. Never use absolute paths.

3. **Not creating the `outputs/` directory.** In Python scripts, always include `os.makedirs("./outputs", exist_ok=True)`.

4. **Installing packages in `run.sh` instead of `pre_run.sh`.** Package installation goes in `pre_run.sh` (the `pre` script). Keep `run.sh` for actual computation only.

5. **Declaring outputs that aren't actually written.** If your `@job.toml` says `outputs = ["result.csv"]`, your script must write exactly `./outputs/result.csv`. Mismatches cause silent failures.

6. **Not reading parameters from JSON.** Don't hardcode values that should come from `params.json` or `global_params.json`. The whole point is that users can modify parameters without editing code.

7. **Using `latest` Docker image tags.** Always pin a specific tag (e.g., `python:3.11-slim`, not `python:latest`) for reproducibility.

### Tips for Efficient Development

- **Test your Python scripts locally first** before wrapping them in Silva. Run them in a virtual environment to catch import errors early.
- **Start with the simplest possible version** — get a 2-node workflow running, then add complexity.
- **Look at the tool's test suite** for example usage and test data.
- **Read the tool's source code**, not just the README. The API often supports more than what's documented.
- **Use `python:3.11-slim` as your default image.** Only switch to something else if you have a reason.
- **Keep pre_run.sh minimal.** Every extra package adds to container startup time.

### Useful Resources

- **Silva documentation:** https://github.com/chiral-data/silva
- **Reference workflows:** https://github.com/chiral-data/collab-workflows
- **Existing workflows** (learn by example): https://github.com/chiral-data/collab-workflows/tree/main/workflows
- **Docker Hub:** https://hub.docker.com (for finding base images)
- **Papers with Code:** https://paperswithcode.com (for discovering tools)
