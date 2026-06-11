---
doc_id: workflow-NNN            # Unique identifier matching the directory name
domain: <domain>                # e.g.: protein-protein-docking, structure-prediction,
                                #   rna-seq, admet, molecular-docking, metabolomics
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  <One or two sentences describing what this workflow does and what
  inputs it takes.>
tags: [<keyword1>, <keyword2>, <keyword3>]
---

<!-- Workflow README template for the collab-workflows repository.
     Please follow the section structure below when writing or
     updating a workflow README. The ## / ### heading hierarchy
     is important — it powers the AI assistant's documentation
     search on the platform. -->

# <Workflow NNN>: <Workflow Title>

<One paragraph: what this workflow does, what scientific question it
addresses, and what tool/algorithm it uses. Keep it under 4 sentences.>

## Overview

<Expand on the title paragraph. Cover:
- The scientific method or algorithm used
- What type of input data it works with
- What the main output is and how researchers use it
- Any key citations (author, year) for the underlying tool

Target: 1-2 paragraphs.>

## When to use this workflow

<Paragraph 1 — when to use:
Describe the scenario where a researcher should choose this workflow.
Be specific about input requirements (file formats, data types, organism
constraints) and the type of question it answers.>

<Paragraph 2 — when NOT to use + alternatives:
Describe scenarios where this workflow is the wrong choice. Name the
alternative workflow(s) by ID (e.g., "use workflow-XXX instead for...").
Mention common mistakes or misconceptions.>

## Architecture and data flow

<ASCII diagram or Mermaid flowchart showing the pipeline stages and key
intermediate files. Example:

```text
input.pdb ──> [01: Validate] ──> [02: Process] ──> [03: Analyze] ──> [04: Report]
                                      |                  |               |
                                 cleaned.pdb        results.json    report.html
```

One sentence below the diagram: "Nodes run sequentially: 01 → 02 → 03 → 04."
or describe any parallel/conditional branching.>

## Input requirements

<What the user needs to prepare before running. Include:
- Required file format(s) (e.g., PDB, FASTA, CSV with SMILES column)
- Any constraints (organism, file size, structure completeness)
- Where to place files (e.g., `input_files/`)

If sample/test input files are included, mention them here.>

## Workflow nodes

### Node 01: <Node Name>

**Goal:** <One sentence: what this node accomplishes.>

**Process:** <2-4 sentences: what the node actually does — algorithm,
tool invocation, key logic. Be specific enough that a researcher
understands the computation without reading the source code.>

**Scientific notes:** <1-3 sentences: the biological or chemical
rationale. Why does this step exist? What domain knowledge informs
the approach? What assumptions does it make?>

**Outputs:** <Bullet list of output files with one-line descriptions.>

### Node 02: <Node Name>

**Goal:** <...>

**Process:** <...>

**Scientific notes:** <...>

**Outputs:** <...>

<!-- Add a ### section for each node in the workflow -->

## Parameters

### <parameter_name>

<!-- For enum parameters, list each value with decision guidance.
     For numeric parameters, explain the tradeoff and recommended ranges. -->

| Value / Range | Description |
|---------------|-------------|
| `<value1>` (default) | <What it does. When to use it.> |
| `<value2>` | <What it does. When to use it instead of the default.> |

**Trade-off:** <One sentence on the key tradeoff, e.g., speed vs accuracy,
sensitivity vs specificity.>

**Test vs production:** <Default is for testing. For publication-quality
results, set to X.>

### <another_parameter>

- **Type:** <integer / float / string / enum>
- **Default:** <value>
- **Description:** <What this parameter controls.>
- **Guidance:** <When and why a researcher would change it from the default.>

## Outputs and interpretation

### <output_file_or_metric>

<Define what this output is. Include:
- What the values mean
- Typical ranges or thresholds (e.g., "pLDDT > 0.7 indicates high confidence")
- Caveats or limitations (e.g., "does not account for disordered regions")

Target: 2-4 sentences per output/metric.>

### <another_output_or_metric>

<!-- Repeat for each significant output file or quality metric -->

## Quick start

### Running with Docker

<Link to the Dockerfile or list the public Docker image(s) that can
be pulled directly.>

```bash
<Docker run command with default/test parameters>
```

### Running on Silva

1. Select this workflow from the workflow list
2. Upload your input files
3. Adjust parameters if needed (see Parameters section)
4. Click Run

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| <param1> | <test_value> | <prod_value> |
| <param2> | <test_value> | <prod_value> |

<One sentence: what test data is included for validation, and what to
expect from a successful test run.>

## Troubleshooting

<!-- Optional — only include if the workflow has known setup or
     configuration issues worth documenting. Users on the platform
     should contact the dev team for runtime errors. -->

**<Error message or symptom>**
<Cause and fix in 1-2 sentences.>

## References

- <Author> et al. "<Title>." *Journal* Volume(Issue):Pages, Year. DOI: <url>
- [<Tool name> documentation](<url>)
