# Workflow 024: ESMFold2 Structure Prediction

Predict all-atom 3D structures from protein, RNA, and DNA sequences using [ESMFold2](https://huggingface.co/biohub/ESMFold2) (local inference, no API required).

## Pipeline

```
00-download → 01-validate-inputs → 02-predict → 03-report
```

| Node | Description | GPU |
|------|-------------|-----|
| 00-download | Download protein sequence from UniProt by accession ID | No |
| 01-validate-inputs | Validate YAML: check required fields and sequence characters | No |
| 02-predict | Run ESMFold2 local inference; output structure + confidence metrics | Yes |
| 03-report | Generate HTML dashboard with 3D viewer and PAE heatmap | No |

## Inputs

Upload a YAML file describing the chains to fold:

```yaml
sequences:
  - id: A
    type: protein
    sequence: "MSKGEELFTGVV..."
  - id: B
    type: rna
    sequence: "CGACACCUGAUUCC"
  - id: C
    type: dna
    sequence: "GGAATCAGGTGTCG"
```

Supported types: `protein`, `rna`, `dna`, `ligand`

Alternatively, provide a UniProt accession ID via `PARAM_UNIPROT_ID` and let 00-download fetch the sequence.

## Outputs

| File | Description |
|------|-------------|
| `structure.cif` | All-atom structure in mmCIF format |
| `confidence.json` | pLDDT (per-residue + mean), pTM, ipTM, PAE matrix |
| `report.html` | Interactive dashboard with 3Dmol.js viewer and PAE heatmap |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `PARAM_UNIPROT_ID` | `P69905` | UniProt ID (00-download only) |
| `PARAM_NUM_LOOPS` | `3` | Refinement iterations (3–5) |
| `PARAM_NUM_SAMPLING_STEPS` | `50` | Diffusion steps (32–400) |
| `PARAM_NUM_DIFFUSION_SAMPLES` | `1` | Structure samples per run |

## Test Input

`input_files/example.yaml` contains an antibody heavy + light chain (PDB 6YIO).

## Citation

If you use ESMFold2 in your work, please cite:

```bibtex
@misc{candido2026language,
  title  = {Language Modeling Materializes a World Model of Protein Biology},
  author = {Candido, Salvatore and Hayes, Thomas and Derry, Alexander and Rao, Roshan
            and Lin, Zeming and Verkuil, Robert and Wu, Bryan and Lee, Jin Sub
            and others},
  year   = {2026},
  url    = {https://biohub.ai/papers/esm_protein.pdf},
  note   = {Preprint}
}
```
