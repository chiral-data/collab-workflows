# Workflow 2: 3D Structure Prediction and Analysis

## Overview
This workflow demonstrates state-of-the-art AI-based structure prediction using Boltz-2 for Nobel Prize-winning therapeutic breakthroughs that transformed modern medicine.

## Examples

### Example 1: COVID-19 mRNA Vaccine Components - Pseudouridine-Modified RNA
**Nobel Prize Connection (2023):**
- Katalin Karikó and Drew Weissman: nucleoside base modifications for mRNA vaccines
- **Key Discovery**: Pseudouridine (Ψ) and N1-methylpseudouridine (m1Ψ) prevent immune activation
- Foundation technology for COVID-19 vaccines that saved millions of lives

**Scientific Context:**
- Modified mRNA sequences with pseudouridine substitutions
- Spike protein RBD sequences with enhanced stability
- Demonstrates RNA structure prediction with chemical modifications

**Available Sequences:**
- Spike protein RBD sequences from published research
- Well-documented pseudouridine modification patterns
- Structure-function relationships of modified nucleotides

### Example 2: FMC63/OKT3 Bispecific T-Cell Engager - Anti-CD19/CD3 Construct
**Nobel Prize Connection (2018):**
- James Allison & Tasuku Honjo: cancer immunotherapy via checkpoint inhibition
- **Key Discovery**: CTLA-4 and PD-1 pathways for unleashing immune system against cancer
- Led to revolution in cancer treatment including T-cell engager therapies

**Scientific Context:**
- Bispecific scFv construct targeting CD19 (B-cells) and CD3 (T-cells)
- Based on FMC63 (anti-CD19) and OKT3-derived (anti-CD3) sequences
- Foundation technology for CAR-T therapies and bispecific antibodies

**Available Sequences:**
- **FMC63 anti-CD19 scFv**: GenBank ID HM852952.1 (FDA-approved in 4 CAR-T therapies)
- **OKT3-derived anti-CD3 scFv**: Published sequences available
- **Linker design**: Flexible glycine-serine connectors

## Pipeline Structure
```
Nobel Prize Sequences → Boltz-2 Prediction → Structure Analysis → Therapeutic Insights
```

## Directory Structure
```
2_structure_prediction/
├── 1_mRNA/
│   ├── sequences/
│   │   └── P0DTC2.fasta           # Downloaded spike protein  
│   ├── inputs/
│   │   ├── spike_rbd.fasta        # Boltz-2 input file
│   │   ├── job_config.json        # Optimized parameters
│   │   └── job_script.sh          # Executable job script
│   └── results/                   # Prediction outputs
├── 2_antibody/
│   ├── sequences/
│   │   └── FMC63-28Z.fasta        # Downloaded FMC63 protein
│   ├── inputs/
│   │   ├── fmc63.fasta            # Boltz-2 input file  
│   │   ├── job_config.json        # Optimized parameters
│   │   └── job_script.sh          # Executable job script
│   └── results/                   # Prediction outputs
├── create_boltz_inputs.py         # Conversion script
└── README.md
```

## Implementation Status ✅
- [x] Downloaded Nobel Prize-related sequences
- [x] Created Boltz-2 format FASTA files  
- [x] Generated optimized job configurations
- [x] Set up executable job scripts
- [x] Ready for container testing

## Research Significance
These examples showcase:
1. **Two Nobel Prizes** (2018 & 2023) that revolutionized medicine
2. **Accessible sequences** from GenBank and published research
3. **Clinical impact** (COVID-19 vaccines, CAR-T therapies)
4. **Structure-function insights** (RNA modifications, antibody engineering)
5. **Therapeutic applications** actively saving lives worldwide

## Sequence Sources & Download Links

### Example 1: COVID-19 mRNA Components

**Recommended Choice: UniProt P0DTC2**
- **Primary**: [UniProt P0DTC2](https://www.uniprot.org/uniprot/P0DTC2) - SARS-CoV-2 Spike protein (Wuhan reference)
  - Complete spike protein sequence in FASTA format
  - RBD region: amino acids 331-524
  - ✅ **Best option**: Clean sequence, well-annotated

**Alternative Options:**
- [NCBI GenBank NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) - Complete genome
- [PDB 6M0J](https://www.rcsb.org/structure/6M0J) - RBD-ACE2 complex structure
- [PDB 7C2L](https://www.rcsb.org/structure/7C2L) - RBD structure

### Example 2: Bispecific T-Cell Engager

**Recommended Choice: FMC63 + UCHT1**
- **FMC63 Anti-CD19**: [GenBank HM852952.1](https://www.ncbi.nlm.nih.gov/nuccore/HM852952.1)
  - FDA-approved in 4 CAR-T therapies (Kymriah, Yescarta, etc.)
  - Whitlow linker region: AA 130-148
  - ✅ **Best option**: Clinically validated, complete sequence

- **UCHT1 Anti-CD3**: More reliable than OKT3 variants
  - Better for T-cell engager applications
  - Published sequences in CAR-T literature
  - ✅ **Best option**: Superior binding characteristics

**Alternative Options:**
- OKT3-derived anti-CD3 sequences (older, more limitations)
- Other anti-CD19 clones (less clinically validated)

### RNA Modifications
- **Pseudouridine (Ψ)**: PubChem database
- **N1-methylpseudouridine (m1Ψ)**: Chemical structure references
- **Modification patterns**: Well-documented in Karikó/Weissman publications

## Parameter Optimization

### Spike RBD (194 amino acids)
- `recycling_steps`: 5 (moderate complexity)
- `diffusion_samples`: 10 (good sampling coverage)  
- `use_msa_server`: true (leverage evolutionary information)

### FMC63 Antibody (489 amino acids)
- `recycling_steps`: 7 (higher for multi-domain structure)
- `diffusion_samples`: 15 (comprehensive sampling for complex antibody)
- `use_msa_server`: true (critical for antibody folding patterns)

**Rationale:** Larger, more complex proteins require more recycling steps and diffusion samples for accurate structure prediction.

## Next Steps
1. **Container Testing** → Test job scripts with `boltz_potter_nvidia_2` 
2. **Structure Analysis** → Compare predictions with known experimental structures
3. **Validation** → Assess prediction quality and biological relevance
4. **Documentation** → Generate analysis reports for Nobel Prize connections

## Container Information
**Container Image**: `chiral.sakuracr.jp/boltz_dok_nvidia_2`
- Boltz-2 AI structure prediction model
- NVIDIA GPU acceleration support
- Includes MSA server connectivity

## Quick Start
```bash
cd /home/roki/container-images-for-potter/sept_workflows/2_structure_prediction

# Test Spike RBD prediction
cd 1_mRNA/inputs && ./job_script.sh

# Test FMC63 prediction  
cd ../../2_antibody/inputs && ./job_script.sh
```

## Container Usage
```bash
# Pull container (if needed)
singularity pull chiral.sakuracr.jp/boltz_dok_nvidia_2

# Run prediction manually
singularity exec --nv chiral.sakuracr.jp/boltz_dok_nvidia_2 \
    python3 -m boltz.main predict spike_rbd.fasta \
    --use_msa_server --output_format pdb
```

## Test Results

### Test 1: Spike RBD Structure Prediction (2025-09-03)

**Command Used:**
```bash
cd /home/ubuntu/chiral/container-images-for-potter/sept_workflows/2_structure_prediction/1_mRNA/inputs
docker run --rm --gpus all -v $(pwd):/workspace -w /workspace chiral.sakuracr.jp/boltz_dok_nvidia_2:latest ./job_script.sh
```

**Issues Found & Fixes:**
1. **FASTA Header Format**: Initial header `>spike_rbd|protein` caused KeyError. Fixed by changing to `>A|protein|` (Boltz requires single letter chain ID)
2. **CCD Data Re-download**: Container downloads CCD data on each run despite being included in build. This is because `/opt/boltz_cache` doesn't persist between container runs.

**Results:**
- ✅ Successfully generated 10 structure models (model_0 through model_9)
- ✅ Generated confidence scores (JSON files)
- ✅ Generated PAE (Predicted Aligned Error) matrices
- ✅ Generated PDE (Predicted Distance Error) matrices
- ✅ Generated pLDDT (predicted Local Distance Difference Test) scores
- **Total Runtime**: ~4 minutes (including CCD download, MSA generation, and structure prediction)
- **GPU Utilization**: Successfully used CUDA GPU
- **Output Files**: 50 files total in `outputs/` directory