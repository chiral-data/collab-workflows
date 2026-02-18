# Boltz-2 Structure Prediction Dashboard Plan

## Overview
Create a comprehensive analysis dashboard for Boltz-2 structure prediction results, similar to the existing AutoDock Vina dashboard. The dashboard will analyze and visualize the quality metrics of predicted protein structures.

## Boltz-2 Output Analysis

### Available Data from Test Results
From the successful spike RBD prediction, we have 50 output files:

#### 1. Structure Files (10 models)
- `spike_rbd_model_0.pdb` through `spike_rbd_model_9.pdb`
- Contains atomic coordinates for predicted structures

#### 2. Confidence Scores (10 JSON files)
- `confidence_spike_rbd_model_0.json` through `confidence_spike_rbd_model_9.json`
- Contains key metrics:
  - `confidence_score`: Overall model confidence (0-1)
  - `ptm`: Predicted Template Modeling score
  - `iptm`: Interface PTM (for complexes)
  - `complex_plddt`: Overall pLDDT score
  - `complex_pde`: Predicted Distance Error
  - `chains_ptm`: Per-chain PTM scores

#### 3. Quality Metrics (30 NumPy files)
- **PAE (Predicted Aligned Error)**: `pae_spike_rbd_model_*.npz` - Per-residue pair confidence
- **PDE (Predicted Distance Error)**: `pde_spike_rbd_model_*.npz` - Distance prediction confidence
- **pLDDT (predicted LDDT)**: `plddt_spike_rbd_model_*.npz` - Per-residue local confidence

## Dashboard Design Plan

### Data Structures
Following the AutoDock Vina dashboard pattern:

```python
@dataclass
class BoltzModel:
    """Data structure for individual Boltz-2 model"""
    model_id: int
    confidence_score: float
    ptm: float
    iptm: float
    plddt_mean: float
    pde_mean: float
    structure_file: str
    
@dataclass
class BoltzResults:
    """Complete Boltz-2 results with all models"""
    models: List[BoltzModel]
    best_model: BoltzModel
    protein_name: str
    job_id: str
    prediction_date: str
```

### Core Components

#### 1. BoltzResultsParser
- Parse confidence JSON files
- Load and analyze NumPy arrays (PAE, PDE, pLDDT)
- Extract metadata from file names
- Calculate summary statistics

#### 2. Quality Metrics Analysis
- **Model Ranking**: Sort by confidence_score, ptm, or plddt_mean
- **Confidence Distribution**: Analyze score variations across models
- **Per-Residue Analysis**: Extract confidence profiles from pLDDT arrays
- **Error Analysis**: PAE and PDE matrix analysis

#### 3. Visualization Components

##### Summary Dashboard
- Model comparison table (similar to Vina binding modes)
- Best model identification
- Quality score distributions
- Model ranking by different metrics

##### Detailed Analysis
- **Confidence Heatmaps**: PAE matrices for top models
- **Per-Residue Plots**: pLDDT profiles along sequence
- **Model Comparison**: Side-by-side quality metrics
- **Structure Quality Assessment**: Distance error analysis

##### Interactive Features
- Model selection dropdown
- Metric filtering and sorting
- 3D structure viewer integration (if feasible)
- Export functionality for structures and data

### Technical Implementation

#### Dependencies (same as AutoDock Vina dashboard)
```python
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import json
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Optional
```

#### Key Functions

1. **Data Loading**
   - `load_confidence_data()`: Parse JSON confidence files
   - `load_quality_arrays()`: Load NumPy PAE/PDE/pLDDT data
   - `parse_structure_files()`: Extract metadata from PDB files

2. **Analysis Functions**
   - `calculate_model_rankings()`: Rank models by various metrics
   - `analyze_confidence_distribution()`: Statistical analysis
   - `extract_residue_profiles()`: Per-residue confidence trends
   - `identify_problematic_regions()`: Low confidence regions

3. **Visualization Functions**
   - `create_model_comparison_table()`: Interactive results table
   - `plot_confidence_distributions()`: Score histograms/box plots
   - `create_pae_heatmap()`: PAE matrix visualization
   - `plot_plddt_profile()`: Per-residue confidence line plots
   - `generate_summary_report()`: Overall quality assessment

### Output Features

#### HTML Dashboard Structure
1. **Header**: Job information, protein name, prediction date
2. **Summary Section**: Best model, overall statistics
3. **Model Comparison Table**: All models with sortable columns
4. **Quality Metrics**: Interactive plots and heatmaps
5. **Detailed Analysis**: Per-model breakdown
6. **Export Options**: Download best structures, quality data

#### Quality Assessment Criteria
Based on Boltz-2 documentation and AlphaFold standards:
- **Very High Confidence**: pLDDT > 90, confidence_score > 0.9
- **High Confidence**: pLDDT > 70, confidence_score > 0.7
- **Medium Confidence**: pLDDT > 50, confidence_score > 0.5
- **Low Confidence**: pLDDT < 50, confidence_score < 0.5

### Integration Plan

#### 1. Copy and Adapt AutoDock Vina Dashboard
- Maintain similar structure and style
- Adapt data structures for Boltz outputs
- Keep same HTML template structure
- Use consistent color schemes and layouts

#### 2. Boltz-Specific Enhancements
- PAE heatmap visualization (unique to structure prediction)
- Per-residue confidence profiles
- Model ensemble analysis
- Structure quality assessment tools

#### 3. Usage Pattern
```bash
cd /path/to/boltz/results
python boltz_dashboard.py outputs/
# Generates: boltz_dashboard_YYYYMMDD_HHMMSS.html
```

## Implementation Priority

1. **Phase 1**: Basic parser and model ranking
2. **Phase 2**: Core visualizations (comparison table, confidence plots)
3. **Phase 3**: Advanced analysis (PAE heatmaps, residue profiles)
4. **Phase 4**: Polish and export features

This plan will create a comprehensive analysis tool that leverages the rich quality metrics from Boltz-2 predictions while maintaining consistency with the existing dashboard ecosystem.