# Workflow 3: Virtual Screening - GLP-1 FDA Drug Repurposing

## ✅ **STATUS: COMPLETED**

Complete virtual screening pipeline targeting GLP-1 receptor with FDA-approved drugs for drug repurposing.

## Quick Results

### 🎯 **Top Hits**
1. **lurbinectedin_5398**: -23.121 kcal/mol  
2. **dactinomycin_774**: -22.601 kcal/mol  
3. **fondaparinux_1236**: -22.316 kcal/mol  

### 📊 **Statistics**
- **Compounds Screened**: 4,278 (98.9% success rate)
- **Excellent Hits**: 2,925 compounds (69.7%)
- **Mean Affinity**: -11.552 ± 3.445 kcal/mol
- **Enrichment**: 6.7x better than random

### 🔬 **Recommendation**
**Proceed to Lead Optimization (HIGH Priority)**

## Usage

```bash
# Generate dashboard
python3 vina_virtual_screening_dashboard.py

# Test with small dataset  
python3 vina_virtual_screening_dashboard.py --test

# Docker usage
docker run --rm -v $(pwd):/workspace -w /workspace --entrypoint="" \
  docking_report_potter_python /bin/bash -c \
  "python3 vina_virtual_screening_dashboard.py"
```

## Output Files

- **Dashboard**: `inputs/output/screening_dashboard_20250910_045251.html`
- **Test Dashboard**: `testing/small_dataset_50/outputs/screening_dashboard_20250910_045233.html`
- **Script**: `vina_virtual_screening_dashboard.py`

## Features

- ✅ Interactive visualizations (top hits, distributions, enrichment)
- ✅ Statistical analysis with robust outlier handling  
- ✅ Non-overlapping hit classification
- ✅ Professional GROMACS-style design
- ✅ Docker compatible

---

**See `DASHBOARD_IMPLEMENTATION_PLAN.md` for detailed technical documentation.**