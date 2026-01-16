# DiffDock-PP Antibody-Antigen Docking

## Overview
This workflow performs protein-protein docking using DiffDock-PP for antibody-antigen complexes.

## Important Findings from Testing (September 2025)

### Key Issues Discovered

1. **Structural Distortions in ML-Generated Poses**
   - DiffDock-PP generated poses contain severe geometric distortions
   - Bond lengths range from 1-5 nm (normal protein bonds: 0.1-0.2 nm)
   - These distortions prevent downstream MD simulation with GROMACS

2. **Critical Bug Fixed in Docking Report**
   - Original script used hardcoded file references (e.g., `receptor_pd1.pdb`)
   - Fixed to use dynamic file references based on run_name from job_config.json
   - This bug caused structural corruption when extracting poses

3. **Confidence Score Interpretation**
   - DiffDock-PP confidence scores: HIGHER is BETTER (opposite of RMSD)
   - Typical range: -20 to +5
   - Best poses typically have positive or near-zero scores

### Compatibility with MD Simulation

**Current Status**: DiffDock-PP poses are NOT directly compatible with GROMACS MD simulation due to:
- Missing atoms in aromatic rings (e.g., incomplete HIS107)
- Extreme bond distortions incompatible with classical force fields
- Both CHARMM27 and AMBER99sb-ildn force fields fail to process the structures

**Potential Solutions**:
1. Use structure repair tools (PDBFixer, OpenMM) before MD
2. Perform gentle energy minimization with strong position restraints
3. Use coarse-grained simulation (MARTINI force field)
4. Extract only the binding interface region for focused MD
5. Use original crystal structures for MD validation

## Workflow Steps

1. **Input Preparation**: Separate receptor and ligand PDB files
2. **DiffDock-PP Execution**: Generate multiple binding poses
3. **Pose Analysis**: Extract and rank poses by confidence score
4. **Structure Validation**: Check for geometric distortions before MD

## File Structure
- `job_config.json`: Configuration for each docking run
- `job_script.sh`: Execution script for Docker container
- `inputs/`: Input PDB files
- `outputs/`: Generated poses and predictions

## Notes
- Always validate ML-generated structures before physics-based simulation
- Consider using the docking poses for interaction analysis rather than MD
- For production MD, consider starting from experimental structures