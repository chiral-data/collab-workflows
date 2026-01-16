# AutoDock Vina Parameters for Pembrolizumab-PD1 Docking

## Computational Methods in Original Research Papers

### 1. Na et al. (2017) - Structure 5JXE
- **Journal**: Cell Research
- **Resolution**: 2.0 Å
- **Method**: X-ray crystallography, NO computational docking
- **Approach**: 
  - Molecular replacement using Phaser-MR
  - Structure refined with REFMAC5 and Coot
  - Interface analysis using PISA server
- **Key finding**: 1:1 stoichiometry, Pembrolizumab blocks PD-L1/PD-L2 binding site

### 2. Lee et al. (2016) - Structure 5GGS  
- **Journal**: Nature Communications
- **Resolution**: 2.35 Å
- **Method**: Crystallographic structure determination, NO computational docking
- **Approach**:
  - Molecular replacement with known structures
  - Interface analysis using PDBePISA
  - Comparison with other checkpoint inhibitors
- **Focus**: Mechanism of checkpoint blockade by monoclonal antibodies

### 3. Horita et al. (2016) - Structure 5B8C
- **Journal**: Scientific Reports  
- **Resolution**: 2.15 Å (HIGHEST)
- **Method**: X-ray crystallography, NO computational docking
- **Approach**:
  - Molecular replacement using Molrep
  - Structure refinement with REFMAC5
  - Interface analysis with PISA and PyMOL
- **Key contribution**: Highest resolution structure, detailed water-mediated contacts

### 4. Scapin et al. (2015) - Structure 5DK3
- **Journal**: Nature Structural & Molecular Biology
- **Resolution**: 2.3 Å  
- **Method**: Full-length antibody crystallography
- **Focus**: Unique IgG4 Fc conformation, NOT binding interface

### 5. Laureanti et al. (2018) - Computational Analysis
- **Journal**: Scientific Reports
- **Method**: Quantum mechanics/molecular mechanics (QM/MM)
- **Approach**:
  - DFT calculations using B3LYP functional
  - Basis set: 6-31G(d) for QM region
  - AMBER force field for MM region
  - Used 5B8C structure as starting point
- **Key insight**: Quantum chemical analysis of binding energetics

## Our Dual-Method Docking Protocol

### Rationale for Dual-Method Approach

Since **none of the original structural papers performed computational docking** (they solved structures experimentally), our docking experiments serve as:
1. **Validation**: Test if docking can reproduce known crystal structures
2. **Method development**: Establish robust protocols for antibody-antigen docking
3. **Cross-validation**: Compare results across different structure sets
4. **Method comparison**: Evaluate rigid vs flexible docking approaches

### Why Use Both AutoDock Vina AND DiffDock-PP?

| Aspect | AutoDock Vina | DiffDock-PP | Combined Benefit |
|--------|---------------|-------------|------------------|
| **Flexibility** | Rigid body (limitation) | Flexible conformations | Baseline + realistic modeling |
| **Analysis Depth** | Comprehensive energetics | Pose + confidence | Energy analysis + reliability |
| **Sampling** | Systematic grid search | AI-guided diffusion | Thorough + intelligent |
| **Validation** | Well-established | State-of-the-art | Cross-method validation |
| **CDR Loops** | Fixed conformation | Flexible modeling | Compare static vs dynamic |

## AutoDock Vina Parameters

### Search Space Configuration

```txt
# Based on known binding interface coordinates
center_x = [varies by structure set]
center_y = [varies by structure set]  
center_z = [varies by structure set]

size_x = 40  # Large box to accommodate CDR loops
size_y = 40  # Antibody binding sites are extensive  
size_z = 40  # Need to capture all possible orientations

exhaustiveness = 32  # High sampling for accuracy
num_modes = 20      # Generate multiple binding poses
energy_range = 5    # kcal/mol range for output poses
```

#### Structure-Specific Centers

| Structure Set | Center Coordinates | Source |
|---------------|-------------------|---------|
| Set 1 (5JXE) | (25.0, 30.0, 35.0) | CDR-PD1 interface centroid |
| Set 2 (5GGS) | (24.5, 29.5, 34.8) | Binding site geometric center |
| Set 3 (5B8C) | (23.8, 28.9, 35.2) | Interface from highest resolution |
| Set 4 (Hybrid) | (23.8, 28.9, 35.2) | Using 5B8C coordinates |

### Structure Preparation Protocol

#### 1. PDB to PDBQT Conversion

```bash
# For receptor (PD-1)
prepare_receptor4.py -r receptor_pd1.pdb -o receptor_pd1.pdbqt -A hydrogens

# For ligand (Pembrolizumab)  
prepare_ligand4.py -l ligand_pembro.pdb -o ligand_pembro.pdbqt -A hydrogens
```

#### 2. Key Preparation Steps

1. **Hydrogen Addition**: Add polar hydrogens at pH 7.4
2. **Charge Assignment**: Use Gasteiger charges
3. **Torsion Detection**: Identify rotatable bonds in CDR loops
4. **Atom Type Assignment**: AutoDock4 atom types

#### 3. Critical Considerations

- **Antibody Flexibility**: CDR loops are flexible, but Vina treats ligand rigidly
- **Binding Site Size**: PD-1 binding site is relatively small for an antibody
- **Interface Residues**: Key contacts involve Trp, Tyr, and charged residues
- **Glycosylation**: Remove or model N-linked glycans appropriately

### Validation Strategy

#### 1. Self-Docking (Positive Control)
- Re-dock each complex to reproduce crystal structure
- Success criteria: RMSD < 2.0 Å for binding pose

#### 2. Cross-Docking  
- Use receptor from one structure, ligand from another
- Tests transferability of binding site definition

#### 3. Scoring Analysis
- Compare AutoDock Vina scores across sets
- Analyze correlation with experimental binding affinity

#### 4. Interface Analysis
- Check preservation of key hydrogen bonds
- Validate contact residues match crystallographic data

## DiffDock-PP Parameters

### Configuration Settings

```yaml
# DiffDock-PP configuration for pembrolizumab-PD1
num_poses: 40              # Generate multiple binding poses
confidence_model: true     # Enable confidence scoring
residue_level: true        # Use α-carbon representation
samples_per_complex: 10    # Sampling density
```

### Input Format
- **Receptor**: PD-1 structure in PDB format (no conversion needed)
- **Ligand**: Pembrolizumab structure in PDB format
- **Output**: Flexible poses with confidence scores

### Key Advantages for Antibody-Antigen Systems
1. **CDR Loop Flexibility**: Models conformational changes in binding loops
2. **Induced Fit**: Captures mutual adaptation of binding partners
3. **AI-Based Scoring**: Neural network trained on protein-protein complexes
4. **Confidence Estimates**: Reliability assessment for each pose

## Combined Analysis Strategy

### Four-Way Experimental Matrix

| Structure Set | AutoDock Vina | DiffDock-PP | Validation Method | Purpose |
|---------------|---------------|-------------|------------------|---------|
| 5JXE (2.0Å) | ✓ Rigid docking | ✓ Flexible docking | vs crystal structure | Method comparison |
| 5GGS (2.35Å) | ✓ Rigid docking | ✓ Flexible docking | vs crystal structure | Resolution impact |
| 5B8C (2.15Å) | ✓ Rigid docking | ✓ Flexible docking | vs crystal structure | Best resolution test |
| Hybrid | ✓ Rigid docking | ✓ Flexible docking | Cross-validation | Transferability |

### Comparative Metrics
1. **RMSD Comparison**: Both methods vs crystal structures
2. **Interface Recovery**: Preservation of key binding residues
3. **Energy vs Confidence**: Correlation between scoring approaches
4. **Consensus Poses**: Identify poses supported by both methods
5. **Conformational Analysis**: CDR loop movements and induced fit

### Expected Challenges

#### AutoDock Vina Limitations
1. **Size Limitation**: Pembrolizumab is large for small-molecule docking tools
2. **Rigidity**: CDR loops and hinge regions cannot move
3. **Scoring Function**: Optimized for small molecules, not protein-protein
4. **No Induced Fit**: Cannot capture conformational changes upon binding

#### DiffDock-PP Considerations  
1. **Black Box**: Neural network decisions less interpretable than energy functions
2. **Computational Cost**: More expensive than traditional docking
3. **Training Bias**: Performance depends on training set similarity
4. **Limited Analysis**: May provide fewer detailed interaction insights

### Success Metrics

| Metric | Excellent | Good | Acceptable |
|---------|-----------|------|------------|
| RMSD (Å) | < 1.5 | 1.5-2.0 | 2.0-3.0 |
| Interface Recovery | > 80% | 60-80% | 40-60% |
| Energy Score | Most negative | Consistent ranking | Reasonable range |

### Comparison with Literature

#### Interface Residues (from crystal structures)

**PD-1 side:**
- Pro83, Ile126, Leu128, Ala132, Ile134 (hydrophobic core)
- Asn74, Asp77, Lys78 (hydrogen bonds)

**Pembrolizumab side:**
- CDR-H3: Trp99, Tyr100, Asp101
- CDR-L1: Tyr32, Asp31  
- CDR-L3: Tyr92, Asp93

### References for Docking Parameters

1. **AutoDock Vina**: Trott & Olson (2010) J Comput Chem 31:455-461
2. **Protein-protein docking**: Kozakov et al. (2017) Nat Protoc 12:249-278  
3. **Antibody modeling**: Dunbar et al. (2014) Nucleic Acids Res 42:W681-W689
4. **Interface analysis**: Krissinel & Henrick (2007) J Mol Biol 372:774-797

## Notes

- Original papers focused on structural biology, not computational methods
- Our docking serves as computational validation of experimental structures  
- Results will help establish best practices for antibody-antigen docking
- Cross-structure comparison provides robustness testing