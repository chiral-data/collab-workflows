#!/bin/bash
# AutoDock Vina Test Script - 30 Compounds
# Based on production job_script.sh - Docker Version

set -e  # Exit on error

# Docker image configuration
DOCKER_IMAGE="chiral.sakuracr.jp/autodock_vina_potter_python:latest"
CURRENT_DIR=$(pwd)

# Configuration
JOB_ID="glp1_test_50_compounds"
RECEPTOR_FILE="7s15_glp1r.pdb"
LIGANDS_DIR="inputs/pdbqt_files"

# Box parameters (same as production)
CENTER_X=69.69
CENTER_Y=69.69
CENTER_Z=60.66
SIZE_X=25.0
SIZE_Y=25.0
SIZE_Z=25.0

# Docking parameters (same as production)
EXHAUSTIVENESS=4
NUM_MODES=9
ENERGY_RANGE=3

# Set output directory
OUTPUT_DIR="outputs"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/poses"
mkdir -p "$OUTPUT_DIR/logs"

echo "==========================================="
echo "Starting Test Virtual Screening Job: $JOB_ID"
echo "==========================================="
echo "Receptor: $RECEPTOR_FILE"
echo "Ligands directory: $LIGANDS_DIR"
echo "Box center: ($CENTER_X, $CENTER_Y, $CENTER_Z)"
echo "Box size: ($SIZE_X, $SIZE_Y, $SIZE_Z)"
echo "Exhaustiveness: $EXHAUSTIVENESS"
echo "Number of modes: $NUM_MODES"
echo "Energy range: $ENERGY_RANGE"
echo ""

# Step 1: Prepare receptor (copy PDBQT directly)
echo "Step 1: Preparing receptor..."
RECEPTOR_BASE=$(basename "$RECEPTOR_FILE" | sed 's/\.[^.]*$//')
# Convert receptor PDB to PDBQT using OpenBabel
obabel -i pdb "$RECEPTOR_FILE" -o pdbqt -O "${RECEPTOR_BASE}.pdbqt" -xr

# Skip Docker test when running inside container
echo "Running inside Docker container - skipping nested Docker test"

# Step 2: Create configuration file for Vina
echo ""
echo "Step 2: Creating Vina configuration file..."
cat > vina_config.txt << EOF
center_x = $CENTER_X
center_y = $CENTER_Y
center_z = $CENTER_Z
size_x = $SIZE_X
size_y = $SIZE_Y
size_z = $SIZE_Z
exhaustiveness = $EXHAUSTIVENESS
num_modes = $NUM_MODES
energy_range = $ENERGY_RANGE
EOF

echo "Configuration file created:"
cat vina_config.txt

# Step 3: Count ligands and prepare for batch processing
echo ""
echo "Step 3: Counting ligands for batch processing..."
LIGAND_COUNT=$(find "$LIGANDS_DIR" -name "*.pdbqt" | wc -l)
echo "Found $LIGAND_COUNT PDBQT ligand files to process"

if [ "$LIGAND_COUNT" -eq 0 ]; then
    echo "Error: No PDBQT files found in $LIGANDS_DIR"
    exit 1
fi

# Step 4: Convert SDF files to proper PDBQT format (using working method)
echo ""
echo "Step 4: Converting SDF files to proper PDBQT format..."
mkdir -p converted_ligands_pdbqt
CONVERTED_COUNT=0
FAILED_COUNT=0

# Use SDF files in current directory for testing
SDF_DIR="."
SDF_COUNT=0
for sdf_file in "$SDF_DIR"/*.sdf; do
    if [ "$SDF_COUNT" -ge 50 ]; then
        break
    fi
    
    if [ -f "$sdf_file" ]; then
        ligand_base=$(basename "$sdf_file" .sdf)
        pdbqt_file="converted_ligands_pdbqt/${ligand_base}.pdbqt"
        
        echo "Converting: $ligand_base"
        
        # Use proper SDF to PDBQT conversion with timeout
        temp_file="temp_${ligand_base}.pdbqt"
        echo -n "  Converting with timeout... "
        if timeout 30s obabel -i sdf "$sdf_file" -o pdbqt -O "$temp_file" -xr 2>/dev/null && [ -f "$temp_file" ]; then
            echo "Success"
            # Add ROOT/ENDROOT/TORSDOF structure as per working AutoDock script
            sed -i '/^TER/d' "$temp_file" 2>/dev/null || true
            sed -i '1i\ROOT' "$temp_file"
            echo 'ENDROOT' >> "$temp_file"
            echo 'TORSDOF 0' >> "$temp_file"
            mv "$temp_file" "$pdbqt_file"
            CONVERTED_COUNT=$((CONVERTED_COUNT + 1))
        else
            echo "Failed (timeout or error)"
            echo "WARNING: Failed to convert $sdf_file" >> "$OUTPUT_DIR/logs/conversion_failures.log"
            FAILED_COUNT=$((FAILED_COUNT + 1))
            # Clean up any leftover temp file
            rm -f "$temp_file" 2>/dev/null || true
        fi
        
        SDF_COUNT=$((SDF_COUNT + 1))
    fi
done

echo "Conversion complete: $CONVERTED_COUNT successful, $FAILED_COUNT failed"

# Step 5: Create ligand list file 
echo ""
echo "Step 5: Creating ligand batch list..."
LIGAND_LIST="ligand_batch_list.txt"
find "converted_ligands_pdbqt" -name "*.pdbqt" > "$LIGAND_LIST"
LIGAND_COUNT=$(cat "$LIGAND_LIST" | wc -l)

echo "Ligand list created with $LIGAND_COUNT entries"

# Step 6: Run batch docking with AutoDock Vina (exact same as production)
echo ""
echo "Step 6: Running batch docking with AutoDock Vina..."

DOCKING_COUNT=0
TOTAL_LIGANDS=$LIGAND_COUNT

while IFS= read -r ligand_file; do
    if [ -f "$ligand_file" ]; then
        ligand_base=$(basename "$ligand_file" .pdbqt)
        output_file="$OUTPUT_DIR/poses/${ligand_base}_out.pdbqt"
        log_file="$OUTPUT_DIR/logs/${ligand_base}.log"
        
        # Run individual docking (direct command, no nested Docker)
        if vina --receptor "${RECEPTOR_BASE}.pdbqt" \
               --ligand "$ligand_file" \
               --config vina_config.txt \
               --out "$output_file" > "$log_file" 2>&1; then
            DOCKING_COUNT=$((DOCKING_COUNT + 1))
        else
            echo "FAILED: $ligand_base" >> "$OUTPUT_DIR/logs/docking_failures.log"
        fi
        
        # Progress indicator (every 10 for small test set)
        if [ $((DOCKING_COUNT % 10)) -eq 0 ] && [ $DOCKING_COUNT -gt 0 ]; then
            echo "Completed docking for $DOCKING_COUNT/$TOTAL_LIGANDS ligands..."
        fi
    fi
done < "$LIGAND_LIST"

echo "Batch docking complete: $DOCKING_COUNT successful dockings"

# Step 7: Extract and summarize results (exact same as production)
echo ""
echo "Step 7: Extracting binding scores and creating summary..."

RESULTS_FILE="$OUTPUT_DIR/screening_results.csv"
echo "Ligand_Name,Best_Score,Second_Score,Third_Score" > "$RESULTS_FILE"

for pose_file in "$OUTPUT_DIR/poses"/*_out.pdbqt; do
    if [ -f "$pose_file" ]; then
        ligand_name=$(basename "$pose_file" _out.pdbqt)
        
        # Extract top 3 scores from PDBQT file
        scores=$(grep "^REMARK VINA RESULT:" "$pose_file" | head -3 | awk '{print $4}' | tr '\n' ',' | sed 's/,$//')
        
        if [ ! -z "$scores" ]; then
            echo "$ligand_name,$scores" >> "$RESULTS_FILE"
        fi
    fi
done

# Sort results by best score (most negative = best)
echo ""
echo "Creating sorted results file..."
head -1 "$RESULTS_FILE" > "$OUTPUT_DIR/screening_results_sorted.csv"
tail -n +2 "$RESULTS_FILE" | sort -t, -k2,2n >> "$OUTPUT_DIR/screening_results_sorted.csv"

# Step 8: Generate summary report
echo ""
echo "Step 8: Generating summary report..."
cat > "$OUTPUT_DIR/test_screening_summary.txt" << EOF
AutoDock Vina Test Virtual Screening Results
===========================================
Job ID: $JOB_ID
Date: $(date)
Target: GLP-1 Receptor FDA Drug Test (30 compounds)

Input Information:
- Receptor: $RECEPTOR_FILE
- Ligand Directory: $LIGANDS_DIR
- Total PDBQT files found: $LIGAND_COUNT

Docking Parameters:
- Box Center: ($CENTER_X, $CENTER_Y, $CENTER_Z)
- Box Size: ($SIZE_X, $SIZE_Y, $SIZE_Z)
- Exhaustiveness: $EXHAUSTIVENESS
- Number of modes: $NUM_MODES
- Energy range: $ENERGY_RANGE

Results Summary:
- Successful dockings: $DOCKING_COUNT
- Failed dockings: $((LIGAND_COUNT - DOCKING_COUNT))
- Success rate: $(( DOCKING_COUNT * 100 / LIGAND_COUNT ))%

Output Files:
- Individual poses: poses/ directory (${DOCKING_COUNT} files)
- Docking logs: logs/ directory
- Results CSV: screening_results.csv
- Sorted results: screening_results_sorted.csv

Top 10 Best Scoring Compounds:
EOF

# Add top 10 results to summary
echo "" >> "$OUTPUT_DIR/test_screening_summary.txt"
head -11 "$OUTPUT_DIR/screening_results_sorted.csv" | tail -10 | while IFS=, read -r name score rest; do
    echo "- $name: $score kcal/mol" >> "$OUTPUT_DIR/test_screening_summary.txt"
done

# List all output files
echo ""
echo "Output files generated:"
ls -la "$OUTPUT_DIR/"

echo ""
echo "==========================================="
echo "Test Virtual Screening job completed!"
echo "==========================================="
echo "Results summary:"
echo "- Total ligands processed: $LIGAND_COUNT"
echo "- Successful dockings: $DOCKING_COUNT"
echo "- Success rate: $(( DOCKING_COUNT * 100 / LIGAND_COUNT ))%"
echo "- Results file: $OUTPUT_DIR/screening_results_sorted.csv"
echo "- Summary report: $OUTPUT_DIR/test_screening_summary.txt"
echo "==========================================="