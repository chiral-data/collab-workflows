#!/bin/bash
# AutoDock Vina Virtual Screening Job Script for Potter Platform
# GLP-1 FDA Drug Repurposing - Batch Docking

set -e  # Exit on error

# Configuration file
CONFIG_FILE="job_config.json"

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file '$CONFIG_FILE' not found"
    exit 1
fi

# Parse JSON config
echo "Reading configuration from $CONFIG_FILE..."
JOB_ID=$(jq -r '.job_id' "$CONFIG_FILE")
RECEPTOR_FILE=$(jq -r '.inputs[0].filename' "$CONFIG_FILE")
LIGANDS_DIR=$(jq -r '.inputs[1].directory' "$CONFIG_FILE")

# Box parameters (from config)
CENTER_X=$(jq -r '.parameters.center_x' "$CONFIG_FILE")
CENTER_Y=$(jq -r '.parameters.center_y' "$CONFIG_FILE")
CENTER_Z=$(jq -r '.parameters.center_z' "$CONFIG_FILE")
SIZE_X=$(jq -r '.parameters.size_x' "$CONFIG_FILE")
SIZE_Y=$(jq -r '.parameters.size_y' "$CONFIG_FILE")
SIZE_Z=$(jq -r '.parameters.size_z' "$CONFIG_FILE")

# Docking parameters
EXHAUSTIVENESS=$(jq -r '.parameters.exhaustiveness // 8' "$CONFIG_FILE")
NUM_MODES=$(jq -r '.parameters.num_modes // 9' "$CONFIG_FILE")
ENERGY_RANGE=$(jq -r '.parameters.energy_range // 3' "$CONFIG_FILE")

# Set output directory
OUTPUT_DIR="/workspace/output"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/poses"
mkdir -p "$OUTPUT_DIR/logs"

echo "==========================================="
echo "Starting Virtual Screening Job: $JOB_ID"
echo "==========================================="
echo "Receptor: $RECEPTOR_FILE"
echo "Ligands directory: $LIGANDS_DIR"
echo "Box center: ($CENTER_X, $CENTER_Y, $CENTER_Z)"
echo "Box size: ($SIZE_X, $SIZE_Y, $SIZE_Z)"
echo "Exhaustiveness: $EXHAUSTIVENESS"
echo "Number of modes: $NUM_MODES"
echo "Energy range: $ENERGY_RANGE"
echo ""

# Step 1: Prepare receptor (convert to PDBQT using OpenBabel)
echo "Step 1: Preparing receptor using OpenBabel..."
RECEPTOR_BASE=$(basename "$RECEPTOR_FILE" | sed 's/\.[^.]*$//')
case "$RECEPTOR_FILE" in
    *.pdb)
        echo "Converting receptor PDB to PDBQT..."
        obabel -i pdb "$RECEPTOR_FILE" -o pdbqt -O "${RECEPTOR_BASE}.pdbqt" -xr
        ;;
    *.pdbqt)
        echo "Using PDBQT receptor directly..."
        cp "$RECEPTOR_FILE" "${RECEPTOR_BASE}.pdbqt"
        ;;
    *)
        echo "Error: Unsupported receptor file format. Use PDB or PDBQT."
        exit 1
        ;;
esac

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
LIGAND_COUNT=$(find "$LIGANDS_DIR" -name "*.sdf" | wc -l)
echo "Found $LIGAND_COUNT SDF ligand files to process"

if [ "$LIGAND_COUNT" -eq 0 ]; then
    echo "Error: No SDF files found in $LIGANDS_DIR"
    exit 1
fi

# Step 4: Convert all ligands to PDBQT and prepare batch file list
echo ""
echo "Step 4: Converting ligands to PDBQT format..."
mkdir -p ligands_pdbqt
LIGAND_LIST="ligand_batch_list.txt"
> "$LIGAND_LIST"  # Clear the file

CONVERTED_COUNT=0
FAILED_COUNT=0

for sdf_file in "$LIGANDS_DIR"/*.sdf; do
    if [ -f "$sdf_file" ]; then
        ligand_base=$(basename "$sdf_file" .sdf)
        pdbqt_file="ligands_pdbqt/${ligand_base}.pdbqt"
        
        echo "Converting: $ligand_base"
        
        # Convert SDF to PDBQT with timeout and proper formatting
        temp_file="temp_${ligand_base}.pdbqt"
        if timeout 30s obabel -i sdf "$sdf_file" -o pdbqt -O "$temp_file" -xr 2>/dev/null && [ -f "$temp_file" ]; then
            # Add ROOT/ENDROOT/TORSDOF structure as per working test script
            sed -i '/^TER/d' "$temp_file" 2>/dev/null || true
            sed -i '1i\ROOT' "$temp_file"
            echo 'ENDROOT' >> "$temp_file"
            echo 'TORSDOF 0' >> "$temp_file"
            mv "$temp_file" "$pdbqt_file"
            echo "$pdbqt_file" >> "$LIGAND_LIST"
            CONVERTED_COUNT=$((CONVERTED_COUNT + 1))
        else
            echo "WARNING: Failed to convert $sdf_file (timeout or error)" >> "$OUTPUT_DIR/logs/conversion_failures.log"
            FAILED_COUNT=$((FAILED_COUNT + 1))
            # Clean up any leftover temp file
            rm -f "$temp_file" 2>/dev/null || true
        fi
        
        # Progress indicator
        if [ $((CONVERTED_COUNT % 100)) -eq 0 ]; then
            echo "Processed $CONVERTED_COUNT ligands so far..."
        fi
    fi
done

echo "Conversion complete: $CONVERTED_COUNT successful, $FAILED_COUNT failed"

# Step 5: Run batch docking with AutoDock Vina
echo ""
echo "Step 5: Running batch docking with AutoDock Vina..."

DOCKING_COUNT=0
TOTAL_LIGANDS=$CONVERTED_COUNT

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
        
        # Progress indicator
        if [ $((DOCKING_COUNT % 50)) -eq 0 ]; then
            echo "Completed docking for $DOCKING_COUNT/$TOTAL_LIGANDS ligands..."
        fi
    fi
done < "$LIGAND_LIST"

echo "Batch docking complete: $DOCKING_COUNT successful dockings"

# Step 6: Extract and summarize results
echo ""
echo "Step 6: Extracting binding scores and creating summary..."

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

# Step 7: Generate summary report
echo ""
echo "Step 7: Generating summary report..."
cat > "$OUTPUT_DIR/virtual_screening_summary.txt" << EOF
AutoDock Vina Virtual Screening Results
======================================
Job ID: $JOB_ID
Date: $(date)
Target: GLP-1 Receptor FDA Drug Repurposing

Input Information:
- Receptor: $RECEPTOR_FILE
- Ligand Library: $LIGANDS_DIR
- Total SDF files found: $LIGAND_COUNT
- Successfully converted to PDBQT: $CONVERTED_COUNT
- Conversion failures: $FAILED_COUNT

Docking Parameters:
- Box Center: ($CENTER_X, $CENTER_Y, $CENTER_Z)
- Box Size: ($SIZE_X, $SIZE_Y, $SIZE_Z)
- Exhaustiveness: $EXHAUSTIVENESS
- Number of modes: $NUM_MODES
- Energy range: $ENERGY_RANGE

Results Summary:
- Successful dockings: $DOCKING_COUNT
- Failed dockings: $((CONVERTED_COUNT - DOCKING_COUNT))

Output Files:
- Individual poses: poses/ directory (${DOCKING_COUNT} files)
- Docking logs: logs/ directory
- Results CSV: screening_results.csv
- Sorted results: screening_results_sorted.csv
- Failure logs: logs/conversion_failures.log, logs/docking_failures.log

Top 10 Best Scoring Compounds:
EOF

# Add top 10 results to summary
echo "" >> "$OUTPUT_DIR/virtual_screening_summary.txt"
head -11 "$OUTPUT_DIR/screening_results_sorted.csv" | tail -10 | while IFS=, read -r name score rest; do
    echo "- $name: $score kcal/mol" >> "$OUTPUT_DIR/virtual_screening_summary.txt"
done

# List all output files
echo ""
echo "Output files generated:"
ls -la "$OUTPUT_DIR/"

echo ""
echo "==========================================="
echo "Virtual Screening job completed successfully!"
echo "==========================================="
echo "Results summary:"
echo "- Total ligands processed: $CONVERTED_COUNT"
echo "- Successful dockings: $DOCKING_COUNT"
echo "- Results file: $OUTPUT_DIR/screening_results_sorted.csv"
echo "- Summary report: $OUTPUT_DIR/virtual_screening_summary.txt"
echo "==========================================="