#!/bin/bash

CONFIG_FILE="job_config.json"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file '$CONFIG_FILE' not found"
    exit 1
fi

JOB_ID=$(jq -r '.job_id' "$CONFIG_FILE")
RECEPTOR_FILE=$(jq -r '.inputs[0].filename' "$CONFIG_FILE")
LIGAND_FILE=$(jq -r '.inputs[1].filename' "$CONFIG_FILE")
NUM_SAMPLES=$(jq -r '.parameters.num_samples // 40' "$CONFIG_FILE")
OUTPUT_DIR=$(jq -r '.output.directory // "/workspace/output"' "$CONFIG_FILE")
RUN_NAME=$(jq -r '.parameters.run_name // "4G6K"' "$CONFIG_FILE")

echo "Starting DiffDock-PP job: $JOB_ID"
echo "Run name: $RUN_NAME, Samples: $NUM_SAMPLES"

mkdir -p "$OUTPUT_DIR"

# Files are expected in /workspace directory
if [ ! -f "/workspace/$RECEPTOR_FILE" ]; then
    echo "Error: Receptor file '/workspace/$RECEPTOR_FILE' not found"
    exit 1
fi

if [ ! -f "/workspace/$LIGAND_FILE" ]; then
    echo "Error: Ligand file '/workspace/$LIGAND_FILE' not found"
    exit 1
fi

# Work in /workspace since /opt/DiffDock-PP is read-only
cd /workspace

# Create required working directories
mkdir -p "/workspace/storage" "/workspace/visualization" "/workspace/ckpts"

COMPLEX_NAME="${RUN_NAME}_complex"

# Create dataset directory structure (DiffDock-PP expects this exact path)
mkdir -p "/workspace/datasets/single_pair_dataset/structures"
cp "/workspace/$RECEPTOR_FILE" "/workspace/datasets/single_pair_dataset/structures/${COMPLEX_NAME}_r_b.pdb"
cp "/workspace/$LIGAND_FILE" "/workspace/datasets/single_pair_dataset/structures/${COMPLEX_NAME}_l_b.pdb"

# Create splits file (DiffDock-PP hardcoded to look here)
echo 'path,split' > "/workspace/datasets/single_pair_dataset/splits_test.csv"
echo "${COMPLEX_NAME},test" >> "/workspace/datasets/single_pair_dataset/splits_test.csv"

# Run DiffDock-PP inference
PYTHONPATH=/opt/DiffDock-PP/src python /opt/DiffDock-PP/src/main_inf.py \
    --mode "test" \
    --config_file /opt/DiffDock-PP/config/single_pair_inference.yaml \
    --data_file "/workspace/datasets/single_pair_dataset/splits_test.csv" \
    --data_path "/workspace/datasets/single_pair_dataset" \
    --run_name "$RUN_NAME" \
    --save_path "/workspace/ckpts/$RUN_NAME" \
    --batch_size 1 \
    --num_folds 1 \
    --num_gpu 1 \
    --gpu 0 \
    --seed 0 \
    --visualization_path "/workspace/visualization/$RUN_NAME" \
    --visualize_n_val_graphs $NUM_SAMPLES \
    --filtering_model_path /opt/DiffDock-PP/checkpoints/confidence_model_dips/fold_0/ \
    --score_model_path /opt/DiffDock-PP/checkpoints/large_model_dips/fold_0/ \
    --num_samples $NUM_SAMPLES \
    --prediction_storage "/workspace/storage/${RUN_NAME}_predictions.pkl"

# Copy results to output
if [ -f "/workspace/storage/${RUN_NAME}_predictions.pkl" ]; then
    cp "/workspace/storage/${RUN_NAME}_predictions.pkl" "$OUTPUT_DIR/"
fi

# Copy pose files to output
if [ -d "/workspace/visualization/epoch-0/${RUN_NAME}_complex" ]; then
    mkdir -p "$OUTPUT_DIR/poses"
    find "/workspace/visualization/epoch-0/${RUN_NAME}_complex" -name "*.pdb" -exec cp {} "$OUTPUT_DIR/poses/" \;
fi

# Clean up temporary directories
echo "Cleaning up temporary directories..."
rm -rf "/workspace/datasets" "/workspace/storage" "/workspace/visualization" "/workspace/ckpts" "/workspace/torchhub"

if [ $? -eq 0 ]; then
    echo "DiffDock-PP inference completed successfully!"
    
    echo "DiffDock-PP Results Summary" > "$OUTPUT_DIR/summary.txt"
    echo "==========================" >> "$OUTPUT_DIR/summary.txt"
    echo "Receptor: $RECEPTOR_FILE" >> "$OUTPUT_DIR/summary.txt"
    echo "Ligand: $LIGAND_FILE" >> "$OUTPUT_DIR/summary.txt"
    echo "Number of poses: $NUM_SAMPLES" >> "$OUTPUT_DIR/summary.txt"
    echo "Complex name: $COMPLEX_NAME" >> "$OUTPUT_DIR/summary.txt"
    echo "" >> "$OUTPUT_DIR/summary.txt"
    echo "Output files:" >> "$OUTPUT_DIR/summary.txt"
    ls "$OUTPUT_DIR/" >> "$OUTPUT_DIR/summary.txt"
else
    echo "DiffDock-PP inference failed!"
    exit 1
fi

echo "Job completed."