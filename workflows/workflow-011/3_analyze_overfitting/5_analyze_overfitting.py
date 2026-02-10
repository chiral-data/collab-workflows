"""
HYBRID MODEL ANALYSIS
Analyzes the trained model from Node 2.
Generates a single 9-panel diagnostic plot with training curves.
"""
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sklearn.metrics
from keras.models import load_model
import tensorflow as tf
import sys
import os
import pickle
import json

# Ensure output directory exists
OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Create log file
log_path = os.path.join(OUTPUT_DIR, 'model_analysis_log.txt')
log_file = open(log_path, 'w')
def log_print(message):
    print(message)
    log_file.write(message + '\n')
    log_file.flush()

# Custom metrics for model loading
def rmse(y_true, y_pred):
    from tensorflow.keras import backend as K
    return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1))

def r_square(y_true, y_pred):
    from tensorflow.keras import backend as K
    SS_res = K.sum(K.square(y_true - y_pred))
    SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
    return (1 - SS_res / (SS_tot + K.epsilon()))

# ---------------------------------------------------------
# MAIN ANALYSIS
# ---------------------------------------------------------
print("="*80)
print("HYBRID MODEL ANALYSIS (JSON EXPORT)")
print("="*80)

# 1. Load Data from Node 1
print("\n1. Loading preprocessed data from Node 1...")
node1_output = "../1_data_prep/outputs/processed_data.pkl"

if not os.path.exists(node1_output):
    print(f"❌ Error: {node1_output} not found! Did Node 1 run successfully?")
    sys.exit(1)

with open(node1_output, 'rb') as f:
    data_bundle = pickle.load(f)

X_train = data_bundle['X_train']
X_test = data_bundle['X_test']
y_train = data_bundle['y_train']
y_test = data_bundle['y_test']

print(f"   ✓ Train set: {X_train.shape}")
print(f"   ✓ Test set: {X_test.shape}")

# 2. Load Model and Training History from Node 2
print("\n2. Loading model and training history from Node 2...")
node2_output = "../2_model_train/outputs"

# Load model
model_path = os.path.join(node2_output, 'hybrid_model.keras')
if not os.path.exists(model_path):
    print(f"❌ Error: {model_path} not found! Did Node 2 run successfully?")
    sys.exit(1)

model = load_model(model_path, custom_objects={'rmse': rmse, 'r_square': r_square})
print(f"   ✓ Loaded model")

# Load training history
history_path = os.path.join(node2_output, 'training_history.pkl')
if not os.path.exists(history_path):
    print(f"❌ Error: {history_path} not found! Did Node 2 save training history?")
    sys.exit(1)

with open(history_path, 'rb') as f:
    history = pickle.load(f)
print(f"   ✓ Loaded training history")

# 3. Evaluate Model
print("\n3. Evaluating model...")

# Predictions
train_pred = model.predict(X_train, verbose=0).flatten()
test_pred = model.predict(X_test, verbose=0).flatten()

# Metrics
train_r2 = sklearn.metrics.r2_score(y_train, train_pred)
test_r2 = sklearn.metrics.r2_score(y_test, test_pred)
train_rmse = math.sqrt(sklearn.metrics.mean_squared_error(y_train, train_pred))
test_rmse = math.sqrt(sklearn.metrics.mean_squared_error(y_test, test_pred))
train_mae = sklearn.metrics.mean_absolute_error(y_train, train_pred)
test_mae = sklearn.metrics.mean_absolute_error(y_test, test_pred)
r2_gap = abs(train_r2 - test_r2) * 100

print(f"   Train R²={train_r2:.4f}, Test R²={test_r2:.4f}, Gap={r2_gap:.2f}%")

# Overfitting Verdict
if r2_gap > 20: verdict = "SEVERE OVERFIT"
elif r2_gap > 10: verdict = "MODERATE OVERFIT"
elif r2_gap > 5: verdict = "MINIMAL OVERFIT"
else: verdict = "NO OVERFITTING"

# Calculate params and ratio
total_params = model.count_params()
data_to_param_ratio = len(X_train) / total_params

# 4. Prepare Data for Export
print("\n4. Exporting data to JSON...")

export_data = {
    "metrics": {
        "train_r2": float(train_r2),
        "test_r2": float(test_r2),
        "train_rmse": float(train_rmse),
        "test_rmse": float(test_rmse),
        "train_mae": float(train_mae),
        "test_mae": float(test_mae),
        "r2_gap": float(r2_gap),
        "total_params": int(total_params),
        "data_ratio": float(data_to_param_ratio),
        "verdict": verdict
    },
    "residuals": {
        "train": {
            "actual": y_train.tolist(),
            "predicted": train_pred.tolist(),
            "residual": (y_train - train_pred).tolist()
        },
        "test": {
            "actual": y_test.tolist(),
            "predicted": test_pred.tolist(),
            "residual": (y_test - test_pred).tolist()
        }
    },
    "history": history  # Already a dict of lists
}

# Fix numpy types in history if any
def convert_numpy(obj):
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    elif isinstance(obj, (np.floating, float)):
        return float(obj)
    elif isinstance(obj, (np.ndarray, list)):
        return [convert_numpy(x) for x in obj]
    return obj

export_data = json.loads(json.dumps(export_data, default=convert_numpy))

json_path = os.path.join(OUTPUT_DIR, 'data.json')
with open(json_path, 'w') as f:
    json.dump(export_data, f)

print(f"✓ Saved: {json_path}")
print("\n" + "="*80)
print("NODE 3 ANALYSIS COMPLETE")
print("="*80)
