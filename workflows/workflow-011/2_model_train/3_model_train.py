"""
QSAR Workflow Node 2: Model Training
Responsibilities:
- Load Processed Data (from Node 1)
- Define Hybrid Model Architecture (Reg + BatchNorm)
- Apply Label Smoothing & Noise Injection
- Train Ensemble of 5 Models
- Evaluate Performance (R2, RMSE)
- Save Final Models
"""
import os
import sys
import json
import math
import pickle
import numpy as np
import matplotlib.pyplot as plt
import sklearn.metrics
import tensorflow as tf
from keras.models import Sequential
from keras.layers import Dense, Dropout, Input, BatchNormalization
from keras.regularizers import l2
from keras.callbacks import EarlyStopping, ReduceLROnPlateau
from numpy.random import seed

# Reproducibility
seed(42)
tf.random.set_seed(42)
np.random.seed(42)

# Ensure output directory exists (Local to this node)
OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# -------------------------
# Custom Metrics
# -------------------------
def rmse(y_true, y_pred):
    from tensorflow.keras import backend as K
    return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1))

def r_square(y_true, y_pred):
    from tensorflow.keras import backend as K
    SS_res = K.sum(K.square(y_true - y_pred))
    SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
    return (1 - SS_res / (SS_tot + K.epsilon()))

# -------------------------
# Utilities
# -------------------------
class FullMetricsCallback(tf.keras.callbacks.Callback):
    """Calculate sklearn metrics on clean data at each epoch"""
    def __init__(self, X_train, y_train, X_test, y_test):
        self.X_train = X_train
        self.y_train = y_train
        self.X_test = X_test
        self.y_test = y_test

    def on_epoch_end(self, epoch, logs=None):
        if logs is None:
            logs = {}
        
        # Predict on clean data
        x_pred = self.model.predict(self.X_train, verbose=0)
        y_pred = self.model.predict(self.X_test, verbose=0)
        
        x_pred = np.nan_to_num(x_pred, nan=0.0)
        y_pred = np.nan_to_num(y_pred, nan=0.0)
        
        # R²
        train_r2 = sklearn.metrics.r2_score(self.y_train, x_pred)
        val_r2 = sklearn.metrics.r2_score(self.y_test, y_pred)
        logs['train_r2_sklearn'] = train_r2
        logs['val_r2_sklearn'] = val_r2
        
        # MAE
        train_mae = sklearn.metrics.mean_absolute_error(self.y_train, x_pred)
        val_mae = sklearn.metrics.mean_absolute_error(self.y_test, y_pred)
        logs['train_mae_sklearn'] = train_mae
        logs['val_mae_sklearn'] = val_mae

        # RMSE
        train_rmse = math.sqrt(sklearn.metrics.mean_squared_error(self.y_train, x_pred))
        val_rmse = math.sqrt(sklearn.metrics.mean_squared_error(self.y_test, y_pred))
        logs['train_rmse_sklearn'] = train_rmse
        logs['val_rmse_sklearn'] = val_rmse


def add_noise_to_features(X, noise_level=0.01):
    return X + np.random.normal(0, noise_level, X.shape)

def smooth_labels(y, epsilon=0.1):
    return y * (1 - epsilon) + epsilon * np.mean(y)

def build_hybrid_model(input_dim):
    """Build Model with Batch Normalization and Strong Regularization"""
    model = Sequential()
    model.add(Input(shape=(input_dim,)))
    
    # Layer 1
    model.add(Dense(192, activation='relu', kernel_regularizer=l2(0.002)))
    model.add(BatchNormalization())
    model.add(Dropout(0.4))
    
    # Layer 2
    model.add(Dense(96, activation='relu', kernel_regularizer=l2(0.002)))
    model.add(BatchNormalization())
    model.add(Dropout(0.3))
    
    # Layer 3
    model.add(Dense(48, activation='relu', kernel_regularizer=l2(0.001)))
    model.add(BatchNormalization())
    model.add(Dropout(0.2))
    
    # Output
    model.add(Dense(1, activation='linear'))
    
    optimizer = tf.keras.optimizers.Adam(learning_rate=0.001, clipnorm=1.0)
    # Add metrics back so they're tracked
    model.compile(loss='mean_squared_error', optimizer=optimizer, metrics=['mae', rmse, r_square])
    return model

# -------------------------
# Main Execution
# -------------------------
if __name__ == "__main__":
    print("="*80)
    print("NODE 2: MODEL TRAINING (ULTIMATE HYBRID)")
    print("="*80)

    # 1. Load Data from Node 1 output
    loaded_path = "inputs/processed_data.pkl"

    if not os.path.exists(loaded_path):
        print(f"❌ Error: {loaded_path} not found! Did Node 1 run successfully?")
        sys.exit(1)

    print(f"Step 1: Loading data from {loaded_path}...")
    with open(loaded_path, 'rb') as f:
        data = pickle.load(f)
        
    X_train = data['X_train']
    y_train = data['y_train']
    X_test = data['X_test']
    y_test = data['y_test']
    # Safely get SMILES if they exist
    smiles_train = data.get('smiles_train', np.array([""] * len(X_train)))
    smiles_test = data.get('smiles_test', np.array([""] * len(X_test)))
    
    print(f"✓ Training set: {X_train.shape}")
    print(f"✓ Test set: {X_test.shape}")

    # ----------------------------------------------------
    # DASHBOARD DATA STRUCTURE
    # ----------------------------------------------------
    dashboard_data = {
        "setup": {
            "s1_architecture": {
                "layers": [
                    {"type": "Input", "units": X_train.shape[1]},
                    {"type": "Dense", "units": 192, "activation": "relu", "l2": 0.002},
                    {"type": "BatchNorm"},
                    {"type": "Dropout", "rate": 0.4},
                    {"type": "Dense", "units": 96, "activation": "relu", "l2": 0.002},
                    {"type": "BatchNorm"},
                    {"type": "Dropout", "rate": 0.3},
                    {"type": "Dense", "units": 48, "activation": "relu", "l2": 0.001},
                    {"type": "BatchNorm"},
                    {"type": "Dropout", "rate": 0.2},
                    {"type": "Output", "units": 1, "activation": "linear"}
                ]
            },
            "s1_hyperparams": {
                "learning_rate": 0.001,
                "batch_size": 256,
                "epochs": 200,
                "optimizer": "Adam"
            }
        },
        "augmentation": {},
        "training": {
            "history": {"loss": [], "val_loss": [], "train_r2": [], "val_r2": [], "lr": []},
            "callbacks": []
        },
        "completion": {},
        "predictions": {},
        "pipeline": {
            "scaler": str(data['scaler']),
            "selector": f"SelectKBest (k={len(data.get('feature_names', []))})"
        }
    }

    # 2. Augmentation Visuals (Capture before training)
    print("\nStep 2: Capturing Augmentation Data...")
    
    # Label Smoothing
    print("  - Applying Label Smoothing...")
    y_train_smooth = smooth_labels(y_train, epsilon=0.05)
    
    # Sample 100 points for visualization
    sample_idx = np.random.choice(len(y_train), 100, replace=False)
    dashboard_data["augmentation"]["smoothing"] = {
        "original": y_train[sample_idx].tolist(),
        "smoothed": y_train_smooth[sample_idx].tolist()
    }
    
    # Noise Injection (Sample)
    print("  - Applying Feature Noise (Sample)...")
    X_sample = X_train[:50] # Take first 50 samples
    X_sample_noisy = add_noise_to_features(X_sample, noise_level=0.02)
    
    dashboard_data["augmentation"]["noise"] = {
        "original_sample": X_sample[0].tolist()[:20], # Top 20 features of sample 0
        "noisy_sample": X_sample_noisy[0].tolist()[:20]
    }

    # 3. Training Prep
    print("\nStep 3: Training Setup...")
    
    # Custom Callback to populate dashboard_data
    class DashboardCallback(tf.keras.callbacks.Callback):
        def on_epoch_end(self, epoch, logs=None):
            dashboard_data["training"]["history"]["loss"].append(float(logs.get("loss")))
            dashboard_data["training"]["history"]["val_loss"].append(float(logs.get("val_loss")))
            dashboard_data["training"]["history"]["train_r2"].append(float(logs.get("train_r2_sklearn", 0)))
            dashboard_data["training"]["history"]["val_r2"].append(float(logs.get("val_r2_sklearn", 0)))
            
            # Get current LR
            lr = self.model.optimizer.learning_rate
            if isinstance(lr, tf.keras.optimizers.schedules.LearningRateSchedule):
                current_lr = float(lr(self.model.optimizer.iterations))
            else:
                current_lr = float(lr.numpy())
            dashboard_data["training"]["history"]["lr"].append(current_lr)

    # 4. Single Model Training
    print("\n" + "="*80)
    print("Step 4: Training Single Hybrid Model")
    print("="*80)
    
    # Use Model 3's seed
    tf.random.set_seed(44)
    np.random.seed(44)
    
    print("\n🔨 Training model...")
    model = build_hybrid_model(X_train.shape[1])
    
    early_stop = EarlyStopping(monitor='val_loss', patience=40, verbose=0, mode='min', restore_best_weights=True)
    reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=15, min_lr=1e-7, verbose=0)
    metrics_callback = FullMetricsCallback(X_train, y_train, X_test, y_test)
    dash_callback = DashboardCallback()
    
    # Noise Injection during training (Full Set)
    X_train_noisy = add_noise_to_features(X_train, noise_level=0.02)
    
    import time
    start_time = time.time()
    
    history = model.fit(
        X_train_noisy, y_train_smooth,
        epochs=200,
        batch_size=256,
        shuffle=True,
        verbose=0,
        validation_data=(X_test, y_test),
        callbacks=[early_stop, reduce_lr, metrics_callback, dash_callback]
    )
    
    end_time = time.time()
    duration = end_time - start_time
    
    # Evaluate
    y_pred = model.predict(X_test, verbose=0)
    r2 = sklearn.metrics.r2_score(y_test, y_pred)
    print(f"  ✓ Model Validation R²: {r2:.4f}")
    
    # Save model
    model_path = os.path.join(OUTPUT_DIR, 'hybrid_model.keras')
    model.save(model_path)
    
    # Save training history
    history_path = os.path.join(OUTPUT_DIR, 'training_history.pkl')
    with open(history_path, 'wb') as f:
        pickle.dump(history.history, f)
    print(f"  ✓ Saved model and history: {model_path}")

    # Save preprocessing pipeline (scaler, selector, feature names)
    pipeline_path = os.path.join(OUTPUT_DIR, 'hybrid_ultimate_pipeline.pkl')
    with open(pipeline_path, 'wb') as f:
        pickle.dump({
            'scaler': data['scaler'],
            'selector': data['selector'],
            'feature_names': data.get('feature_names', [])
        }, f)
    print(f"  ✓ Saved pipeline: {pipeline_path}")

    # 5. Final Evaluation & Data Export
    print("\nStep 5: Final Evaluation & Dashboard Export...")
    x_pred_train = model.predict(X_train, verbose=0).flatten()
    y_pred_test = model.predict(X_test, verbose=0).flatten()
    
    train_r2 = sklearn.metrics.r2_score(y_train, x_pred_train)
    test_r2 = sklearn.metrics.r2_score(y_test, y_pred_test)
    train_rmse = math.sqrt(sklearn.metrics.mean_squared_error(y_train, x_pred_train))
    test_rmse = math.sqrt(sklearn.metrics.mean_squared_error(y_test, y_pred_test))
    
    print("\n" + "="*80)
    print("FINAL MODEL RESULTS")
    print("="*80)
    print(f"✓ Train R²: {train_r2:.4f} (RMSE: {train_rmse:.4f})")
    print(f"✓ Test R²:  {test_r2:.4f} (RMSE: {test_rmse:.4f})")
    
    # Populate Completion Data
    dashboard_data["completion"] = {
        "duration_sec": duration,
        "epochs_run": len(history.epoch),
        "final_metrics": {
            "train_r2": train_r2, "test_r2": test_r2,
            "train_rmse": train_rmse, "test_rmse": test_rmse
        }
    }
    
    dashboard_data["predictions"] = {
        "train": {
            "actual": y_train.tolist(),
            "pred": x_pred_train.tolist(),
            "smiles": smiles_train.tolist() if len(smiles_train) == len(y_train) else []
        },
        "test": {
            "actual": y_test.tolist(),
            "pred": y_pred_test.tolist(),
            "smiles": smiles_test.tolist() if len(smiles_test) == len(y_test) else []
        }
    }
    
    json_path = os.path.join(OUTPUT_DIR, "data.json")
    with open(json_path, "w") as f:
        json.dump(dashboard_data, f)
    print(f"✓ Saved Dashboard Data to: {json_path}")

    print("="*80)
    print("NODE 2 COMPLETION SUCCESSFUL")
    print("="*80)
