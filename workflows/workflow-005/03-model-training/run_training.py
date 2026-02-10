import numpy as np
import os
import json
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras import backend as K
import sklearn.metrics

# Reproducibility
tf.random.set_seed(1)
np.random.seed(3)

def rmse(y_true, y_pred):
    return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1))

def r_square(y_true, y_pred):
    SS_res = K.sum(K.square(y_true - y_pred))
    SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
    return (1 - SS_res / (SS_tot + K.epsilon()))

def find_file(filename, search_path):
    for root, dirs, files in os.walk(search_path):
        if filename in files:
            return os.path.join(root, filename)
    return None

import re

def get_param(param_name, default_if_missing=None):
    # 1. Try Environment Variable
    val = os.environ.get(param_name.upper()) or os.environ.get(f"PARAMS_{param_name.upper()}")
    if val is not None:
        return val

    # 2. Try Parsing .chiral/job.toml
    try:
        toml_path = os.path.join(".chiral", "job.toml")
        if os.path.exists(toml_path):
            with open(toml_path, "r") as f:
                content = f.read()
            section_match = re.search(f'\\[params\\.{param_name}\\]', content, re.IGNORECASE)
            if section_match:
                start_index = section_match.end()
                default_match = re.search(r'default\s*=\s*(.+)', content[start_index:])
                if default_match:
                    raw_val = default_match.group(1).strip()
                    if raw_val.startswith('"') and raw_val.endswith('"'):
                         return raw_val[1:-1].replace('\\"', '"')
                    elif raw_val.startswith("'") and raw_val.endswith("'"):
                         return raw_val[1:-1]
                    else:
                         return raw_val
    except Exception as e:
        print(f"Warning: Failed to parse job.toml for {param_name}: {e}")
        
    if default_if_missing is not None:
        return default_if_missing
        
    raise ValueError(f"Parameter {param_name} missing (Env Var and job.toml default)")

def run():
    print("--- Node 3: Model Training ---")
    os.makedirs("outputs", exist_ok=True)

    # Silva shared workspace: read directly from previous node's outputs
    input_file = os.path.join("..", "02-feature-engineering", "outputs", "processed_data.npz")
    
    if not os.path.exists(input_file):
        print(f"Path {input_file} not found, searching for processed_data.npz...")
        input_file = find_file("processed_data.npz", "..")

    if not input_file:
        print("Processed data not found.")
        # Debug
        for root, dirs, files in os.walk(".."):
            if "processed_data.npz" in files:
                print(f"  Found in: {root}")
        return

    print(f"Found input: {input_file}")

    data = np.load(input_file)
    X_train, y_train = data['X_train'], data['y_train']
    X_test, y_test = data['X_test'], data['y_test']

    # Build Model
    model = Sequential()
    model.add(Dense(600, input_dim=X_train.shape[1], activation='relu'))
    model.add(Dense(100, activation='relu'))
    model.add(Dense(100, activation='relu'))
    model.add(Dense(1, activation='linear'))
    
    model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mae','mape', rmse, r_square])

    early_stop = EarlyStopping(monitor='val_r_square', patience=200, mode='max', restore_best_weights=True)

    epochs = int(get_param("epochs"))
    batch_size = int(get_param("batch_size"))
    history = model.fit(X_train, y_train, epochs=epochs, batch_size=batch_size, 
                        validation_data=(X_test, y_test), callbacks=[early_stop], verbose=1)

    # Save Model
    model.save("outputs/model.h5")
    
    # Evaluation for JSON
    y_pred = model.predict(X_test).flatten()
    r2_score = sklearn.metrics.r2_score(y_test, y_pred)
    
    # Clean numpy types for JSON
    json_data = {
        "params": {"epochs": epochs, "batch_size": batch_size, "optimizer": "adam"},
        "metrics": {
            "r2_score": float(r2_score),
            "mae": float(sklearn.metrics.mean_absolute_error(y_test, y_pred))
        },
        "history": {
            "loss": [float(x) for x in history.history['loss']],
            "val_loss": [float(x) for x in history.history['val_loss']],
            "val_r_square": [float(x) for x in history.history['val_r_square']]
        },
        "predictions": {
            "y_test": [float(x) for x in y_test[:500]], # Limit size for JSON
            "y_pred": [float(x) for x in y_pred[:500]]
        }
    }
    
    with open("outputs/data.json", "w") as f:
        json.dump(json_data, f)
    print("Training Complete.")

if __name__ == "__main__":
    run()