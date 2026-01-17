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

def run():
    print("--- Node 3: Model Training ---")
    os.makedirs("outputs", exist_ok=True)
    
    input_file = "../02_Feature_Engineering/outputs/processed_data.npz"
    if not os.path.exists(input_file):
        print("Processed data not found.")
        return

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

    history = model.fit(X_train, y_train, epochs=200, batch_size=400, 
                        validation_data=(X_test, y_test), callbacks=[early_stop], verbose=1)

    # Save Model
    model.save("outputs/model.h5")
    
    # Evaluation for JSON
    y_pred = model.predict(X_test).flatten()
    r2_score = sklearn.metrics.r2_score(y_test, y_pred)
    
    # Clean numpy types for JSON
    json_data = {
        "params": {"epochs": 200, "batch_size": 400, "optimizer": "adam"},
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