"""
QSAR Prediction Web Application - Flask Backend
Provides API endpoints for real-time SMILES prediction
"""

from flask import Flask, jsonify, request, send_from_directory
from flask_cors import CORS
import os
import json
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, MACCSkeys, Draw
from rdkit.ML.Descriptors import MoleculeDescriptors
import tensorflow as tf
from tensorflow import keras
import pickle
import io
import base64

# Disable RDKit warnings
RDLogger.DisableLog('rdApp.*')

app = Flask(__name__, static_folder='static', static_url_path='')
CORS(app)

# Global variables for loaded model
model = None
scaler = None
ad_stats = None
all_predictions = []

def load_model_components():
    """Load model, scaler, and AD stats on server startup"""
    global model, scaler, ad_stats
    
    print("Loading model components...")
    
    # Load trained model
    model_path = "outputs/model.h5"
    if os.path.exists(model_path):
        model = keras.models.load_model(model_path, compile=False)
        model.compile(loss='mean_squared_error', optimizer='adam')
        print(f"✓ Model loaded from {model_path}")
    else:
        print(f"ERROR: Model not found at {model_path}")
        return False
    
    # Load scaler
    scaler_path = "outputs/scaler.pkl"
    if os.path.exists(scaler_path):
        with open(scaler_path, 'rb') as f:
            scaler = pickle.load(f)
        print(f"✓ Scaler loaded from {scaler_path}")
    else:
        print(f"ERROR: Scaler not found at {scaler_path}")
        return False
    
    # Load AD stats
    ad_path = "outputs/ad_stats.json"
    if os.path.exists(ad_path):
        with open(ad_path, 'r') as f:
            ad_stats = json.load(f)
        print(f"✓ AD stats loaded from {ad_path}")
    else:
        print(f"ERROR: AD stats not found at {ad_path}")
        return False
    
    return True

def load_existing_predictions():
    """Load existing predictions from data.json"""
    global all_predictions
    
    data_path = "outputs/data.json"
    if os.path.exists(data_path):
        with open(data_path, 'r') as f:
            all_predictions = json.load(f)
        print(f"✓ Loaded {len(all_predictions)} existing predictions")
    else:
        print("No existing predictions found")
        all_predictions = []

def compute_descriptors(smiles_str):
    """Compute RDKit descriptors for a SMILES string"""
    mol = Chem.MolFromSmiles(smiles_str)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles_str}")
    
    # Compute RDKit 2D descriptors
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
    rdkit_desc = calc.CalcDescriptors(mol)
    
    # Compute MACCS fingerprints
    maccs = list(MACCSkeys.GenMACCSKeys(mol).ToBitString())
    
    # Combine into feature vector (skip first column which is SMILES)
    features = list(rdkit_desc) + [int(x) for x in maccs]
    
    return np.array(features).reshape(1, -1)

def check_applicability_domain(features, ad_stats):
    """Check if compound is within applicability domain"""
    mins = np.array(ad_stats['min'])
    maxs = np.array(ad_stats['max'])
    
    # Check if all features are within training range
    within_domain = np.all((features >= mins) & (features <= maxs))
    
    return 'IN' if within_domain else 'OUT'

def generate_structure_image(smiles_str):
    """Generate molecular structure image as base64"""
    mol = Chem.MolFromSmiles(smiles_str)
    if mol is None:
        return None
    
    img = Draw.MolToImage(mol, size=(300, 300))
    buf = io.BytesIO()
    img.save(buf, format='PNG')
    buf.seek(0)
    img_str = base64.b64encode(buf.getvalue()).decode('utf-8')
    
    return f"data:image/png;base64,{img_str}"

@app.route('/')
def index():
    """Serve the main dashboard"""
    return send_from_directory('static', 'index.html')

@app.route('/api/compounds', methods=['GET'])
def get_compounds():
    """Get all existing predictions"""
    in_count = sum(1 for p in all_predictions if p['ad_status'] == 'IN')
    out_count = len(all_predictions) - in_count
    
    return jsonify({
        'success': True,
        'compounds': all_predictions,
        'total': len(all_predictions),
        'in_domain': in_count,
        'out_domain': out_count
    })

@app.route('/api/predict', methods=['POST'])
def predict():
    """
    Predict binding affinity for a new SMILES string
    
    Request JSON:
    {
        "smiles": "CCO",
        "drug_name": "Ethanol"  // optional
    }
    
    Response JSON:
    {
        "success": true,
        "drug_name": "Ethanol",
        "smiles": "CCO",
        "prediction": -3.8327,
        "ad_status": "IN",
        "image": "data:image/png;base64,..."
    }
    """
    try:
        data = request.get_json()
        
        if not data or 'smiles' not in data:
            return jsonify({
                'success': False,
                'error': 'Missing SMILES string in request'
            }), 400
        
        smiles_str = data['smiles'].strip()
        # Don't auto-generate "Compound X" name - leave empty/None if not provided
        drug_name = data.get('drug_name', '').strip()
        
        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            return jsonify({
                'success': False,
                'error': f'Invalid SMILES: {smiles_str}'
            }), 400
        
        # Compute descriptors
        features = compute_descriptors(smiles_str)
        
        # Scale features
        features_scaled = scaler.transform(features)
        
        # Make prediction
        prediction = model.predict(features_scaled, verbose=0)[0][0]
        
        # Check applicability domain
        ad_status = check_applicability_domain(features, ad_stats)
        
        # Generate structure image
        image = generate_structure_image(smiles_str)
        
        result = {
            'success': True,
            'drug_name': drug_name,
            'smiles': smiles_str,
            'prediction': float(prediction),
            'ad_status': ad_status,
            'image': image
        }
        
        return jsonify(result)
        
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/health', methods=['GET'])
def health():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy',
        'model_loaded': model is not None,
        'scaler_loaded': scaler is not None,
        'predictions': len(all_predictions)
    })

@app.route('/api/shutdown', methods=['POST'])
def shutdown():
    """Shuts down the server (useful for closing the app from the browser)"""
    def kill_server():
        # Give time for the response to be sent
        import time, os, signal
        time.sleep(1)
        os.kill(os.getpid(), signal.SIGINT)
    
    # Run shutdown in a separate thread so we can return a response 200 first
    import threading
    threading.Thread(target=kill_server).start()
    return jsonify({'success': True, 'message': 'Server is shutting down...'})

def start_server(host='0.0.0.0', port=5000):
    """Start the Flask development server"""
    print("=" * 60)
    print("🚀 Starting QSAR Prediction Web Application")
    print("=" * 60)
    
    # Load model components
    if not load_model_components():
        print("ERROR: Failed to load model components")
        return
    
    # Load existing predictions
    load_existing_predictions()
    
    print("")
    print("✓ Server ready!")
    print(f"✓ Access dashboard at: http://localhost:{port}")
    print(f"✓ API endpoint: http://localhost:{port}/api/predict")
    print(f"✓ Total predictions: {len(all_predictions)}")
    print("=" * 60)
    print("")
    
    # Start server
    app.run(host=host, port=port, debug=False)

if __name__ == '__main__':
    start_server()
