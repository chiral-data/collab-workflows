import json
import os
import sys

def generate_report():
    print("Generating Refined Dashboard (Stages 1, 2, 4, 8) for Node 2...")
    json_path = "outputs/data.json"
    if not os.path.exists(json_path):
        print(f"Error: {json_path} not found. Run 3_model_train.py first.")
        return

    with open(json_path, "r") as f:
        data = json.load(f)

    # ----------------------------------------
    # Unpack Data
    # ----------------------------------------
    setup = data.get("setup", {})
    aug = data.get("augmentation", {})
    training = data.get("training", {})
    # completion = data.get("completion", {}) # Removed Stage 5
    # preds = data.get("predictions", {})     # Removed Stage 6
    pipeline_meta = data.get("pipeline", {})
    
    # Safely get units for Mermaid
    try:
        input_units = setup.get('s1_architecture', {}).get('layers', [])[0].get('units', '?')
    except:
        input_units = '?'

    # Safely get pipeline strings and escape quotes
    scaler_str = str(pipeline_meta.get('scaler', 'StandardScaler')).replace('"', "'")
    selector_str = str(pipeline_meta.get('selector', 'SelectKBest')).replace('"', "'")

    # ----------------------------------------
    # HTML Template
    # ----------------------------------------
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>QSAR Node 2: Enhanced Model Dashboard</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.min.js"></script>
        <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
        <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">
        <style>
            :root {{
                --primary: #1e3a8a;       
                --primary-dark: #172554;  
                --accent: #059669;        
                --bg-body: #f3f4f6;       
                --bg-card: #ffffff;
                --text-main: #1f2937;
                --text-muted: #6b7280;
                --border: #e5e7eb;
                --shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
            }}

            body {{ 
                font-family: 'Inter', sans-serif; 
                background: var(--bg-body); 
                color: var(--text-main); 
                margin: 0; 
                padding: 0; 
                height: 100vh; 
                display: flex; 
                flex-direction: column; 
                overflow: hidden;
            }}

            /* NAVBAR */
            .navbar {{ background: var(--bg-card); border-bottom: 1px solid var(--border); z-index: 100; }}
            .nav-top {{ padding: 15px 30px; display: flex; justify-content: space-between; align-items: center; border-bottom: 1px solid var(--border); }}
            .brand h1 {{ margin: 0; font-size: 1.25rem; font-weight: 700; color: var(--primary); }}
            .brand span {{ background: var(--accent); color: white; padding: 2px 8px; border-radius: 4px; font-size: 0.75rem; }}
            
            .nav-links {{ display: flex; gap: 0px; background: #fff; overflow-x: auto; }}
            .nav-item {{
                padding: 15px 20px;
                cursor: pointer;
                font-size: 0.85rem;
                font-weight: 500;
                color: var(--text-muted);
                border-bottom: 3px solid transparent;
                white-space: nowrap;
                transition: all 0.2s;
            }}
            .nav-item:hover {{ color: var(--primary); background: #f9fafb; }}
            .nav-item.active {{ color: var(--primary); border-bottom-color: var(--primary); font-weight: 600; }}
            .nav-item i {{ margin-right: 6px; color: var(--accent); }}

            /* CONTENT */
            .main-content {{ flex: 1; overflow-y: auto; padding: 30px; max-width: 1600px; width: 100%; margin: 0 auto; box-sizing: border-box; }}
            
            .stage-view {{ display: none; animation: fadeIn 0.4s ease-out; }}
            .stage-view.active {{ display: block; }}
            @keyframes fadeIn {{ from {{ opacity: 0; transform: translateY(5px); }} to {{ opacity: 1; transform: translateY(0); }} }}

            /* CARDS */
            .card {{ background: var(--bg-card); border: 1px solid var(--border); border-radius: 12px; box-shadow: var(--shadow); padding: 25px; margin-bottom: 30px; }}
            .card h3 {{ margin-top: 0; font-size: 1.1rem; color: var(--primary); border-bottom: 2px solid #f0fdf4; padding-bottom: 10px; display: inline-block; }}
            
            .plot-box {{ width: 100%; height: 500px; border-radius: 8px; overflow: hidden; }}
            .grid-2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 25px; }}
            
            /* TABLES */
            table {{ width: 100%; border-collapse: collapse; margin: 15px 0; }}
            th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid var(--border); }}
            th {{ background: #f9fafb; font-weight: 600; font-size: 0.85rem; color: var(--text-muted); }}
            
            /* MERMAID */
            .mermaid {{ display: flex; justify-content: center; background: #fafafa; padding: 20px; border-radius: 8px; }}
        </style>
    </head>
    <body>

    <div class="navbar">
        <div class="nav-top">
            <div class="brand">
                <h1>QSAR Node 2 <span>Training Dashboard</span></h1>
            </div>
            <div style="font-size: 0.8rem; color: #9ca3af;">Ultimate Hybrid Model</div>
        </div>
        <div class="nav-links">
            <div class="nav-item active" onclick="showStage('s1')"><i class="fas fa-layer-group"></i> 1. Setup & Arch</div>
            <div class="nav-item" onclick="showStage('s2')"><i class="fas fa-magic"></i> 2. Augmentation</div>
            <div class="nav-item" onclick="showStage('s4')"><i class="fas fa-history"></i> 3. Callbacks</div>
            <div class="nav-item" onclick="showStage('s8')"><i class="fas fa-project-diagram"></i> 4. Pipeline</div>
        </div>
    </div>

    <div class="main-content">
        
        <!-- STAGE 1: SETUP & ARCHITECTURE -->
        <div id="s1" class="stage-view active">
            <h2>Stage 1: Model Architecture & Pre-Training Setup</h2>
            <div class="grid-2">
                <div class="card">
                    <h3>Architecture Blueprint</h3>
                    <div class="mermaid">
                        graph TD
                        I["Input ({input_units} units)"] --> L1["Dense 192 <br/>(ReLU, L2=0.002)"]
                        L1 --> BN1["Batch Norm"]
                        BN1 --> D1["Dropout 0.4"]
                        D1 --> L2["Dense 96 <br/>(ReLU, L2=0.002)"]
                        L2 --> BN2["Batch Norm"]
                        BN2 --> D2["Dropout 0.3"]
                        D2 --> L3["Dense 48 <br/>(ReLU, L2=0.001)"]
                        L3 --> BN3["Batch Norm"]
                        BN3 --> D3["Dropout 0.2"]
                        D3 --> O["Output <br/>(Linear)"]
                        
                        style I fill:#e0e7ff,stroke:#3730a3
                        style O fill:#dcfce7,stroke:#166534
                    </div>
                </div>
                <div class="card">
                    <h3>Hyperparameter Configuration</h3>
                    <table>
                        <tr><th>Parameter</th><th>Value</th></tr>
                        <tr><td>Optimizer</td><td>{setup.get('s1_hyperparams', {}).get('optimizer', 'Adam')}</td></tr>
                        <tr><td>Learning Rate</td><td>{setup.get('s1_hyperparams', {}).get('learning_rate', 0.001)}</td></tr>
                        <tr><td>Batch Size</td><td>{setup.get('s1_hyperparams', {}).get('batch_size', 256)}</td></tr>
                        <tr><td>Epochs</td><td>{setup.get('s1_hyperparams', {}).get('epochs', 200)}</td></tr>
                        <tr><td>Regularizers</td><td>L2 + BatchNormalization + Dropout</td></tr>
                    </table>
                </div>
            </div>
        </div>

        <!-- STAGE 2: DATA AUGMENTATION -->
        <div id="s2" class="stage-view">
            <h2>Stage 2: Data Augmentation Effects</h2>
            <div class="grid-2">
                <div class="card">
                    <h3>Label Smoothing (Refinement)</h3>
                    <p>Comparison of Original Targets vs Smoothed Targets (epsilon=0.05)</p>
                    <div id="smoothPlot" class="plot-box"></div>
                </div>
                <div class="card">
                    <h3>Feature Noise Injection</h3>
                    <p>Comparison of Original Feature Vector vs Noised Version (std=0.02)</p>
                    <div id="noisePlot" class="plot-box"></div>
                </div>
            </div>
        </div>

        <!-- STAGE 4: CALLBACK ACTIONS (Renumbered to 3 in UI) -->
        <div id="s4" class="stage-view">
            <h2>Stage 3: Callback Actions & Events</h2>
            <div class="card">
                <h3>Event Timeline</h3>
                <p>Visualization of ModelCheckpoint, ReduceLROnPlateau, and EarlyStopping events.</p>
                <div id="timelinePlot" class="plot-box" style="height: 300px;"></div>
            </div>
        </div>

        <!-- STAGE 8: PIPELINE METADATA (Renumbered to 4 in UI) -->
        <div id="s8" class="stage-view">
             <h2>Stage 4: Pipeline Metadata</h2>
             <div class="card">
                 <h3>Preprocessing Flow</h3>
                 <div class="mermaid">
                    graph LR
                    Raw["Raw Input"] --> Scaler["{scaler_str}"]
                    Scaler --> Select["{selector_str}"]
                    Select --> Model["Hybrid Model"]
                    
                    style Raw fill:#f3f4f6,stroke:#9ca3af
                    style Scaler fill:#dbeafe,stroke:#1e40af
                    style Select fill:#dbeafe,stroke:#1e40af
                    style Model fill:#dcfce7,stroke:#166534
                 </div>
             </div>
        </div>

    </div>

    <script>
        mermaid.initialize({{ startOnLoad: false, theme: 'neutral' }});
        
        // DATA
        const aug = {json.dumps(aug)};
        const history = {json.dumps(training.get('history', {}))};
        
        // Render Mermaid on Load
        document.addEventListener('DOMContentLoaded', () => {{
            mermaid.run({{ querySelector: '.stage-view.active .mermaid' }});
        }});

        // NAV LOGIC
        function showStage(id) {{
            document.querySelectorAll('.nav-item').forEach(el => el.classList.remove('active'));
            event.currentTarget.classList.add('active');
            document.querySelectorAll('.stage-view').forEach(el => el.classList.remove('active'));
            document.getElementById(id).classList.add('active');
            
            // Re-render Mermaid for the newly visible stage
            setTimeout(() => {{
                const mermaidContainer = document.querySelector('#' + id + ' .mermaid');
                if (mermaidContainer && !mermaidContainer.getAttribute('data-processed')) {{
                    mermaid.run({{ nodes: [mermaidContainer] }});
                }}
                window.dispatchEvent(new Event('resize'));
            }}, 50);
        }}

        // --- PLOTS --- //

        // 1. Augmentation (Smoothing)
        if (document.getElementById('smoothPlot')) {{
            Plotly.newPlot('smoothPlot', [
                {{ x: aug.smoothing.original, type: 'histogram', name: 'Original', opacity: 0.6 }},
                {{ x: aug.smoothing.smoothed, type: 'histogram', name: 'Smoothed', opacity: 0.6 }}
            ], {{ barmode: 'overlay', title: 'Label Smoothing Distribution' }});
        }}

        // 2. Augmentation (Noise)
        if (document.getElementById('noisePlot')) {{
            Plotly.newPlot('noisePlot', [
                {{ y: aug.noise.original_sample, type: 'bar', name: 'Original Feature Vector' }},
                {{ y: aug.noise.noisy_sample, type: 'bar', name: 'Noised Vector' }}
            ], {{ title: 'Feature Noise Injection Effect (Sample)' }});
        }}

        // 3. Timeline
        if (document.getElementById('timelinePlot') && history.loss) {{
            const epochs = history.loss.map((_, i) => i + 1);
            const lrDrops = [];
            for(let i=1; i<history.lr.length; i++) {{
                if(history.lr[i] < history.lr[i-1]) lrDrops.push(i);
            }}
            
            Plotly.newPlot('timelinePlot', [
            {{ x: [0, epochs.length], y: [1, 1], mode: 'lines', line: {{width: 4, color:'#ddd'}}, showlegend:false }},
            {{ x: lrDrops, y: Array(lrDrops.length).fill(1), mode: 'markers', marker: {{size: 15, color: 'orange'}}, name: 'LR Reduction' }},
            {{ x: [epochs.length], y: [1], mode: 'markers', marker: {{size: 15, color: 'red', symbol: 'x'}}, name: 'Stop' }}
            ], {{ 
                yaxis: {{ showgrid: false, zeroline: false, showticklabels: false }},
                xaxis: {{ title: 'Epochs' }},
                height: 200
            }});
        }}
    </script>
    </body>
    </html>
    """
    
    output_path = "outputs/report.html"
    with open(output_path, "w") as f:
        f.write(html_content)
    print(f"Refined Dashboard generated at {output_path}")

if __name__ == "__main__":
    generate_report()