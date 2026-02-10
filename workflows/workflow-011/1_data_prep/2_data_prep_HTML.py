import json
import os
import sys

def generate_report():
    print("Generating Scientific Dashboard for Node 1...")
    json_path = "outputs/data.json"
    if not os.path.exists(json_path):
        print(f"Error: {json_path} not found. Run 1_data_prep.py first.")
        return

    with open(json_path, "r") as f:
        data = json.load(f)

    # ----------------------------------------
    # Unpack Data
    # ----------------------------------------
    stats = data.get('summary', {})
    aug_viz = data.get('augmentation', {})
    sankey = data.get('sankey_counts', {})
    feat_sel = data.get('feature_selection', {})
    pca_data = data.get('pca_analysis', {})
    dist_data = data.get('distributions', {})
    corr_data = data.get('correlation', {})

    # ----------------------------------------
    # HTML Template
    # ----------------------------------------
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>QSAR Node 1: Scientific Dashboard</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
        <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">
        <style>
            :root {{
                /* PREMIUM SCIENTIFIC PALETTE */
                --primary: #1e3a8a;       /* Deep Ocean Blue */
                --primary-dark: #172554;  /* Darker Navy */
                --accent: #059669;        /* Scientific Teal/Green */
                --accent-hover: #047857;
                --bg-body: #f3f4f6;       /* Clean Light Gray */
                --bg-card: #ffffff;
                --text-main: #1f2937;
                --text-muted: #6b7280;
                --border-color: #e5e7eb;
                --shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06);
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

            /* TOP NAVIGATION BAR */
            .navbar {{
                background: var(--bg-card);
                border-bottom: 1px solid var(--border-color);
                box-shadow: 0 2px 4px rgba(0,0,0,0.05);
                display: flex;
                flex-direction: column;
                z-index: 100;
            }}
            .nav-top {{
                padding: 15px 30px;
                display: flex;
                justify-content: space-between;
                align-items: center;
                border-bottom: 1px solid var(--border-color);
            }}
            .brand {{ display: flex; align-items: center; gap: 15px; }}
            .brand h1 {{ margin: 0; font-size: 1.25rem; font-weight: 700; color: var(--primary); letter-spacing: -0.5px; }}
            .brand span {{ background: var(--accent); color: white; padding: 2px 8px; border-radius: 4px; font-size: 0.75rem; font-weight: 600; text-transform: uppercase; }}
            
            .nav-links {{
                display: flex;
                gap: 5px;
                padding: 0 30px;
                background: #fff;
                overflow-x: auto;
            }}
            .nav-item {{
                padding: 15px 20px;
                cursor: pointer;
                font-size: 0.9rem;
                font-weight: 500;
                color: var(--text-muted);
                border-bottom: 3px solid transparent;
                transition: all 0.2s;
                white-space: nowrap;
            }}
            .nav-item:hover {{ color: var(--primary); background: #f9fafb; }}
            .nav-item.active {{ color: var(--primary); border-bottom-color: var(--primary); font-weight: 600; }}
            .nav-item i {{ margin-right: 8px; color: var(--accent); }}

            /* MAIN CONTENT AREA */
            .main-content {{
                flex: 1;
                overflow-y: auto;
                padding: 30px;
                max-width: 1600px;
                width: 100%;
                margin: 0 auto;
                box-sizing: border-box;
            }}

            /* ANIMATIONS */
            .stage-view {{ display: none; animation: slideUp 0.4s ease-out; }}
            .stage-view.active {{ display: block; }}
            @keyframes slideUp {{ from {{ opacity: 0; transform: translateY(10px); }} to {{ opacity: 1; transform: translateY(0); }} }}

            /* HEADERS */
            .section-header {{ margin-bottom: 25px; display: flex; justify-content: space-between; align-items: flex-end; }}
            .section-header h2 {{ margin: 0; font-size: 1.8rem; font-weight: 700; color: var(--primary-dark); }}
            .section-header p {{ margin: 5px 0 0; color: var(--text-muted); font-size: 1rem; }}

            /* CARDS */
            .card {{
                background: var(--bg-card);
                border: 1px solid var(--border-color);
                border-radius: 12px;
                box-shadow: var(--shadow);
                padding: 25px;
                margin-bottom: 30px;
                transition: transform 0.2s, box-shadow 0.2s;
            }}
            .card:hover {{ box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1); }}

            /* CONTROLS BAR */
            .controls-bar {{
                background: white;
                padding: 12px 20px;
                background: linear-gradient(to right, #ffffff, #f9fafb);
                border: 1px solid var(--border-color);
                border-left: 4px solid var(--accent);
                border-radius: 8px;
                display: flex;
                gap: 20px;
                align-items: center;
                margin-bottom: 20px;
                box-shadow: 0 1px 2px rgba(0,0,0,0.05);
            }}
            .control-group {{ display: flex; align-items: center; gap: 10px; }}
            .control-group label {{ font-size: 0.85rem; font-weight: 600; text-transform: uppercase; color: var(--text-muted); letter-spacing: 0.5px; }}
            
            input[type="range"] {{ height: 6px; background: #e5e7eb; border-radius: 5px; appearance: none; }}
            input[type="range"]::-webkit-slider-thumb {{ appearance: none; width: 18px; height: 18px; background: var(--accent); border-radius: 50%; cursor: pointer; border: 2px solid white; box-shadow: 0 0 2px rgba(0,0,0,0.2); }}
            
            input[type="text"] {{ padding: 8px 12px; border: 1px solid var(--border-color); border-radius: 6px; font-family: 'Inter'; font-size: 0.9rem; transition: border 0.2s; }}
            input[type="text"]:focus {{ border-color: var(--accent); outline: none; }}

            /* GRIDS */
            .grid-2 {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(500px, 1fr)); gap: 25px; }}
            
            .plot-box {{ width: 100%; height: 500px; border-radius: 8px; overflow: hidden; }}

            /* AUGMENTATION STYLES */
            .aug-grid {{ display: flex; justify-content: center; align-items: center; gap: 30px; padding: 20px; }}
            .aug-item {{ text-align: center; background: #f8fafc; padding: 15px; border-radius: 8px; border: 1px solid var(--border-color); }}
            .aug-item img {{ mix-blend-mode: multiply; max-width: 200px; }}
            
        </style>
    </head>
    <body>

    <!-- TOP NAVIGATION -->
    <div class="navbar">
        <div class="nav-top">
            <div class="brand">
                <h1>QSAR Pipeline <span>Node 1</span></h1>
                <div style="height: 20px; width: 1px; background: #e5e7eb;"></div>
                <div style="font-size: 0.9rem; color: var(--text-muted);">Scientific Data Preparation Report</div>
            </div>
            <div style="font-size: 0.8rem; color: #9ca3af;">Generated with Plotly & Python</div>
        </div>
        
        <div class="nav-links">
            <div class="nav-item active" onclick="showStage('stage-summary')"><i class="fas fa-home"></i> Summary</div>
            <div class="nav-item" onclick="showStage('stage-1')"><i class="fas fa-flask"></i> Augmentation</div>
            <div class="nav-item" onclick="showStage('stage-2')"><i class="fas fa-dna"></i> Descriptors</div>
            <div class="nav-item" onclick="showStage('stage-3')"><i class="fas fa-filter"></i> Feature Selection</div>
            <div class="nav-item" onclick="showStage('stage-5')"><i class="fas fa-columns"></i> Splits</div>
            <div class="nav-item" onclick="showStage('stage-6')"><i class="fas fa-search-minus"></i> Outliers & PCA</div>
            <div class="nav-item" onclick="showStage('stage-7')"><i class="fas fa-project-diagram"></i> Correlations</div>
        </div>
    </div>

    <!-- MAIN CONTENT -->
    <div class="main-content">

        <!-- SUMMARY -->
        <div id="stage-summary" class="stage-view active">
            <div class="section-header">
                <div>
                    <h2>Pipeline Overview</h2>
                    <p>High-level Data Flow Summary</p>
                </div>
            </div>
            <div class="card">
                <div id="sankeyPlot" class="plot-box"></div>
            </div>
        </div>

        <!-- STAGE 1: AUGMENTATION -->
        <div id="stage-1" class="stage-view">
            <div class="section-header">
                <h2>Data Augmentation</h2>
                <p>Generating synthetic SMILES variants to enhance robustness</p>
            </div>
            <div class="card">
                <div class="aug-grid">
                    <div class="aug-item">
                        <div style="font-weight:600; color:var(--primary); margin-bottom:5px;">Original</div>
                        <img src="{aug_viz.get('original_img', '')}">
                    </div>
                    <div style="font-size: 2rem; color: var(--text-muted);"><i class="fas fa-arrow-right"></i></div>
                    {''.join([f'''
                    <div class="aug-item">
                         <div style="font-weight:600; color:var(--accent); margin-bottom:5px;">Variant {i+1}</div>
                        <img src="{v['img']}">
                    </div>
                    ''' for i, v in enumerate(aug_viz.get('variants', []))])}
                </div>
                <div style="text-align:center; margin-top:15px; font-family:monospace; color:var(--text-muted);">
                    {aug_viz.get('original_smiles', '')}
                </div>
            </div>
        </div>

        <!-- STAGE 2: DESCRIPTORS -->
        <div id="stage-2" class="stage-view">
            <div class="section-header">
                <h2>Descriptor Calculation</h2>
                <p>Feature Engineering Composition</p>
            </div>
            <div class="card" style="display:flex; justify-content:center; align-items:center; flex-direction:column;">
                <div id="descPiePlot" class="plot-box" style=" max-width:800px;"></div>
            </div>
        </div>

        <!-- STAGE 3: FEATURES -->
        <div id="stage-3" class="stage-view">
            <div class="section-header">
                <h2>Feature Selection</h2>
                <p>Random Forest Importance Analysis</p>
            </div>
            
            <div class="controls-bar">
                <div class="control-group">
                    <label>Top N Features:</label>
                    <input type="range" id="featSlider" min="5" max="50" value="20" oninput="updateFeatPlot(this.value)">
                    <span id="featCountVal" style="font-weight:bold; color:var(--primary); min-width:30px;">20</span>
                </div>
            </div>

            <div class="grid-2">
                <div class="card">
                    <h3 style="margin-top:0; color:var(--primary);">Cumulative Importance</h3>
                    <div id="cumulativePlot" style="height:400px;"></div>
                </div>
                <div class="card">
                    <h3 style="margin-top:0; color:var(--primary);">Top Feature Importance</h3>
                    <div id="featPlot" style="height:400px;"></div>
                </div>
            </div>
        </div>

         <!-- STAGE 5: SPLIT -->
        <div id="stage-5" class="stage-view">
            <div class="section-header">
                <h2>Train/Test Split</h2>
                <p>Distribution Stability Check</p>
            </div>
            
            <div class="controls-bar">
                <label style="margin-right:15px;">Overlay Controls:</label>
                <div class="control-group">
                    <input type="checkbox" id="chkTrain" checked onclick="updateSplitPlot()" style="height:16px; width:16px;"> <label for="chkTrain" style="cursor:pointer;">Train Set</label>
                </div>
                <div class="control-group" style="margin-left:15px;">
                    <input type="checkbox" id="chkTest" checked onclick="updateSplitPlot()" style="height:16px; width:16px;"> <label for="chkTest" style="cursor:pointer;">Test Set</label>
                </div>
            </div>

            <div class="card">
                <div id="splitDistPlot" class="plot-box"></div>
            </div>
        </div>

        <!-- STAGE 6: OUTLIERS -->
        <div id="stage-6" class="stage-view">
            <div class="section-header">
                <h2>Outlier Detection</h2>
                <p>Isolation Forest Decision Boundary & PCA</p>
            </div>
            
            <div class="controls-bar">
                <div class="control-group">
                    <label>Max Docking Score:</label>
                    <input type="range" id="scoreRange" min="-15" max="0" step="0.5" value="-12.0" oninput="applyFilters()">
                    <span id="scoreVal" style="font-weight:bold; color:var(--primary); min-width:40px;">-12.0</span>
                </div>
                <div style="width: 1px; height: 24px; background: var(--border-color); margin: 0 15px;"></div>
                <div class="control-group">
                    <label><i class="fas fa-search"></i> SMILES:</label>
                    <input type="text" id="smilesSearch" placeholder="Structure search..." onkeyup="applyFilters()">
                </div>
                <div style="margin-left:auto; font-size:0.85rem; font-weight:600; color:var(--accent);">
                    Displayed: <span id="filterCount">0</span>
                </div>
            </div>

            <div class="grid-2">
                <div class="card">
                    <div id="outlierPcaPlot" class="plot-box"></div>
                </div>
                <div class="card">
                    <div id="variancePlot" class="plot-box"></div>
                </div>
            </div>
        </div>

        <!-- STAGE 7: CORRELATION -->
        <div id="stage-7" class="stage-view">
            <div class="section-header">
                <h2>Correlation Matrix</h2>
                <p>Feature Collinearity Heatmap</p>
            </div>
            <div class="card" style="overflow:hidden;">
                <div id="corrPlot" style="height: 700px;"></div>
            </div>
        </div>

    </div>

    <script>
        // ----------------------------------------------------
        // VARIABLES & DATA
        // ----------------------------------------------------
        const sankeyStats = {json.dumps(sankey)};
        const featSel = {json.dumps(feat_sel)};
        const descStats = {json.dumps(stats)};
        const pcaInliers = {json.dumps(pca_data.get('inliers', {}))};
        const pcaOutliers = {json.dumps(pca_data.get('outliers', {}))};
        const pcaVar = {json.dumps(pca_data.get('variance_ratio', []))};
        const trainY = {json.dumps(dist_data.get('distributions', {}).get('train_labels', dist_data.get('train_labels', [])))};
        const testY = {json.dumps(dist_data.get('distributions', {}).get('test_labels', dist_data.get('test_labels', [])))};
        const corrData = {json.dumps(corr_data)};

        // Colors
        const C_PRIMARY = '#1e3a8a';
        const C_ACCENT = '#059669';
        const C_TRAIN = '#f59e0b';
        const C_TEST = '#3b82f6';
        const C_OUTLIE = '#ef4444';

        // ----------------------------------------------------
        // NAVIGATION
        // ----------------------------------------------------
        function showStage(id) {{
            document.querySelectorAll('.nav-item').forEach(el => el.classList.remove('active'));
            event.currentTarget.classList.add('active');
            
            document.querySelectorAll('.stage-view').forEach(el => el.classList.remove('active'));
            document.getElementById(id).classList.add('active');
            window.dispatchEvent(new Event('resize'));
        }}
        
        // ----------------------------------------------------
        // CHARTS
        // ----------------------------------------------------
        
        // 1. Sankey
        Plotly.newPlot('sankeyPlot', [{{
            type: "sankey", orientation: "h",
            node: {{
                pad: 15, thickness: 20, line: {{ color: "black", width: 0.5 }},
                label: ["Loaded", "Augment", "Features", "Train", "Test", "Train Final", "Test Final", "Outliers"],
                color: ["#60a5fa", "#34d399", "#a78bfa", "#fbbf24", "#818cf8", "#22c55e", "#f472b6", "#f87171"]
            }},
            link: {{
                source: [0, 1, 2, 2, 3, 4, 3],
                target: [1, 2, 3, 4, 5, 6, 7],
                value: [
                    sankeyStats.loaded, sankeyStats.augmented, sankeyStats.split_train, sankeyStats.split_test,
                    sankeyStats.final_train, sankeyStats.final_test, descStats.outliers_removed_train
                ]
            }}
        }}], {{ margin: {{ t: 20, b: 20, l: 20, r: 20 }}, font: {{ family: 'Inter' }} }});

        // 2. Pie
        Plotly.newPlot('descPiePlot', [{{
            labels: ['RDKit 2D', 'MACCS Keys'],
            values: [descStats.desc_rdkit, descStats.desc_maccs],
            type: 'pie', marker: {{ colors: [C_PRIMARY, '#8b5cf6'] }}, hole: 0.6,
            textinfo: 'label+percent', hoverinfo: 'label+value+percent'
        }}], {{ 
            title: {{ text: 'Feature Composition', font: {{ size: 20, family: 'Inter' }} }},
            showlegend: false
        }});

        // 3. Features
        const allFeats = featSel.top_features; 
        function updateFeatPlot(n) {{
            n = parseInt(n);
            document.getElementById('featCountVal').innerText = n;
            const slice = allFeats.slice(0, n);
            
            Plotly.newPlot('featPlot', [{{
                x: slice.map(d => d.name), y: slice.map(d => d.importance),
                type: 'bar', marker: {{ color: C_PRIMARY, opacity: 0.8 }}
            }}], {{ margin: {{ t: 10, b: 100, l: 40, r: 10 }}, xaxis: {{ tickangle: -45 }} }});
        }}
        updateFeatPlot(20);
        
        Plotly.newPlot('cumulativePlot', [{{
            y: featSel.cumulative_importance, type: 'scatter', mode: 'lines', fill: 'tozeroy', 
            line: {{ color: C_ACCENT, width: 3 }}
        }}], {{ margin: {{ t: 20, b: 30, l: 40, r: 20 }}, yaxis: {{range: [0, 1.05]}} }});

        // 5. Splits
        function updateSplitPlot() {{
            const showTrain = document.getElementById('chkTrain').checked;
            const showTest = document.getElementById('chkTest').checked;
            
            const traces = [];
            if(showTrain) traces.push({{ x: trainY, type: 'histogram', name: 'Train Set', opacity: 0.6, marker: {{ color: C_TRAIN }} }});
            if(showTest) traces.push({{ x: testY, type: 'histogram', name: 'Test Set', opacity: 0.6, marker: {{ color: C_TEST }} }});
            
            Plotly.newPlot('splitDistPlot', traces, {{ 
                barmode: 'overlay', margin: {{ t: 20, b: 30 }},
                legend: {{ orientation: 'h', y: 1.1 }}
            }});
        }}
        updateSplitPlot();

        // 6. Outliers
        Plotly.newPlot('variancePlot', [{{
            x: ['PC1', 'PC2'], y: pcaVar, type: 'bar', text: pcaVar.map(v => (v*100).toFixed(1) + '%'),
            textposition: 'auto', marker: {{ color: '#94a3b8' }}
        }}], {{ title: 'Explained Variance Ratio', margin: {{ t: 40, b: 30 }} }});

        function applyFilters() {{
            const scoreThresh = parseFloat(document.getElementById('scoreRange').value);
            const query = document.getElementById('smilesSearch').value.toLowerCase();
            document.getElementById('scoreVal').innerText = scoreThresh;

            const indices = []; 
            for(let i=0; i<pcaInliers.c.length; i++) {{
                const score = pcaInliers.c[i];
                const smi = pcaInliers.s ? pcaInliers.s[i].toLowerCase() : "";
                
                if (score <= scoreThresh && (query === "" || smi.includes(query))) {{
                    indices.push(i);
                }}
            }}
            
            document.getElementById('filterCount').innerText = indices.length + " / " + pcaInliers.c.length;

            const currentX = indices.map(i => pcaInliers.x[i]);
            const currentY = indices.map(i => pcaInliers.y[i]);
            const currentC = indices.map(i => pcaInliers.c[i]);
            const currentS = indices.map(i => pcaInliers.s ? pcaInliers.s[i] : "");

             Plotly.newPlot('outlierPcaPlot', [
                {{
                    x: currentX, y: currentY, mode: 'markers', name: 'Inliers',
                    text: currentS, marker: {{ color: currentC, colorscale: 'Viridis', size: 6, showscale: true }}
                }},
                {{
                    x: pcaOutliers.x, y: pcaOutliers.y, mode: 'markers', name: 'Outliers',
                    marker: {{ color: C_OUTLIE, size: 6, symbol: 'x' }}
                }}
            ], {{
                title: 'PCA Analysis (Train Set)',
                xaxis: {{ title: 'PC1' }}, yaxis: {{ title: 'PC2' }},
                margin: {{ t: 40, b: 30, l: 40, r: 10 }}, hovermode: 'closest',
                showlegend: true
            }});
        }}
        applyFilters();

        // 7. Correlation
         Plotly.newPlot('corrPlot', [{{
            z: corrData.z, x: corrData.x, y: corrData.y, type: 'heatmap', colorscale: 'RdBu', zmin: -1, zmax: 1
        }}], {{ margin: {{ t: 20, b: 100, l: 150, r: 20 }} }});
        
    </script>
    </body>
    </html>
    """

    output_path = "outputs/report.html"
    with open(output_path, "w") as f:
        f.write(html_content)
    print(f"Scientific SPA Dashboard generated at {output_path}")

if __name__ == "__main__":
    generate_report()