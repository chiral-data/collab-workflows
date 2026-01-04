"""HTML Generator for Node 07: Metabolic Network Clustering
Generates interactive HTML visualization: Single-Panel Dashboard
"""

import json
import os

def generate_clustering_html(json_filename='clustering_data.json'):
    """Generate HTML for Metabolic Network Clustering visualization (Single Panel)"""
    
    # Determine script directory to locate output folder correctly
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    
    # Read JSON data to embed directly
    # Try absolute path first
    json_path = os.path.join(SCRIPT_DIR, 'outputs', os.path.basename(json_filename))
    
    if not os.path.exists(json_path):
        # Fallback to relative path
        json_path = os.path.join('outputs', json_filename)
        if not os.path.exists(json_path):
            json_path = json_filename
        
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            json_content = f.read()
    except Exception as e:
        print(f"Warning: Could not read JSON file {json_path} for embedding: {e}")
        json_content = '{}'

    html = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Metabolic Network Clustering</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.27.0/plotly.min.js" crossorigin="anonymous"></script>
    <style>
        :root {
            --primary: #2c3e50;
            --secondary: #667eea;
            --bg-light: #f5f5f5;
            --text-muted: #7f8c8d;
        }
        
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            margin: 0;
            padding: 20px;
            min-height: 100vh;
        }
        
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 12px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        }
        
        header {
            text-align: center;
            margin-bottom: 30px;
            border-bottom: 2px solid #eee;
            padding-bottom: 20px;
        }
        
        h1 { margin: 0; color: var(--primary); }
        .subtitle { color: var(--text-muted); margin-top: 5px; }
        
        .controls {
            display: flex;
            justify-content: center;
            align-items: center;
            gap: 20px;
            margin-bottom: 30px;
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
        }
        
        label {
            font-weight: 600;
            color: var(--primary);
            font-size: 1.1em;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        
        select {
            padding: 10px 20px;
            border-radius: 6px;
            border: 1px solid #ddd;
            font-size: 18px;
            min-width: 300px;
            background: white;
            cursor: pointer;
            transition: all 0.2s;
            color: var(--primary);
            font-weight: 500;
        }
        
        select:hover { border-color: var(--secondary); }
        select:focus { outline: none; border-color: var(--secondary); box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.2); }
        
        .plot-card {
            background: white;
            border-radius: 10px;
            padding: 20px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.05);
            border: 1px solid #eee;
            min-height: 900px;
            display: flex;
            flex-direction: column;
            align-items: center;
        }
        
        .card-header {
            text-align: center;
            margin-bottom: 20px;
            font-size: 1.4em;
            font-weight: 600;
            color: var(--primary);
        }

        #main-plot {
            width: 100%;
        }

    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🔬 Metabolic Network Clustering</h1>
            <p class="subtitle">Interactive Clustered Correlation Analysis</p>
        </header>

        <div class="controls">
            <label for="cohort-select">Select Dataset:</label>
            <select id="cohort-select" onchange="updateView()"></select>
        </div>

        <div class="plot-card">
            <div class="card-header" id="plot-title">Loading...</div>
            <div id="main-plot"></div>
        </div>
        
    </div>

    <script>
        // Data injected from Python
        const DATA = __JSON_DATA__;
        
        // Configuration
        const COHORTS = Object.keys(DATA.cohorts).sort().reverse(); // Ensure MS+MG (usually M..) comes before Control if possible, or just explicit sort
        
        // Colorscale (YlGnBu matching static fig)
        const COLORSCALE = [
            [0.0, '#ffffd9'],
            [0.2, '#c7e9b4'],
            [0.4, '#41b6c4'],
            [0.6, '#1d91c0'],
            [0.8, '#225ea8'],
            [1.0, '#0c2c84']
        ];

        function init() {
            const select = document.getElementById('cohort-select');
            
            // Populate dropdowns
            COHORTS.forEach(c => {
                select.add(new Option(c, c));
            });

            // Set Defaults
            if (COHORTS.includes('MS+MG')) select.value = 'MS+MG';
            
            updateView();
        }

        function updateView() {
            const val = document.getElementById('cohort-select').value;
            renderHeatmap('main-plot', val);
            document.getElementById('plot-title').textContent = val + " Clustered Matrix";
        }

        function renderHeatmap(divId, cohortName) {
            const data = DATA.cohorts[cohortName];
            if (!data) return;

            const traces = [];
            const layout = {
                height: 850,
                width: 1000, // Fixed width for better centering
                showlegend: false,
                margin: { t: 50, r: 50, b: 50, l: 50 },
                xaxis: { 
                    tickangle: -45, 
                    tickfont: { size: 10 },
                    domain: [0.15, 1],
                    anchor: 'y'
                },
                yaxis: { 
                    automargin: true, 
                    autorange: 'reversed', 
                    tickfont: { size: 10 },
                    domain: [0, 0.85],
                    anchor: 'x',
                    side: 'right'
                },
                // Axes for dendrograms
                xaxis2: {
                    domain: [0.15, 1],
                    anchor: 'y2',
                    showgrid: false,
                    zeroline: false,
                    showticklabels: false
                },
                yaxis2: {
                    domain: [0.85, 1],
                    anchor: 'x2',
                    showgrid: false,
                    zeroline: false,
                    showticklabels: false
                },
                xaxis3: {
                    domain: [0, 0.15],
                    anchor: 'y3',
                    showgrid: false,
                    zeroline: false,
                    showticklabels: false,
                    autorange: 'reversed'
                },
                yaxis3: {
                    domain: [0, 0.85],
                    anchor: 'x3',
                    showgrid: false,
                    zeroline: false,
                    showticklabels: false,
                    autorange: 'reversed'
                }
            };
            
            const config = {
                responsive: true,
                displayModeBar: true,
                modeBarButtons: [['toImage', 'zoom2d', 'pan2d', 'resetScale2d']],
                displaylogo: false
            };

            // 1. Main Heatmap
            traces.push({
                z: data.z,
                x: data.x,
                y: data.y,
                type: 'heatmap',
                colorscale: COLORSCALE,
                zmin: -0.2,
                zmax: 1.0,
                xaxis: 'x',
                yaxis: 'y',
                xgap: 1,
                ygap: 1,
                hovertemplate: '<b>%{x}</b><br><b>%{y}</b><br>Corr: %{z:.2f}<extra></extra>',
                colorbar: {
                    len: 0.2,
                    y: 1.0,
                    yanchor: 'top',
                    x: -0.15,
                    xanchor: 'right',
                    title: 'Corr',
                    tickmode: 'array',
                    tickvals: [-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
                }
            });

            // 2. Dendrograms (if available)
            if (data.dendrogram) {
                // Top Dendrogram (Columns) -> xaxis2, yaxis2
                if (data.dendrogram.col) {
                    const d = data.dendrogram.col;
                    for (let i = 0; i < d.icoords.length; i++) {
                        traces.push({
                            x: d.icoords[i],
                            y: d.dcoords[i],
                            mode: 'lines',
                            type: 'scatter',
                            xaxis: 'x2',
                            yaxis: 'y2',
                            line: { color: 'black', width: 1 },
                            hoverinfo: 'none'
                        });
                    }
                }

                // Left Dendrogram (Rows) -> xaxis3, yaxis3
                if (data.dendrogram.row) {
                    const d = data.dendrogram.row;
                    for (let i = 0; i < d.icoords.length; i++) {
                        // Swap for vertical dendrogram
                        traces.push({
                            x: d.dcoords[i], // Height
                            y: d.icoords[i], // Leaves
                            mode: 'lines',
                            type: 'scatter',
                            xaxis: 'x3',
                            yaxis: 'y3',
                            line: { color: 'black', width: 1 },
                            hoverinfo: 'none'
                        });
                    }
                }
            }

            Plotly.newPlot(divId, traces, layout, config);
        }

        // Run
        document.addEventListener('DOMContentLoaded', init);

    </script>
</body>
</html>'''

    html = html.replace('__JSON_DATA__', json_content)
    return html
