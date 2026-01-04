"""HTML Generator for Node 07: Metabolic Network Clustering
Generates interactive HTML visualization: Enhanced Single-Panel Dashboard with Pairwise Search
"""

import json
import os

def generate_clustering_html(json_filename='clustering_data.json'):
    """Generate HTML for Metabolic Network Clustering visualization (Enhanced + Pairwise)"""
    
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
            --highlight-a: #ff00ff; /* Magenta */
            --highlight-b: #00ffff; /* Cyan */
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
            max-width: 1400px;
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
        
        /* Controls Section */
        .controls-area {
            background: #f8f9fa;
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            display: flex;
            flex-wrap: wrap;
            gap: 20px;
            justify-content: center;
            align-items: flex-start;
            box-shadow: inset 0 2px 4px rgba(0,0,0,0.02);
            border: 1px solid #eee;
        }

        .control-group {
            display: flex;
            flex-direction: column;
            gap: 8px;
            min-width: 200px;
        }

        .control-label {
            font-weight: 600;
            color: var(--primary);
            font-size: 0.85em;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            display: flex;
            justify-content: space-between;
        }
        
        .badge-a { color: var(--highlight-a); }
        .badge-b { color: var(--highlight-b); }

        /* Inputs */
        select, input[type="text"] {
            padding: 10px;
            border-radius: 6px;
            border: 1px solid #ddd;
            font-size: 14px;
            background: white;
            transition: all 0.2s;
            width: 100%;
            box-sizing: border-box; 
        }

        select:focus, input[type="text"]:focus {
            outline: none;
            border-color: var(--secondary);
            box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.2);
        }

        /* Slider */
        input[type=range] {
            -webkit-appearance: none;
            width: 100%;
            background: transparent;
            margin-top: 5px;
        }
        
        input[type=range]::-webkit-slider-thumb {
            -webkit-appearance: none;
            height: 16px;
            width: 16px;
            border-radius: 50%;
            background: var(--secondary);
            cursor: pointer;
            margin-top: -6px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.3);
        }
        
        input[type=range]::-webkit-slider-runnable-track {
            width: 100%;
            height: 4px;
            cursor: pointer;
            background: #ddd;
            border-radius: 2px;
        }

        /* Toggles */
        .toggle-wrapper {
            display: flex;
            align-items: center;
            gap: 10px;
            cursor: pointer;
            height: 40px; /* Align with inputs */
        }

        /* Plot Area */
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
            position: relative;
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

        .highlight-msg {
            position: absolute;
            top: 20px;
            right: 20px;
            padding: 8px 16px;
            background: #d4edda;
            color: #155724;
            border-radius: 4px;
            font-size: 0.9em;
            display: none;
            animation: fadeIn 0.3s;
            z-index: 10;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }

        @keyframes fadeIn { from { opacity: 0; } to { opacity: 1; } }

    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🔬 Metabolic Network Clustering</h1>
            <p class="subtitle">Interactive Clustered Correlation Analysis</p>
        </header>

        <div class="controls-area">
            <!-- Dataset Select -->
            <div class="control-group">
                <div class="control-label">Select Cohort</div>
                <select id="cohort-select" onchange="updateView()"></select>
            </div>

            <!-- Filter Slider -->
            <div class="control-group">
                <div class="control-label">
                    <span>Filter Weak Correlations</span>
                    <span id="slider-val">|r| < 0.0</span>
                </div>
                <input type="range" id="filter-slider" min="0" max="0.9" step="0.05" value="0" oninput="updateSliderVal(this.value)" onchange="updateView()">
            </div>

            <!-- Search A -->
            <div class="control-group">
                <div class="control-label">
                    <span>Search Metabolite A</span>
                    <span class="badge-a">●</span>
                </div>
                <input type="text" id="search-box-1" list="aa-list" placeholder="Highlight A..." onchange="handleSearch(1)">
            </div>

            <!-- Search B -->
            <div class="control-group">
                <div class="control-label">
                    <span>Search Metabolite B</span>
                    <span class="badge-b">●</span>
                </div>
                <input type="text" id="search-box-2" list="aa-list" placeholder="Highlight B..." onchange="handleSearch(2)">
                <datalist id="aa-list"></datalist>
            </div>

            <!-- Toggle Dendrogram -->
            <div class="control-group" style="justify-content: flex-end;">
                 <label class="toggle-wrapper">
                    <input type="checkbox" id="dendro-toggle" checked onchange="updateView()">
                    <span>Show Dendrograms</span>
                </label>
            </div>
        </div>

        <div class="plot-card">
            <div class="card-header" id="plot-title">Loading...</div>
            <div id="main-plot"></div>
            <div id="msg-box" class="highlight-msg"></div>
        </div>
        
    </div>

    <script>
        // Data injected from Python
        const DATA = __JSON_DATA__;
        
        // Configuration
        const COHORTS = Object.keys(DATA.cohorts).sort().reverse();
        
        const COLORSCALE = [
            [0.0, '#ffffd9'],
            [0.2, '#c7e9b4'],
            [0.4, '#41b6c4'],
            [0.6, '#1d91c0'],
            [0.8, '#225ea8'],
            [1.0, '#0c2c84']
        ];

        // State
        let currentHighlights = { 1: null, 2: null };

        function init() {
            const select = document.getElementById('cohort-select');
            
            // Populate dropdowns
            COHORTS.forEach(c => {
                select.add(new Option(c, c));
            });

            // Set Defaults
            if (COHORTS.includes('MS+MG')) select.value = 'MS+MG';

            // Populate Datalist ONCE (All cohorts share same metabolites)
            const firstCohort = COHORTS[0];
            const data = DATA.cohorts[firstCohort];
            if (data && data.x) {
                const list = document.getElementById('aa-list');
                list.innerHTML = '';
                data.x.forEach(aa => {
                    const opt = document.createElement('option');
                    opt.value = aa;
                    list.appendChild(opt);
                });
            }
            
            updateView();
        }

        function updateSliderVal(val) {
            document.getElementById('slider-val').textContent = `|r| < ${val}`;
        }

        function handleSearch(id) {
            const box = document.getElementById(`search-box-${id}`);
            const term = box.value.trim();
            
            if (!term) {
                if (currentHighlights[id]) {
                     currentHighlights[id] = null;
                     updateView();
                }
                return;
            }
            
            // Validate
            const data = DATA.cohorts[document.getElementById('cohort-select').value];
            if (data && data.x.includes(term)) {
                currentHighlights[id] = term;
                updateView();
            } else {
                if (currentHighlights[id] !== null) {
                    currentHighlights[id] = null;
                    updateView(); 
                }
                showMessage(`"${term}" not found in this cohort`, true);
            }
        }

        function showMessage(msg, isError=false) {
            const box = document.getElementById('msg-box');
            box.textContent = msg;
            box.style.background = isError ? '#f8d7da' : '#d4edda';
            box.style.color = isError ? '#721c24' : '#155724';
            box.style.display = 'block';
            setTimeout(() => { box.style.display = 'none'; }, 3000);
        }

        function updateView() {
            const cohortName = document.getElementById('cohort-select').value;
            const threshold = parseFloat(document.getElementById('filter-slider').value);
            const showDendro = document.getElementById('dendro-toggle').checked;

            renderHeatmap('main-plot', cohortName, threshold, showDendro);
            document.getElementById('plot-title').textContent = cohortName + " Clustered Matrix";
        }

        function renderHeatmap(divId, cohortName, threshold, showDendro) {
            const data = DATA.cohorts[cohortName];
            if (!data) return;

            // Apply Filter
            let z_filtered = data.z.map(row => row.slice());
            
            if (threshold > 0) {
                for(let i=0; i<z_filtered.length; i++) {
                    for(let j=0; j<z_filtered[i].length; j++) {
                        const val = z_filtered[i][j];
                        if (val !== null && Math.abs(val) < threshold) {
                            z_filtered[i][j] = null;
                        }
                    }
                }
            }

            const traces = [];
            const shapes = [];

            // Highlighting Logic
            [1, 2].forEach(id => {
                const term = currentHighlights[id];
                if (term) {
                    const idx = data.x.indexOf(term);
                    if (idx !== -1) {
                        const color = id === 1 ? '#ff00ff' : '#00ffff'; // Magenta vs Cyan
                        
                        // Highlight Row
                        shapes.push({
                            type: 'rect',
                            xref: 'paper', yref: 'y',
                            x0: 0, x1: 1,
                            y0: idx - 0.5, y1: idx + 0.5,
                            line: { color: color, width: 2 },
                            fillcolor: color,
                            opacity: 0.1
                        });
                        // Highlight Col
                        shapes.push({
                            type: 'rect',
                            xref: 'x', yref: 'paper',
                            x0: idx - 0.5, x1: idx + 0.5,
                            y0: 0, y1: 1,
                            line: { color: color, width: 2 },
                            fillcolor: color,
                            opacity: 0.1
                        });
                    }
                }
            });

            // Intersection Highlight?
            // If both are selected, the overlapping cell is naturally mixed color.
            // We could add a strong border to the intersection cell(s).
            if (currentHighlights[1] && currentHighlights[2]) {
                const idx1 = data.x.indexOf(currentHighlights[1]);
                const idx2 = data.x.indexOf(currentHighlights[2]);
                
                if (idx1 !== -1 && idx2 !== -1) {
                    // Cell (idx1, idx2) and (idx2, idx1)
                    [ {r: idx1, c: idx2}, {r: idx2, c: idx1} ].forEach(cell => {
                         shapes.push({
                            type: 'rect',
                            xref: 'x', yref: 'y',
                            x0: cell.c - 0.5, x1: cell.c + 0.5,
                            y0: cell.r - 0.5, y1: cell.r + 0.5,
                            line: { color: '#000000', width: 2 }, // Black border for intersection
                            fillcolor: 'rgba(0,0,0,0)'
                        });
                    });
                }
            }


            const layout = {
                height: 850,
                width: 1000,
                showlegend: false,
                margin: { t: 50, r: 50, b: 50, l: 50 },
                shapes: shapes,
                xaxis: { 
                    tickangle: -45, 
                    tickfont: { size: 10 },
                    domain: showDendro ? [0.15, 1] : [0, 1],
                    anchor: 'y'
                },
                yaxis: { 
                    automargin: true, 
                    autorange: 'reversed', 
                    tickfont: { size: 10 },
                    domain: showDendro ? [0, 0.85] : [0, 1],
                    anchor: 'x',
                    side: 'right'
                },
                xaxis2: { // Top Dendro
                    domain: showDendro ? [0.15, 1] : [0, 0], 
                    anchor: 'y2', showgrid: false, zeroline: false, showticklabels: false
                },
                yaxis2: { 
                    domain: showDendro ? [0.85, 1] : [1, 1], 
                    anchor: 'x2', showgrid: false, zeroline: false, showticklabels: false
                },
                xaxis3: { // Left Dendro
                    domain: showDendro ? [0, 0.15] : [0, 0], 
                    anchor: 'y3', showgrid: false, zeroline: false, showticklabels: false, autorange: 'reversed'
                },
                yaxis3: { 
                    domain: showDendro ? [0, 0.85] : [0, 0], 
                    anchor: 'x3', showgrid: false, zeroline: false, showticklabels: false, autorange: 'reversed'
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
                z: z_filtered,
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

            // 2. Dendrograms
            if (showDendro && data.dendrogram) {
                if (data.dendrogram.col) {
                    const d = data.dendrogram.col;
                    for (let i = 0; i < d.icoords.length; i++) {
                        traces.push({
                            x: d.icoords[i], y: d.dcoords[i],
                            mode: 'lines', type: 'scatter',
                            xaxis: 'x2', yaxis: 'y2',
                            line: { color: 'black', width: 1 },
                            hoverinfo: 'none'
                        });
                    }
                }
                if (data.dendrogram.row) {
                    const d = data.dendrogram.row;
                    for (let i = 0; i < d.icoords.length; i++) {
                        traces.push({
                            x: d.dcoords[i], y: d.icoords[i],
                            mode: 'lines', type: 'scatter',
                            xaxis: 'x3', yaxis: 'y3',
                            line: { color: 'black', width: 1 },
                            hoverinfo: 'none'
                        });
                    }
                }
            }

            Plotly.react(divId, traces, layout, config);
        }

        // Run
        document.addEventListener('DOMContentLoaded', init);

    </script>
</body>
</html>'''

    html = html.replace('__JSON_DATA__', json_content)
    return html
