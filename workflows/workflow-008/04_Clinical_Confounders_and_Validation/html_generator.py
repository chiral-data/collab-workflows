"""HTML Generator for Node 04: Clinical Confounders Dashboard
Generates interactive HTML visualization with Tabs, Filters, and Trendlines.
"""

import json
import os

def generate_confounders_html(json_filename='confounders_data.json'):
    """Generate HTML for Clinical Confounders and Validation visualization"""
    
    json_path = os.path.join('outputs', json_filename)
    if not os.path.exists(json_path):
        json_path = json_filename
        
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            json_content = f.read()
    except Exception as e:
        print(f"Warning: Could not read JSON file {json_path}: {e}")
        json_content = '{}'

    html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Clinical Confounders Dashboard</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {
            --primary: #2c3e50;
            --accent: #3498db;
            --bg: #f8f9fa;
            --card-bg: #ffffff;
            --text: #333;
            --control-color: #2ecc71; /* Green for Control */
            --ms-color: #e74c3c;      /* Red for MS */
        }
        
        body { font-family: 'Segoe UI', Roboto, sans-serif; margin: 0; padding: 0; background: var(--bg); color: var(--text); height: 100vh; display: flex; flex-direction: column; }
        
        /* Header */
        header { background: var(--card-bg); padding: 15px 30px; border-bottom: 1px solid #ddd; display: flex; justify-content: space-between; align-items: center; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }
        h1 { margin: 0; font-size: 1.2rem; color: var(--primary); }
        .timestamp { font-size: 0.8rem; color: #777; }

        /* Controls Bar */
        .controls { background: #fff; padding: 10px 30px; border-bottom: 1px solid #eee; display: flex; gap: 20px; align-items: center; }
        .tab-group { display: flex; gap: 5px; background: #eee; padding: 4px; border-radius: 6px; }
        .tab-btn { border: none; background: none; padding: 8px 16px; cursor: pointer; border-radius: 4px; font-weight: 500; color: #666; transition: 0.2s; }
        .tab-btn.active { background: white; color: var(--accent); box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
        
        .toggle-container { display: flex; align-items: center; gap: 8px; font-size: 0.9rem; margin-left: auto; }
        input[type="checkbox"] { accent-color: var(--accent); cursor: pointer; transform: scale(1.2); }

        /* Main Content Grid */
        #app { flex: 1; overflow-y: auto; padding: 20px 30px; }
        
        .grid-container {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
            gap: 20px;
            padding-bottom: 50px;
        }

        .chart-card {
            background: var(--card-bg);
            border-radius: 8px;
            padding: 10px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            transition: transform 0.2s;
            border: 1px solid transparent;
            height: 320px;
            display: flex;
            flex-direction: column;
        }
        
        .chart-card:hover { transform: translateY(-2px); box-shadow: 0 5px 15px rgba(0,0,0,0.1); }
        .chart-card.highlight { border-color: var(--accent); }
        
        .card-header { display: flex; justify-content: space-between; padding: 0 10px; margin-bottom: 5px; font-size: 0.9rem; font-weight: bold; color: var(--primary); }
        .stats-badge { font-size: 0.75rem; background: #eee; padding: 2px 6px; border-radius: 4px; font-weight: normal; }
        .sig-badge { background: #fee; color: #c0392b; font-weight: bold; }
        
        .plot-div { flex: 1; width: 100%; }

        .loading { text-align: center; padding: 50px; font-size: 1.2rem; color: #666; }
    </style>
    <script>
        const RAW_DATA = __JSON_DATA__;
    </script>
</head>
<body>
    <header>
        <h1>Clinical Confounders Analysis</h1>
        <span class="timestamp" id="timestamp">Generated: Just now</span>
    </header>

    <div class="controls">
        <div class="tab-group">
            <button class="tab-btn active" onclick="switchTab('fig2a')">Age (MS vs Ctrl)</button>
            <button class="tab-btn" onclick="switchTab('fig2b')">Duration (MS)</button>
            <button class="tab-btn" onclick="switchTab('fig2c')">EDSS (MS)</button>
        </div>

        <div class="toggle-container">
            <input type="checkbox" id="sigFilter" onchange="renderCurrentView()">
            <label for="sigFilter" title="Show only correlations with p < 0.05">Significant Only (p<0.05)</label>
        </div>
        
        <div class="toggle-container">
            <input type="checkbox" id="trendlineToggle" checked onchange="renderCurrentView()">
            <label for="trendlineToggle">Show Trendlines</label>
        </div>
    </div>

    <div id="app">
        <div class="loading">Initializing Dashboard...</div>
    </div>

    <script>
        let currentTab = 'fig2a';
        
        // Colors matching typical static plots
        const COLORS = {
            MS: '#e74c3c',      // Red
            Control: '#2ecc71', // Green
            MS_Line: '#c0392b',
            Ctrl_Line: '#27ae60'
        };

        function init() {
            document.getElementById('timestamp').textContent = new Date().toLocaleString();
            renderCurrentView();
        }

        function switchTab(tabKey) {
            currentTab = tabKey;
            // Update buttons
            document.querySelectorAll('.tab-btn').forEach(btn => {
                btn.classList.remove('active');
                if(btn.textContent.includes(getTabName(tabKey))) btn.classList.add('active');
            });
            renderCurrentView();
        }
        
        function getTabName(key) {
            if(key === 'fig2a') return 'Age';
            if(key === 'fig2b') return 'Duration';
            if(key === 'fig2c') return 'EDSS';
            return '';
        }

        function renderCurrentView() {
            const container = document.getElementById('app');
            const dataSection = RAW_DATA[currentTab];
            const traces = dataSection.traces;
            const isSigOnly = document.getElementById('sigFilter').checked;
            const showLines = document.getElementById('trendlineToggle').checked;
            
            // Clean container
            container.innerHTML = '<div class="grid-container" id="grid"></div>';
            const grid = document.getElementById('grid');

            // Process each AA
            traces.forEach((item, index) => {
                const aaName = item.aa;
                
                // Determine significance (if any p-value < 0.05)
                let msP = item.MS?.stats?.p_value ?? 1;
                let ctrlP = item.Control?.stats?.p_value ?? 1;
                let isSig = (msP < 0.05) || (ctrlP < 0.05 && currentTab === 'fig2a');

                if (isSigOnly && !isSig) return; // Skip if filter is on

                // Create Card
                const card = document.createElement('div');
                card.className = 'chart-card';
                
                // Header with P-values
                let pText = `P(MS): ${msP.toFixed(3)}`;
                if(currentTab === 'fig2a') pText += ` | P(Ctrl): ${ctrlP.toFixed(3)}`;
                
                card.innerHTML = `
                    <div class="card-header">
                        <span>${aaName}</span>
                        <span class="stats-badge ${isSig ? 'sig-badge' : ''}">${pText}</span>
                    </div>
                    <div class="plot-div" id="plot_${currentTab}_${index}"></div>
                `;
                grid.appendChild(card);

                // Prepare Plotly Data
                const plotData = [];
                
                // Helper to add traces
                const addGroup = (groupName, groupData, color, lineColor) => {
                    if(!groupData) return;
                    
                    // Scatter points
                    plotData.push({
                        x: groupData.x,
                        y: groupData.y,
                        mode: 'markers',
                        type: 'scatter',
                        name: groupName,
                        marker: { color: color, size: 5, opacity: 0.6 }
                    });

                    // Trendline
                    if(showLines && groupData.stats) {
                        plotData.push({
                            x: groupData.stats.line_x,
                            y: groupData.stats.line_y,
                            mode: 'lines',
                            type: 'scatter',
                            name: `${groupName} Fit`,
                            line: { color: lineColor, width: 2 },
                            showlegend: false,
                            hoverinfo: 'skip'
                        });
                    }
                };

                addGroup('MS', item.MS, COLORS.MS, COLORS.MS_Line);
                if(item.Control) addGroup('Control', item.Control, COLORS.Ctrl_Line);

                // Layout
                const layout = {
                    margin: { t: 10, r: 10, b: 30, l: 40 },
                    xaxis: { title: dataSection.variable, showgrid: true, gridcolor: '#eee' },
                    yaxis: { title: 'Level', showgrid: true, gridcolor: '#eee' },
                    paper_bgcolor: 'rgba(0,0,0,0)',
                    plot_bgcolor: 'rgba(0,0,0,0)',
                    showlegend: (index === 0), // Only show legend on first chart to save space
                    legend: { x: 0, y: 1, bgcolor: 'rgba(255,255,255,0.5)' }
                };

                Plotly.newPlot(`plot_${currentTab}_${index}`, plotData, layout, {responsive: true, displayModeBar: false});
            });
        }

        // Run
        setTimeout(init, 100);
    </script>
</body>
</html>"""
    
    # Inject JSON data
    html = html.replace('__JSON_DATA__', json_content)
    
    return html
